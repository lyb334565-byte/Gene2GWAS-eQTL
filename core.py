"""Gene2GWAS-eQTL: API access, evidence normalization and export (MVP).

No live network calls are performed when this module is imported.
GWAS Catalog data are CURATED associations mapped to/near genes, not causal proof.
GTEx results are the precomputed significant single-tissue eQTL subset.
"""
from __future__ import annotations

import io
import json
import re
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any

import pandas as pd
import requests
from openpyxl.styles import PatternFill, Font
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

GWAS_API = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"
GTEX_API = "https://gtexportal.org/api/v2"
ENSEMBL_API = "https://rest.ensembl.org"
USER_AGENT = "Gene2GWAS-eQTL/1.0 (academic research tool)"
GENCODE_BY_DATASET = {"gtex_v10": "v39", "gtex_v8": "v26"}
GENOME_BUILD = "GRCh38/hg38"

GWAS_COLUMNS = [
    "Gene", "Trait", "Reported_trait", "rsID", "Chromosome", "Position",
    "P_value", "Effect_type", "Effect_size", "Risk_allele", "Mapped_genes",
    "Study_accession", "PMID", "Association_ID", "Source_URL",
]
EQTL_COLUMNS = [
    "Gene", "GENCODE_ID", "Tissue", "rsID", "Variant_ID", "Chromosome",
    "Position", "P_value", "NES", "Dataset", "Source_URL",
]
OVERLAP_COLUMNS = [
    "Gene", "Chromosome", "Position", "GWAS_trait", "GWAS_P", "GTEx_tissue",
    "eQTL_P", "eQTL_NES", "GWAS_rsID", "GTEx_rsID", "GWAS_study",
    "GTEx_variant_ID", "Evidence_note",
]


class APIError(RuntimeError):
    """Network, HTTP, or unexpected upstream API payload error."""


@dataclass
class Pages:
    rows: list[dict[str, Any]]
    total: int | None
    truncated: bool
    pages_fetched: int


def make_session() -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=3, backoff_factor=0.7,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["GET"]),
        respect_retry_after_header=True,
    )
    session.mount("https://", HTTPAdapter(max_retries=retry))
    session.headers.update({"User-Agent": USER_AGENT, "Accept": "application/json"})
    return session


def get_json(session: requests.Session, url: str, params: dict | None = None) -> Any:
    try:
        response = session.get(url, params=params, timeout=(8, 35))
        response.raise_for_status()
        return response.json()
    except requests.HTTPError as exc:
        status = exc.response.status_code if exc.response is not None else "unknown"
        raise APIError(f"接口返回 HTTP {status}：{url}") from exc
    except requests.RequestException as exc:
        raise APIError(f"网络连接失败：{url}；{str(exc)[:180]}") from exc
    except ValueError as exc:
        raise APIError(f"接口未返回合法 JSON：{url}") from exc


def _collection(payload: Any, key: str) -> list[dict]:
    """Handle GWAS v2 HAL or plain content as well as GTEx 'data'.

    Do not silently treat an unknown response format as 'zero results'.
    """
    if isinstance(payload, list):
        return [r for r in payload if isinstance(r, dict)]
    if not isinstance(payload, dict):
        raise APIError("接口响应结构无法识别：不是 JSON 对象/数组")
    for parent in [payload.get("_embedded"), payload]:
        if not isinstance(parent, dict):
            continue
        for k in [key, "data", "results", "items", "content"]:
            if k in parent and isinstance(parent[k], list):
                return [r for r in parent[k] if isinstance(r, dict)]
    page = payload.get("page", payload.get("paging_info", {}))
    if isinstance(page, dict) and any(page.get(k) == 0 for k in ("totalElements", "totalNumberOfItems", "totalPages", "numberOfPages")):
        return []
    raise APIError("接口返回了未识别的数据结构；请检查 API 版本或项目 README")


def _total(payload: dict, is_gtex: bool) -> int | None:
    page = payload.get("paging_info" if is_gtex else "page", {})
    if not isinstance(page, dict):
        return None
    for key in ("totalNumberOfItems", "totalElements", "total_items"):
        try:
            if page.get(key) is not None:
                return int(page[key])
        except (TypeError, ValueError):
            pass
    return None


def _has_next(payload: dict, page_index: int, page_size: int, n_received: int, is_gtex: bool) -> bool:
    meta = payload.get("paging_info" if is_gtex else "page", {})
    if isinstance(meta, dict):
        for key in ("numberOfPages", "totalPages"):
            if meta.get(key) is not None:
                return page_index + 1 < int(meta[key])
        total = _total(payload, is_gtex)
        if total is not None:
            return (page_index + 1) * page_size < total
    links = payload.get("_links", {})
    if isinstance(links, dict) and "next" in links:
        return bool(links["next"])
    return n_received >= page_size


def fetch_gwas(gene_symbol: str, extended: bool = False, max_pages: int = 3,
               page_size: int = 100, session: requests.Session | None = None) -> Pages:
    own_session = session is None
    session = session or make_session()
    rows: list[dict] = []
    total: int | None = None
    truncated = False
    pages_fetched = 0
    try:
        for page in range(max_pages):
            params = {
                "mapped_gene": gene_symbol, "extended_geneset": str(extended).lower(),
                "sort": "p_value", "direction": "asc", "page": page, "size": page_size,
            }
            payload = get_json(session, GWAS_API, params)
            page_rows = _collection(payload, "associations")
            rows.extend(page_rows)
            pages_fetched += 1
            if total is None:
                total = _total(payload, is_gtex=False)
            has_next = _has_next(payload, page, page_size, len(page_rows), is_gtex=False)
            if not has_next:
                break
            if page + 1 == max_pages:
                truncated = True
            time.sleep(0.12)
        return Pages(rows, total, truncated, pages_fetched)
    finally:
        if own_session:
            session.close()


def _versioned_gene_match(row: dict, query: str) -> bool:
    expected = query.strip().upper()
    symbol = str(row.get("geneSymbol", "")).upper()
    gencode_id = str(row.get("gencodeId", "")).upper()
    return expected == symbol or expected == gencode_id or expected == gencode_id.split(".")[0]


def fetch_gtex_gene(query: str, dataset: str = "gtex_v10",
                    session: requests.Session | None = None) -> dict | None:
    """Resolve against the GENCODE annotation VERSION used by the GTEx dataset.

    Ensembl IDs without a GENCODE suffix cannot be used directly for eQTL query.
    """
    if dataset not in GENCODE_BY_DATASET:
        raise ValueError("不支持的 GTEx 数据版本")
    own_session = session is None
    session = session or make_session()
    params = {
        "geneId": query, "gencodeVersion": GENCODE_BY_DATASET[dataset],
        "genomeBuild": GENOME_BUILD, "itemsPerPage": 250, "page": 0,
    }
    try:
        payload = get_json(session, f"{GTEX_API}/reference/gene", params)
        matches = [x for x in _collection(payload, "data") if _versioned_gene_match(x, query)]
        if not matches:
            payload = get_json(session, f"{GTEX_API}/reference/geneSearch", params)
            matches = [x for x in _collection(payload, "data") if _versioned_gene_match(x, query)]
        if not matches:
            return None
        # Prefer canonical exact-symbol match and a complete versioned GENCODE ID.
        matches.sort(key=lambda r: (str(r.get("geneSymbol", "")).upper() != query.upper(),
                                    "." not in str(r.get("gencodeId", ""))))
        chosen = dict(matches[0])
        if not re.match(r"^ENSG\d+\.\d+$", str(chosen.get("gencodeId", ""))):
            raise APIError("GTEx 返回的 GENCODE ID 缺少版本号；无法安全查询 eQTL")
        return chosen
    finally:
        if own_session:
            session.close()


def fetch_ensembl(query: str, session: requests.Session | None = None) -> dict | None:
    own_session = session is None
    session = session or make_session()
    clean = query.strip()
    url = (f"{ENSEMBL_API}/lookup/id/{clean.split('.')[0]}" if re.fullmatch(r"ENSG\d+(?:\.\d+)?", clean, re.I)
           else f"{ENSEMBL_API}/lookup/symbol/homo_sapiens/{clean}")
    try:
        try:
            data = get_json(session, url, {"content-type": "application/json"})
        except APIError as exc:
            if "HTTP 400" in str(exc) or "HTTP 404" in str(exc):
                return None
            raise
        if not isinstance(data, dict) or not data.get("id"):
            return None
        return data
    finally:
        if own_session:
            session.close()


def fetch_gtex_eqtl(gencode_id: str, tissues: list[str] | None, dataset: str = "gtex_v10",
                    max_pages: int = 3, page_size: int = 500,
                    session: requests.Session | None = None) -> Pages:
    own_session = session is None
    session = session or make_session()
    rows: list[dict] = []
    total: int | None = None
    truncated = False
    pages_fetched = 0
    try:
        for page in range(max_pages):
            params = {
                "gencodeId": gencode_id, "datasetId": dataset,
                "itemsPerPage": page_size, "page": page,
            }
            if tissues:  # FastAPI supports repeated query params for an array of strings
                params["tissueSiteDetailId"] = tissues
            payload = get_json(session, f"{GTEX_API}/association/singleTissueEqtl", params)
            page_rows = _collection(payload, "data")
            rows.extend(page_rows)
            pages_fetched += 1
            if total is None:
                total = _total(payload, is_gtex=True)
            has_next = _has_next(payload, page, page_size, len(page_rows), is_gtex=True)
            if not has_next:
                break
            if page + 1 == max_pages:
                truncated = True
            # GTEx specifically requests sequential queries for low throughput use.
            time.sleep(0.2)
        return Pages(rows, total, truncated, pages_fetched)
    finally:
        if own_session:
            session.close()


def _get(row: dict, *keys: str) -> Any:
    for key in keys:
        val = row.get(key)
        if val is not None and val != "":
            return val
    return None


def _text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, dict):
        for key in ("name", "efo_trait", "disease_trait", "shortForm", "label", "trait", "accessionId", "rs_id", "rsId", "value"):
            if value.get(key) is not None:
                return _text(value[key])
        return json.dumps(value, ensure_ascii=False, default=str)
    if isinstance(value, (list, tuple)):
        return "; ".join(_text(x) for x in value if x is not None)
    return str(value)


def _number(value: Any) -> float | None:
    if value is None or isinstance(value, (dict, list)):
        return None
    try:
        result = float(str(value).strip().replace(",", ""))
        if pd.isna(result):
            return None
        return result
    except (TypeError, ValueError, OverflowError):
        return None


def _location(row: dict) -> tuple[str, int | None]:
    loc = _get(row, "locations", "location", "genomicLocations")
    if isinstance(loc, list) and loc:
        loc = loc[0]
    if isinstance(loc, dict):
        ch = _get(loc, "chromosomeName", "chromosome", "chr")
        pos = _get(loc, "chromosomePosition", "position", "pos")
        return str(ch or "").removeprefix("chr"), int(pos) if str(pos or "").isdigit() else None
    if isinstance(loc, str):
        m = re.search(r"(?:chr)?([0-9XYM]+):(\d+)", loc, re.I)
        if m:
            return m.group(1), int(m.group(2))
    chrom = _get(row, "chromosome", "chr")
    pos = _get(row, "position", "pos")
    return str(chrom or "").removeprefix("chr"), int(pos) if str(pos or "").isdigit() else None


def normalize_gwas(rows: list[dict], queried_gene: str) -> pd.DataFrame:
    result = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        chrom, pos = _location(row)
        study = row.get("study") if isinstance(row.get("study"), dict) else {}
        pub = row.get("publication") if isinstance(row.get("publication"), dict) else {}
        pubinfo = study.get("publicationInfo", {}) if isinstance(study.get("publicationInfo"), dict) else {}
        pval = _number(_get(row, "p_value", "pvalue", "pValue"))
        if pval is None:
            mant = _number(_get(row, "pvalueMantissa", "p_value_mantissa"))
            exponent = _number(_get(row, "pvalueExponent", "p_value_exponent"))
            if mant is not None and exponent is not None:
                pval = mant * 10 ** exponent
        beta = _get(row, "beta_num", "betaNum", "beta")
        odds = _get(row, "or_per_copy_num", "orPerCopyNum", "odds_ratio", "or_value")
        rsid = _text(_get(row, "rs_id", "rs_ids", "rsId", "rsID", "snps", "snpIds", "variant_id"))
        url_id = _get(row, "association_id", "associationId", "id")
        result.append({
            "Gene": queried_gene,
            "Trait": _text(_get(row, "efo_traits", "efoTraits", "efo_trait", "efoTrait", "trait", "mapped_traits")),
            "Reported_trait": _text(_get(row, "disease_trait", "diseaseTrait", "reported_trait")),
            "rsID": rsid,
            "Chromosome": chrom, "Position": pos,
            "P_value": pval,
            "Effect_type": "Beta" if beta is not None else ("OR" if odds is not None else ""),
            "Effect_size": _number(beta if beta is not None else odds),
            "Risk_allele": _text(_get(row, "snp_effect_allele", "strongestAllele", "strongest_allele", "risk_allele")),
            "Mapped_genes": _text(_get(row, "mapped_genes", "mappedGenes", "mapped_gene", "reported_genes")),
            "Study_accession": _text(_get(row, "accession_id", "study_accession", "studyId") or _get(study, "accessionId", "accession_id")),
            "PMID": _text(_get(row, "pubmed_id", "pubmedId") or _get(pub, "pubmedId", "pubmed_id") or _get(pubinfo, "pubmedId", "pubmed_id")),
            "Association_ID": _text(url_id),
            "Source_URL": f"https://www.ebi.ac.uk/gwas/rest/api/v2/associations/{url_id}" if url_id is not None else "",
        })
    return pd.DataFrame(result, columns=GWAS_COLUMNS)


def normalize_eqtl(rows: list[dict], queried_gene: str) -> pd.DataFrame:
    result = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        result.append({
            "Gene": _text(_get(row, "geneSymbol", "geneSymbolUpper")) or queried_gene,
            "GENCODE_ID": _text(row.get("gencodeId")),
            "Tissue": _text(row.get("tissueSiteDetailId")),
            "rsID": _text(_get(row, "snpId", "snpIdUpper")),
            "Variant_ID": _text(row.get("variantId")),
            "Chromosome": _text(row.get("chromosome")),
            "Position": row.get("pos"),
            "P_value": _number(row.get("pValue")),
            "NES": _number(row.get("nes")),
            "Dataset": _text(row.get("datasetId")),
            "Source_URL": "https://gtexportal.org/home/gene/" + queried_gene,
        })
    return pd.DataFrame(result, columns=EQTL_COLUMNS)


def _rs_ids(value: Any) -> set[str]:
    return {x.lower() for x in re.findall(r"\brs\d+\b", str(value), flags=re.I)}


def _chromosome_key(value: Any) -> str:
    """Normalize chromosome labels (e.g. chr1/1, chrX/X) for coordinate matching."""
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    chrom = re.sub(r"^chr", "", str(value).strip(), flags=re.I).upper()
    if not chrom:
        return ""
    if chrom in {"M", "MT", "MITO", "MITOCHONDRIAL"}:
        return "MT"
    if chrom in {"X", "Y", "XY"}:
        return chrom
    if chrom.isdigit():
        return str(int(chrom))
    return chrom


def _position_key(value: Any) -> int | None:
    """Normalize integer genomic coordinates; reject missing/non-integer positions."""
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    text = str(value).strip().replace(",", "")
    if not text:
        return None
    try:
        number = float(text)
        if not (number > 0):
            return None
        position = int(number)
        if position != number:
            return None
        return position
    except (TypeError, ValueError, OverflowError):
        return None


def _coordinate_key(chromosome: Any, position: Any) -> tuple[str, int] | None:
    chrom = _chromosome_key(chromosome)
    pos = _position_key(position)
    return (chrom, pos) if chrom and pos is not None else None


def match_positions(gwas: pd.DataFrame, eqtl: pd.DataFrame, max_matches: int = 20000) -> pd.DataFrame:
    """Match exact chromosome + position pairs. This is not colocalization or allele harmonization."""
    if gwas.empty or eqtl.empty:
        return pd.DataFrame(columns=OVERLAP_COLUMNS)

    index: dict[tuple[str, int], list[dict]] = {}
    for row in eqtl.to_dict("records"):
        coordinate = _coordinate_key(row.get("Chromosome"), row.get("Position"))
        if coordinate is not None:
            index.setdefault(coordinate, []).append(row)

    matches = []
    seen = set()
    for gr in gwas.to_dict("records"):
        coordinate = _coordinate_key(gr.get("Chromosome"), gr.get("Position"))
        if coordinate is None:
            continue
        chromosome, position = coordinate
        for er in index.get(coordinate, []):
            key = (
                gr.get("Association_ID"), gr.get("Study_accession"),
                gr.get("Trait"), gr.get("Reported_trait"), gr.get("P_value"),
                chromosome, position, er.get("Tissue"), er.get("Variant_ID"),
                er.get("P_value"),
            )
            if key in seen:
                continue
            seen.add(key)
            matches.append({
                "Gene": gr.get("Gene") or er.get("Gene"),
                "Chromosome": chromosome,
                "Position": position,
                "GWAS_trait": gr.get("Trait") or gr.get("Reported_trait"),
                "GWAS_P": gr.get("P_value"),
                "GTEx_tissue": er.get("Tissue"),
                "eQTL_P": er.get("P_value"),
                "eQTL_NES": er.get("NES"),
                "GWAS_rsID": gr.get("rsID"),
                "GTEx_rsID": er.get("rsID"),
                "GWAS_study": gr.get("Study_accession"),
                "GTEx_variant_ID": er.get("Variant_ID"),
                "Evidence_note": "Exact chromosome-position overlap; assumes GRCh38/hg38; NOT colocalization or allele harmonization",
            })
            if len(matches) >= max_matches:
                return pd.DataFrame(matches, columns=OVERLAP_COLUMNS)
    return pd.DataFrame(matches, columns=OVERLAP_COLUMNS)

def gene_info_table(query: str, gtex_gene: dict | None, ensembl: dict | None) -> pd.DataFrame:
    result = [{"Field": "Query", "Value": query}]
    fields = (
        ("Gene symbol", (gtex_gene or {}).get("geneSymbol") or (ensembl or {}).get("display_name")),
        ("GTEx GENCODE ID (versioned)", (gtex_gene or {}).get("gencodeId")),
        ("Ensembl Gene ID", (ensembl or {}).get("id") or str((gtex_gene or {}).get("gencodeId", "")).split(".")[0]),
        ("Chromosome", (gtex_gene or {}).get("chromosome") or (ensembl or {}).get("seq_region_name")),
        ("Start (GRCh38)", (gtex_gene or {}).get("start") or (ensembl or {}).get("start")),
        ("End (GRCh38)", (gtex_gene or {}).get("end") or (ensembl or {}).get("end")),
        ("Annotation", (gtex_gene or {}).get("description") or (ensembl or {}).get("description")),
    )
    result.extend({"Field": key, "Value": _text(val)} for key, val in fields)
    return pd.DataFrame(result)


def workbook_bytes(tables: dict[str, pd.DataFrame], metadata: dict[str, Any]) -> bytes:
    out = io.BytesIO()
    metadata_df = pd.DataFrame([{"Parameter": k, "Value": _text(v)} for k, v in metadata.items()])
    with pd.ExcelWriter(out, engine="openpyxl") as writer:
        metadata_df.to_excel(writer, sheet_name="Query_Info", index=False)
        for name, frame in tables.items():
            frame.to_excel(writer, sheet_name=name[:31], index=False)
        for sheet in writer.book:
            sheet.freeze_panes = "A2"
            sheet.auto_filter.ref = sheet.dimensions
            for cell in sheet[1]:
                cell.font = Font(name="Calibri", size=11, bold=True, color="FFFFFF")
                cell.fill = PatternFill("solid", fgColor="16324F")
            for col in sheet.columns:
                letter = col[0].column_letter
                max_length = max((len(str(c.value or "")) for c in list(col)[:300]), default=10)
                sheet.column_dimensions[letter].width = min(max(max_length + 2, 12), 42)
    return out.getvalue()


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
