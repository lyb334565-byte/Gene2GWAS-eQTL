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
GWAS_VARIANT_API = "https://www.ebi.ac.uk/gwas/rest/api/v2/single-nucleotide-polymorphism"
GTEX_API = "https://gtexportal.org/api/v2"
ENSEMBL_API = "https://rest.ensembl.org"
USER_AGENT = "Gene2GWAS-eQTL/1.0 (academic research tool)"
GENCODE_BY_DATASET = {"gtex_v10": "v39", "gtex_v8": "v26"}
GENOME_BUILD = "GRCh38/hg38"

GWAS_COLUMNS = [
    "Gene", "Trait", "Reported_trait", "rsID", "rsID_Source", "Chromosome", "Position",
    "P_value", "Effect_type", "Effect_size", "Risk_allele", "Mapped_genes",
    "Study_accession", "PMID", "Association_ID", "Source_URL",
]
EQTL_COLUMNS = [
    "Gene", "GENCODE_ID", "Tissue", "rsID", "Variant_ID", "Chromosome",
    "Position", "P_value", "NES", "Dataset", "Source_URL",
]
OVERLAP_COLUMNS = [
    "Gene", "rsID", "GWAS_trait", "GWAS_P", "GTEx_tissue", "eQTL_P", "eQTL_NES",
    "GWAS_study", "GTEx_variant_ID", "Evidence_note",
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
    variant_keys = [
        "singleNucleotidePolymorphisms", "singleNucleotidePolymorphism",
        "single_nucleotide_polymorphisms", "single_nucleotide_polymorphism",
        "variants",
    ] if key in {"variants", "singleNucleotidePolymorphisms", "singleNucleotidePolymorphism"} else []
    candidate_keys = list(dict.fromkeys([key, *variant_keys, "data", "results", "items", "content"]))
    for parent in [payload.get("_embedded"), payload]:
        if not isinstance(parent, dict):
            continue
        for k in candidate_keys:
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


def _parse_location_object(obj: Any, depth: int = 0) -> tuple[str, int | None]:
    """Extract a chromosome/position pair from GWAS Catalog v2 or nested variant data."""
    if depth > 5:
        return "", None
    if isinstance(obj, str):
        match = re.search(r"(?:chr)?([0-9XYM]+)[:_-](\d+)", obj, re.I)
        if match:
            return match.group(1).removeprefix("chr"), int(match.group(2))
        return "", None
    if isinstance(obj, list):
        for item in obj:
            chrom, pos = _parse_location_object(item, depth + 1)
            if chrom and pos is not None:
                return chrom, pos
        return "", None
    if not isinstance(obj, dict):
        return "", None

    chrom = _get(obj, "chromosomeName", "chromosome_name", "chromosome", "chr", "chr_name")
    pos = _get(obj, "chromosomePosition", "chromosome_position", "position", "pos", "bp_location", "bpLocation")
    if chrom is not None and pos is not None and str(pos).isdigit():
        chrom_text = str(chrom).strip().removeprefix("chr")
        return chrom_text, int(pos)

    for key in ("locations", "location", "genomicLocations", "genomic_locations", "loci", "locus",
                "variant", "variants", "singleNucleotidePolymorphism", "single_nucleotide_polymorphism"):
        if key in obj:
            found_chrom, found_pos = _parse_location_object(obj[key], depth + 1)
            if found_chrom and found_pos is not None:
                return found_chrom, found_pos
    return "", None


def _location(row: dict) -> tuple[str, int | None]:
    """Read coordinates from direct fields or nested v2 locus/variant structures."""
    return _parse_location_object(row)


def _extract_rsids(value: Any) -> list[str]:
    """Extract explicit dbSNP rsIDs from a value while preserving their first-seen order."""
    found: list[str] = []

    def walk(item: Any) -> None:
        if isinstance(item, dict):
            for child in item.values():
                walk(child)
        elif isinstance(item, (list, tuple, set)):
            for child in item:
                walk(child)
        elif isinstance(item, str):
            found.extend(re.findall(r"\brs\d+\b", item, flags=re.I))

    walk(value)
    return list(dict.fromkeys(x.lower() for x in found))


def _record_rsids(row: dict) -> list[str]:
    """Inspect variant-specific fields in a GWAS association/variant record."""
    keys = (
        "rs_id", "rs_ids", "rsId", "rsID", "rsid", "snps", "snpIds", "snp_id", "snpId",
        "variant_id", "variantId", "variant", "variants", "loci", "locus",
        "singleNucleotidePolymorphism", "single_nucleotide_polymorphism",
        "risk_allele", "riskAllele", "strongest_allele", "strongestAllele", "id",
        "strongest_risk_alleles", "strongestRiskAlleles",
    )
    found: list[str] = []
    for key in keys:
        if key in row and row[key] is not None:
            found.extend(_extract_rsids(row[key]))
    return list(dict.fromkeys(found))


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
        rsids = _record_rsids(row)
        rsid = "; ".join(rsids)
        rsid_source = "GWAS Catalog association record" if rsid else "Unresolved"
        url_id = _get(row, "association_id", "associationId", "id")
        result.append({
            "Gene": queried_gene,
            "Trait": _text(_get(row, "efo_traits", "efoTraits", "efo_trait", "efoTrait", "trait", "mapped_traits")),
            "Reported_trait": _text(_get(row, "disease_trait", "diseaseTrait", "reported_trait")),
            "rsID": rsid, "rsID_Source": rsid_source,
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


def fetch_gwas_variants(gene_symbol: str, extended: bool = False, max_pages: int = 3,
                         page_size: int = 100, session: requests.Session | None = None) -> Pages:
    """Fetch GWAS Catalog variant records as a fallback for associations without an rsID."""
    own_session = session is None
    session = session or make_session()
    rows: list[dict] = []
    total: int | None = None
    truncated = False
    pages_fetched = 0
    try:
        for page in range(max_pages):
            params = {
                "mapped_gene": gene_symbol,
                "extended_geneset": str(extended).lower(),
                "page": page, "size": page_size,
            }
            payload = get_json(session, GWAS_VARIANT_API, params)
            page_rows = _collection(payload, "variants")
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


def _coordinate_key(chromosome: Any, position: Any) -> tuple[str, int] | None:
    chrom = str(chromosome or "").strip().lower().removeprefix("chr")
    if not chrom or position is None:
        return None
    try:
        pos = int(position)
    except (TypeError, ValueError, OverflowError):
        return None
    if pos <= 0:
        return None
    if chrom in {"m", "mt"}:
        chrom = "mt"
    return chrom, pos


def enrich_gwas_rsids(gwas: pd.DataFrame, variant_rows: list[dict]) -> pd.DataFrame:
    """Fill missing GWAS rsIDs only when GRCh38 coordinate mapping is unique.

    GWAS Catalog and GTEx use GRCh38 coordinates. Ambiguous coordinate-to-rsID
    mappings are deliberately left unresolved rather than guessing.
    """
    if gwas.empty or not variant_rows:
        return gwas.copy()
    coordinate_ids: dict[tuple[str, int], set[str]] = {}
    for row in variant_rows:
        chrom, pos = _location(row)
        key = _coordinate_key(chrom, pos)
        if key is None:
            continue
        rsids = _record_rsids(row)
        if not rsids:
            continue
        coordinate_ids.setdefault(key, set()).update(rsids)

    out = gwas.copy()
    if "rsID_Source" not in out.columns:
        out["rsID_Source"] = "Unresolved"
    for idx, row in out.iterrows():
        if _rs_ids(row.get("rsID", "")):
            if not str(row.get("rsID_Source", "")).strip() or row.get("rsID_Source") == "Unresolved":
                out.at[idx, "rsID_Source"] = "GWAS Catalog association record"
            continue
        key = _coordinate_key(row.get("Chromosome"), row.get("Position"))
        if key is None:
            continue
        candidates = coordinate_ids.get(key, set())
        if len(candidates) == 1:
            out.at[idx, "rsID"] = next(iter(candidates))
            out.at[idx, "rsID_Source"] = "GWAS Catalog variant coordinate match (GRCh38)"
        elif len(candidates) > 1:
            out.at[idx, "rsID_Source"] = "Ambiguous coordinate mapping; left unresolved"
    return out


def match_rsids(gwas: pd.DataFrame, eqtl: pd.DataFrame, max_matches: int = 20000) -> pd.DataFrame:
    """Exact rsID overlap only. This is NOT formal colocalization or allele harmonization."""
    if gwas.empty or eqtl.empty:
        return pd.DataFrame(columns=OVERLAP_COLUMNS)
    index: dict[str, list[dict]] = {}
    for r in eqtl.to_dict("records"):
        for snp in _rs_ids(r.get("rsID")):
            index.setdefault(snp, []).append(r)
    matches = []
    seen = set()
    for gr in gwas.to_dict("records"):
        for snp in _rs_ids(gr.get("rsID")):
            for er in index.get(snp, []):
                key = (gr.get("Association_ID"), gr.get("Study_accession"),
                       gr.get("Trait"), gr.get("P_value"), snp,
                       er.get("Tissue"), er.get("Variant_ID"))
                if key in seen:
                    continue
                seen.add(key)
                matches.append({
                    "Gene": gr.get("Gene"), "rsID": snp,
                    "GWAS_trait": gr.get("Trait") or gr.get("Reported_trait"),
                    "GWAS_P": gr.get("P_value"), "GTEx_tissue": er.get("Tissue"),
                    "eQTL_P": er.get("P_value"), "eQTL_NES": er.get("NES"),
                    "GWAS_study": gr.get("Study_accession"),
                    "GTEx_variant_ID": er.get("Variant_ID"),
                    "Evidence_note": "Exact rsID overlap only; NOT colocalization",
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
