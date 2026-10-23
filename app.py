"""Run with: streamlit run app.py"""
from __future__ import annotations

import re

import pandas as pd
import plotly.express as px
import streamlit as st

from core import (
    APIError, EQTL_COLUMNS, GENOME_BUILD, GWAS_COLUMNS, OVERLAP_COLUMNS,
    enrich_gwas_rsids, fetch_ensembl, fetch_gtex_eqtl, fetch_gtex_gene,
    fetch_gwas, fetch_gwas_variants, gene_info_table, match_rsids, normalize_eqtl, normalize_gwas,
    utc_now, workbook_bytes,
)

st.set_page_config(
    page_title="Gene2GWAS-eQTL | 遗传学证据查询",
    page_icon="🧬", layout="wide", initial_sidebar_state="expanded",
)
st.markdown("""
<style>
.block-container {padding-top: 1.4rem; padding-bottom: 3rem;}
h1 {letter-spacing: -.035em;}
[data-testid="stMetric"] {background: #F3F6FA; border-radius: 12px; padding: .9rem 1rem;}
.small-note {font-size:.87rem; color:#586779;}
</style>
""", unsafe_allow_html=True)

CV_TISSUES = [
    "Artery_Aorta", "Artery_Coronary", "Artery_Tibial",
    "Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Whole_Blood",
]
TISSUES = CV_TISSUES + [
    "Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland",
    "Brain_Cortex", "Brain_Cerebellum", "Colon_Transverse", "Esophagus_Mucosa",
    "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Pancreas",
    "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg",
    "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Thyroid",
]


@st.cache_data(ttl=3600, show_spinner=False)
def run_query(query: str, dataset: str, tissues: tuple[str, ...],
              all_tissues: bool, extended: bool, max_gwas_pages: int,
              max_eqtl_pages: int) -> dict:
    result = {
        "queried_at": utc_now(), "query": query, "dataset": dataset,
        "tissues": "所有GTEx组织" if all_tissues else ", ".join(tissues),
        "gwas": pd.DataFrame(columns=GWAS_COLUMNS),
        "eqtl": pd.DataFrame(columns=EQTL_COLUMNS),
        "gene_info": pd.DataFrame(),
        "overlap": pd.DataFrame(columns=OVERLAP_COLUMNS),
        "gwas_pages": 0, "eqtl_pages": 0,
        "gwas_total": None, "eqtl_total": None,
        "gwas_truncated": False, "eqtl_truncated": False,
        "warnings": [], "errors": {}, "status": {},
    }
    ensembl = None
    gtex_gene = None
    try:
        ensembl = fetch_ensembl(query)
        result["status"]["Ensembl"] = "ok" if ensembl else "not_found"
    except APIError as exc:
        result["errors"]["Ensembl"] = str(exc)
        result["status"]["Ensembl"] = "error"

    try:
        gtex_gene = fetch_gtex_gene(query, dataset=dataset)
        result["status"]["GTEx_gene"] = "ok" if gtex_gene else "not_found"
    except (APIError, ValueError) as exc:
        result["errors"]["GTEx_gene"] = str(exc)
        result["status"]["GTEx_gene"] = "error"

    # GWAS lookup is independent of GTEx / Ensembl availability.
    gwas_gene = ((gtex_gene or {}).get("geneSymbol")
                 or (ensembl or {}).get("display_name") or query)
    result["gene_info"] = gene_info_table(query, gtex_gene, ensembl)
    try:
        g = fetch_gwas(gwas_gene, extended=extended, max_pages=max_gwas_pages)
        gwas_df = normalize_gwas(g.rows, gwas_gene)
        result["gwas_pages"] = g.pages_fetched
        result["gwas_total"] = g.total
        result["gwas_truncated"] = g.truncated
        result["status"]["GWAS"] = "ok" if len(g.rows) else "empty"

        # Some v2 association payloads omit rsID at the association level. Query
        # the Catalog's SNP resource and recover IDs only for unique GRCh38
        # coordinate matches; ambiguous mappings are not guessed.
        missing_rsid = ~gwas_df["rsID"].fillna("").astype(str).str.contains(r"\brs\d+\b", case=False, regex=True)
        coord_available = (gwas_df["Chromosome"].fillna("").astype(str).str.strip().ne("")
                           & pd.to_numeric(gwas_df["Position"], errors="coerce").notna())
        if (missing_rsid & coord_available).any():
            try:
                variants = fetch_gwas_variants(
                    gwas_gene, extended=extended, max_pages=max_gwas_pages,
                )
                before = int(gwas_df["rsID"].fillna("").astype(str).str.contains(r"\brs\d+\b", case=False, regex=True).sum())
                gwas_df = enrich_gwas_rsids(gwas_df, variants.rows)
                after = int(gwas_df["rsID"].fillna("").astype(str).str.contains(r"\brs\d+\b", case=False, regex=True).sum())
                result["status"]["GWAS_rsID_lookup"] = "ok" if after > before else "no_additional_matches"
                result["warnings"].append(f"rsID补充查询：坐标唯一匹配新增 {after - before} 条；未能安全映射的记录保留为空。")
            except APIError as exc:
                result["status"]["GWAS_rsID_lookup"] = "error"
                result["warnings"].append(f"GWAS rsID补充查询失败，保留原始记录：{exc}")
        result["gwas"] = gwas_df
    except APIError as exc:
        result["errors"]["GWAS"] = str(exc)
        result["status"]["GWAS"] = "error"

    if gtex_gene:
        try:
            e = fetch_gtex_eqtl(
                gtex_gene["gencodeId"],
                tissues=None if all_tissues else list(tissues), dataset=dataset,
                max_pages=max_eqtl_pages,
            )
            result["eqtl"] = normalize_eqtl(e.rows, gwas_gene)
            result["eqtl_pages"] = e.pages_fetched
            result["eqtl_total"] = e.total
            result["eqtl_truncated"] = e.truncated
            result["status"]["GTEx_eQTL"] = "ok" if len(e.rows) else "empty"
        except APIError as exc:
            result["errors"]["GTEx_eQTL"] = str(exc)
            result["status"]["GTEx_eQTL"] = "error"
    elif result["status"]["GTEx_gene"] == "not_found":
        result["status"]["GTEx_eQTL"] = "not_found"
    else:
        result["status"]["GTEx_eQTL"] = "error"

    result["overlap"] = match_rsids(result["gwas"], result["eqtl"])
    return result


def filter_p(df: pd.DataFrame, threshold: float | None) -> pd.DataFrame:
    if threshold is None or df.empty:
        return df
    return df[pd.to_numeric(df["P_value"], errors="coerce") <= threshold].copy()


def dataframe_view(df: pd.DataFrame, height: int = 420) -> None:
    if df.empty:
        st.info("该筛选条件下没有可展示的记录。请尝试放宽筛选条件。")
        return
    columns = {}
    if "Source_URL" in df.columns:
        columns["Source_URL"] = st.column_config.LinkColumn("来源链接", display_text="打开", help="打开数据提供方页面")
    st.dataframe(df, use_container_width=True, hide_index=True, height=height,
                 column_config=columns)


st.title("🧬 Gene2GWAS-eQTL")
st.caption("输入人类基因 → GWAS Catalog v2 关联证据 + GTEx v2 显著 eQTL → Excel 导出")
st.info("**重要说明：** GWAS 的基因映射不等于因果基因；相同 rsID 同时出现在 GWAS/eQTL 中，也不等于已证实共定位或相同效应等位基因。", icon="🔬")

with st.sidebar:
    st.header("查询设置")
    with st.form("search_form"):
        gene = st.text_input("基因名 / Ensembl Gene ID", value="BACH1",
                             placeholder="例如 BACH1 / NFE2L2 / ENSG...",
                             help="第一版每次查询一个人类基因；支持标准 HGNC gene symbol 或 Ensembl gene ID。")
        dataset = st.selectbox("GTEx 数据版本", ["gtex_v10", "gtex_v8"], index=0,
                               help="GTEx v10 对应 GENCODE v39；v8 对应 GENCODE v26。")
        all_tissues = st.checkbox("查询 GTEx 所有组织", value=False,
                                 help="可能返回大量数据；可提高 eQTL 最大页数。")
        tissues = st.multiselect("组织范围", TISSUES, default=CV_TISSUES,
                                 disabled=all_tissues,
                                 help="默认 6 个心血管相关组织。")
        extended = st.checkbox("GWAS 扩展基因集", value=False,
                               help="extended_geneset=true 会检索更广泛的相邻/映射基因集合。")
        st.markdown("**获取数量限制**")
        max_gwas_pages = st.slider("GWAS 最大页数（每页100条）", 1, 20, 3)
        max_eqtl_pages = st.slider("eQTL 最大页数（每页500条）", 1, 20, 3)
        submitted = st.form_submit_button("🔎 开始查询", type="primary", use_container_width=True)
    if st.button("🗑️ 清除缓存", use_container_width=True):
        run_query.clear()
        st.session_state.pop("result", None)
        st.success("已清除缓存。")
    st.divider()
    st.caption("数据源：GWAS Catalog、GTEx、Ensembl。采用顺序请求，不适合一次抓取成千上万个基因。")

if submitted:
    query = gene.strip()
    valid = bool(re.fullmatch(r"(?:[A-Za-z][A-Za-z0-9._-]{0,39}|ENSG\d+(?:\.\d+)?)", query))
    if not valid:
        st.error("请输入单个有效基因标识，例如 BACH1、NFE2L2 或 ENSG00000156273。")
    elif not all_tissues and not tissues:
        st.error("请至少选择一个 GTEx 组织，或勾选“查询所有组织”。")
    else:
        with st.spinner("正在顺序查询 Ensembl、GTEx 与 GWAS Catalog，请留意数据库访问速度…"):
            st.session_state["result"] = run_query(
                query, dataset, tuple(tissues), all_tissues, extended,
                max_gwas_pages, max_eqtl_pages,
            )

if "result" not in st.session_state:
    st.markdown("### 开始使用")
    st.write("左侧保留默认基因 **BACH1**，点击 **开始查询** 即可。也可以输入 `NFE2L2`、`NOS3`、`KLF2` 等人类基因。")
    a, b, c = st.columns(3)
    a.markdown("**GWAS**  \n疾病/性状、rsID、P值、效应")
    b.markdown("**GTEx eQTL**  \n组织、rsID、NES、P值")
    c.markdown("**Excel**  \n完整查询记录及数据源说明")
    st.stop()

r = st.session_state["result"]
st.markdown(f"### 查询结果：`{r['query']}`")
st.caption(f"查询时间：{r['queried_at']} · {r['dataset']} · 基因组坐标：{GENOME_BUILD}")
for source, message in r["errors"].items():
    st.error(f"{source} 查询失败（**不是阴性结果**）：{message}")
for message in r.get("warnings", []):
    st.warning(message)
if r["status"].get("GTEx_gene") == "not_found":
    st.warning("GTEx 所选数据版本中未能精确匹配该基因的版本化 GENCODE ID；eQTL 查询未执行。GWAS 仍独立查询。")
if r["status"].get("GWAS") == "empty":
    st.info("GWAS Catalog 当前基因映射条件下没有返回记录；这不等于该基因与疾病没有遗传学关联。")
if r["status"].get("GTEx_eQTL") == "empty":
    st.info("GTEx 在所选组织/版本下没有返回预计算的**显著**单组织 eQTL；不代表不存在调控关系。")
if r["gwas_truncated"]:
    st.warning(f"GWAS 共获取 {len(r['gwas']):,} 条，已达到分页上限；以下分析和导出仅基于**已获取的子集**。如需更多数据，请在左侧提高页数。")
if r["eqtl_truncated"]:
    st.warning(f"GTEx 共获取 {len(r['eqtl']):,} 条，已达到分页上限；以下分析和导出仅基于**已获取的子集**（可能不是最显著的前N条）。")

m1, m2, m3, m4 = st.columns(4)
m1.metric("GWAS 已获取记录", f"{len(r['gwas']):,}")
m2.metric("GTEx 已获取 eQTL", f"{len(r['eqtl']):,}")
m3.metric("eQTL 涉及组织", r["eqtl"]["Tissue"].nunique() if not r["eqtl"].empty else 0)
m4.metric("相同 rsID 匹配行", f"{len(r['overlap']):,}")

st.markdown("#### 筛选已获取的结果")
c1, c2 = st.columns(2)
with c1:
    gwas_limit = st.selectbox("GWAS P 值", ["全部", "≤ 5×10⁻⁸", "≤ 1×10⁻⁵", "≤ 0.05"], index=0)
with c2:
    tissue_values = sorted(x for x in r["eqtl"]["Tissue"].dropna().unique().tolist() if x)
    selected_tissues = st.multiselect("展示的 eQTL 组织（留空代表全部）", tissue_values, default=[])
threshold = {"全部": None, "≤ 5×10⁻⁸": 5e-8, "≤ 1×10⁻⁵": 1e-5, "≤ 0.05": 0.05}[gwas_limit]
gdf = filter_p(r["gwas"], threshold)
edf = r["eqtl"]
if selected_tissues:
    edf = edf[edf["Tissue"].isin(selected_tissues)].copy()
ov = match_rsids(gdf, edf)

tab1, tab2, tab3, tab4, tab5 = st.tabs(["GWAS 结果", "GTEx eQTL", "rsID 重叠", "基因信息", "下载与说明"])

with tab1:
    rsid_count = int(gdf["rsID"].fillna("").astype(str).str.contains(r"\brs\d+\b", case=False, regex=True).sum()) if "rsID" in gdf.columns else 0
    st.caption(f"显示 {len(gdf):,} / {len(r['gwas']):,} 条已获取 GWAS 关联；其中 {rsid_count:,} 条包含可识别 rsID。")
    dataframe_view(gdf)
    plotdf = gdf.copy()
    if not plotdf.empty:
        plotdf["Trait_label"] = plotdf["Trait"].where(plotdf["Trait"].astype(bool), plotdf["Reported_trait"])
        plotdf["Trait_label"] = plotdf["Trait_label"].replace("", "（未返回性状名称）")
        plotdf["P_value"] = pd.to_numeric(plotdf["P_value"], errors="coerce")
        plotdf = plotdf[plotdf["P_value"].notna() & (plotdf["P_value"] >= 0)].copy()
        if not plotdf.empty:
            plotdf["-log10(P)"] = -plotdf["P_value"].clip(lower=1e-300).map(__import__("math").log10)
            plotdf = (plotdf.groupby("Trait_label", as_index=False)["-log10(P)"].max()
                      .sort_values("-log10(P)", ascending=False).head(18))
            fig = px.bar(plotdf.sort_values("-log10(P)"), x="-log10(P)", y="Trait_label", orientation="h",
                         title="各性状最显著的 GWAS 关联（已获取结果）",
                         labels={"Trait_label": "疾病 / 表型"}, color_discrete_sequence=["#1F6E8C"])
            fig.update_layout(height=max(330, min(760, 32 * len(plotdf) + 130)), yaxis_title="")
            st.plotly_chart(fig, use_container_width=True)

with tab2:
    st.caption(f"显示 {len(edf):,} / {len(r['eqtl']):,} 条已获取的 GTEx **显著 cis-eQTL**；NES 方向与 GTEx 所定义的替代等位基因有关。")
    dataframe_view(edf)
    if not edf.empty:
        plotdf = edf.dropna(subset=["NES", "Tissue"]).copy()
        if not plotdf.empty:
            fig = px.box(plotdf, x="Tissue", y="NES", points=False,
                         title="不同组织的显著 eQTL NES 分布（已获取结果）",
                         color_discrete_sequence=["#1F6E8C"])
            fig.update_layout(height=480, xaxis_tickangle=-40)
            st.plotly_chart(fig, use_container_width=True)

with tab3:
    st.info("这里仅按 **完全相同的 rsID** 匹配两个结果表。程序会先读取关联记录中的 rsID；缺失时尝试通过 GWAS Catalog 变异记录按 GRCh38 染色体坐标补充，只有唯一匹配才会填入。仍无法唯一映射的记录会保留为空。该分析未进行 LD、等位基因协调、fine-mapping 或 coloc 共定位。")
    dataframe_view(ov)
    if ov.empty and (r["gwas_truncated"] or r["eqtl_truncated"]):
        st.warning("当前两组数据有截断，未发现 rsID 重叠不能解释为全量数据无重叠。")

with tab4:
    dataframe_view(r["gene_info"], height=330)
    st.caption("基因坐标优先来自与所选 GTEx 版本一致的 GENCODE 注释；Ensembl 作为补充来源。")

with tab5:
    st.markdown("**下载当前页面筛选后的结果（含元数据、已查询页数及截断状态）**")
    tables = {
        "Gene_Info": r["gene_info"],
        "GWAS_Results": gdf,
        "GTEx_eQTL": edf,
        "rsID_Overlap": ov,
        "Source_Status": pd.DataFrame([
            {"Source": k, "Status": v, "Error": r["errors"].get(k, "")}
            for k, v in r["status"].items()
        ]),
        "Sources": pd.DataFrame([
            {"Source": "GWAS Catalog v2", "URL": "https://www.ebi.ac.uk/gwas/docs/programmatic-access/rest-api/"},
            {"Source": "GTEx v2", "URL": "https://gtexportal.org/api/v2/docs"},
            {"Source": "Ensembl REST", "URL": "https://rest.ensembl.org/documentation/info/symbol_lookup"},
        ]),
    }
    meta = {
        "App": "Gene2GWAS-eQTL 1.0 (MVP)", "Gene": r["query"],
        "Queried_at": r["queried_at"], "GTEx_dataset": r["dataset"],
        "GTEx_tissues_requested": r["tissues"], "Genome_build": GENOME_BUILD,
        "GWAS_total_upstream": r["gwas_total"], "GWAS_pages_fetched": r["gwas_pages"],
        "GWAS_truncated": r["gwas_truncated"],
        "GTEx_total_upstream": r["eqtl_total"], "GTEx_pages_fetched": r["eqtl_pages"],
        "GTEx_truncated": r["eqtl_truncated"],
        "GWAS_P_filter_for_export": gwas_limit,
        "GTEx_tissue_filter_for_export": ", ".join(selected_tissues) or "no additional filter",
        "Scientific_caveat": "Gene mapping and rsID overlap do not establish causality or colocalization.",
    }
    excel = workbook_bytes(tables, meta)
    st.download_button("⬇️ 下载 Excel 工作簿", data=excel,
                       file_name=f"Gene2GWAS_eQTL_{re.sub('[^A-Za-z0-9_-]', '_', r['query'])}_{r['dataset']}.xlsx",
                       mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                       use_container_width=True, type="primary")
    st.markdown("""
**数据边界**
- GWAS Catalog 为人工整理的关联记录，不是某性状所有 SNP 的全基因组 summary statistics。
- GWAS 按基因映射查询可能包括附近位点，不能直接推断作用基因。
- GTEx `/singleTissueEqtl` 返回预计算的显著 eQTL，并不包含全部未显著检测结果。
- 超过页数上限时，工作簿只包括已获取并通过当前筛选的记录。保留页数和截断标记便于复现。
- 重叠表只是同 rsID 匹配，不提供因果结论。若需要正式共定位，应另取匹配的全区域 summary statistics。
""")
    st.markdown("[GWAS Catalog 接口文档](https://www.ebi.ac.uk/gwas/docs/programmatic-access/rest-api/) · [GTEx v2 接口文档](https://gtexportal.org/api/v2/docs) · [Ensembl REST 文档](https://rest.ensembl.org/documentation/info/symbol_lookup)")
