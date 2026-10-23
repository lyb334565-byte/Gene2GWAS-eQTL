# Gene2GWAS-eQTL | Streamlit Community Cloud

A Streamlit browser app that queries human gene-associated GWAS Catalog v2 records and GTEx v2 significant cis-eQTL results, displays tables and plots, matches shared chromosome-position coordinates (GRCh38/hg38), and exports an Excel workbook.

> For academic exploratory use. A GWAS mapped gene is not a causal-gene assignment; chromosome-position overlap is **not** formal GWAS–eQTL colocalization and does not guarantee allele concordance. Not finding records is not evidence of no association.

## Deploy without installing Python on other computers

1. Log in to [GitHub](https://github.com/) and create a **public** repository (e.g. `Gene2GWAS-eQTL`).
2. Upload these **four extracted files** into the repository **root**, not within another folder:
   - `app.py` — Streamlit interface and Excel download
   - `core.py` — Ensembl / GWAS Catalog / GTEx API calls and result normalization
   - `requirements.txt` — Python dependencies
   - `README.md` — this file
3. Open [Streamlit Community Cloud](https://share.streamlit.io/), log in via GitHub and grant repository access.
4. Choose **Create app** → **Yup, I have an app**. Select repository, `main` branch and `app.py` main-file path.
5. Open **Advanced settings**, select Python **3.13** if offered (matches the tested local installation), save and select **Deploy**.
6. Once the app opens, test the `BACH1` gene and check the GWAS and GTEx status separately; use **下载与说明** to export Excel.
7. Use the deployed `https://...streamlit.app` URL from any authorized computer; no local Python installation is needed.

Official references:
- [Community Cloud deployment](https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/deploy)
- [Dependencies](https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/app-dependencies)
- [App sharing](https://docs.streamlit.io/deploy/streamlit-community-cloud/share-your-app)

## Notes

- Cloud deployment is a **separate Linux-hosted application**. Do not upload a Windows `.venv`, `.bat` files, passwords or research data you do not want made public.
- The application requires internet access to `www.ebi.ac.uk`, `gtexportal.org` and `rest.ensembl.org`; upstream APIs may change or be unavailable. API/network errors are displayed separately from empty results.
- This app currently queries only one human gene at a time. The GTEx result is limited to significant cis-eQTL and a configurable number of pages, so returned tables may be incomplete.
- Public deployments may receive requests from anyone; avoid treating this free research demo as a production-scale API harvesting service.

## Local usage (optional)

```shell
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

You need Python 3.10–3.13 locally. Windows MGLTools Python 2.5 should remain unchanged and cannot run this app.
