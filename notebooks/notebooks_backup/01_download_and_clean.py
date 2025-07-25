# %% [markdown]
# # SynTact — Data-prep pipeline  ❱  TCGA-OV  +  GTEx (v8)
# 
# Cleans RNA-seq TPM matrices, maps Ensembl → HGNC symbols, merges tumour and
# normal expression, adds a plasma-membrane flag (HPA), and saves a feature
# matrix ready for target-ranking.
# 
# **Directory assumptions** (notebooks folder):
# 
# ```
# data/
# ├─ raw/
# │  ├─ tcga_ov_expression.tsv.gz            # STAR-TPM matrix from UCSC Xena
# │  ├─ GTEx_gene_tpm_v10.gct.gz                 # GTEx v8 gene TPM matrix
# │  ├─ gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz
# │  └─ hpa_subcellular_location.tsv         # HPA v24 sub-cell locations
# └─ processed/                              # outputs land here
# ```
# 
# Run this file as a notebook in VS Code / Cursor (`# %%` cells) **or** execute
# with `python 01_download_and_clean.py`.


# %%
import gzip
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd

RAW = Path("../data/raw")
PROCESSED = Path("../data/processed")
PROCESSED.mkdir(exist_ok=True, parents=True)

# --------------------------------------------------------------------------------------
# utility: build / load gene-id ➜ HGNC symbol map  (cached for speed)
# --------------------------------------------------------------------------------------
CACHE = PROCESSED / "gene_map.pkl"
if CACHE.exists():
    gene_map: dict[str, str] = pickle.loads(CACHE.read_bytes())
else:
    gtf_path = RAW / "gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz"
    gene_map: dict[str, str] = {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if parts[2] != "gene":
                continue
            attr = parts[8]
            ens_id = re.search(r'gene_id "([^"]+)"', attr).group(1).split(".")[0]
            symbol = re.search(r'gene_name "([^"]+)"', attr).group(1)
            gene_type = re.search(r'gene_type "([^"]+)"', attr).group(1)
            if gene_type == "protein_coding":
                gene_map[ens_id] = symbol
    CACHE.write_bytes(pickle.dumps(gene_map))
print(f"🗺️  gene_map entries: {len(gene_map):,}")

# --------------------------------------------------------------------------------------
# 1️⃣  TCGA-OV  — load STAR-TPM, map → symbols, save
# --------------------------------------------------------------------------------------
print("\n▶ TCGA-OV …")
if not (PROCESSED / "tcga_ov_expression_cleaned.csv").exists():
    tcga_path = RAW / "tcga_ov_expression.tsv.gz"
    tcga = pd.read_csv(tcga_path, sep="\t", index_col=0, compression="gzip")
    tcga.index = tcga.index.str.split(".").str[0]
    tcga = tcga[tcga.index.isin(gene_map)]
    tcga.index = tcga.index.map(gene_map)
    tcga = tcga[~tcga.index.duplicated(keep="first")]
    tcga.to_csv(PROCESSED / "tcga_ov_expression_cleaned.csv")
    print("   saved cleaned TCGA :", tcga.shape)
else:
    tcga = pd.read_csv(PROCESSED / "tcga_ov_expression_cleaned.csv", index_col=0)
    print("   loaded cached TCGA :", tcga.shape)

# --------------------------------------------------------------------------------------
# 2️⃣  GTEx v8  — load .gct.gz, map → symbols, save (with cache guard)
# --------------------------------------------------------------------------------------
print("\n▶ GTEx v8 …")

gtex_clean_path = PROCESSED / "gtex_expression_cleaned.csv"
gtex = None

# ——— Try loading the cache —————————————
if gtex_clean_path.exists():
    gtex_tmp = pd.read_csv(gtex_clean_path, index_col=0)
    # If it has far fewer rows than expected, drop and rebuild
    if gtex_tmp.shape[0] < 10000:
        print("   ⚠️  Cached GTEx too small (", gtex_tmp.shape, ") — rebuilding …")
        gtex_clean_path.unlink()
    else:
        gtex = gtex_tmp
        print("   loaded cached GTEx :", gtex.shape)

# ——— Build cache if needed —————————————————
if gtex is None:
    gct = RAW / "GTEx_gene_tpm_v10.gct.gz"
    gtex = pd.read_csv(gct, sep="\t", skiprows=2, header=0, index_col=0)
    gtex = gtex.drop(columns=["Description"])
    gtex.index = gtex.index.str.split(".").str[0]
    gtex = gtex[gtex.index.isin(gene_map)]
    gtex.index = gtex.index.map(gene_map)
    gtex = gtex[~gtex.index.duplicated(keep="first")]
    gtex.to_csv(gtex_clean_path)
    print("   saved cleaned GTEx :", gtex.shape)


# --------------------------------------------------------------------------------------
# 3️⃣  Surfaceome list — union of HPA + local CSPA Excel  (no manual overrides)
# --------------------------------------------------------------------------------------
print("\n▶ Building union surfaceome …")

# --- a) HPA plasma-membrane annotations ----------------------------------
hpa = pd.read_csv(RAW / "hpa_subcellular_location.tsv", sep="\t")
loc_cols = ["Main location", "Enhanced", "Supported", "Approved", "Uncertain"]

hpa_mask = hpa[loc_cols].apply(
    lambda row: row.astype(str).str.contains("Plasma membrane", na=False).any(),
    axis=1,
)
hpa_genes = set(hpa.loc[hpa_mask, "Gene name"].str.upper())
print(f"   HPA surface genes          : {len(hpa_genes):,}")

# --- b) CSPA Excel  ------------------------------------------------------
cspa_xlsx = RAW / "S2_File.xlsx"
if not cspa_xlsx.exists():
    raise FileNotFoundError(
        "❌  CSPA Excel not found in data/raw/. "
        "Place S2_File.xlsx there first."
    )

# First sheet = human surfaceome table
cspa_df = pd.read_excel(cspa_xlsx, sheet_name=0)

# Auto-detect a gene-symbol column
candidate_cols = [
    c for c in cspa_df.columns
    # requires the S2_File.xlsx to be downloaded from the CSPA site, should be column H
    if ("entrez gene symbol" in c.lower())
]
if not candidate_cols:
    raise ValueError("Could not locate a gene-symbol column in CSPA sheet")

primary_col = next(col for col in candidate_cols if cspa_df[col].notna().sum() > 0)
print("   Using CSPA gene column     :", primary_col)

cspa_genes = set(cspa_df[primary_col].dropna().astype(str).str.upper())
print(f"   CSPA surface genes         : {len(cspa_genes):,}")

# Known surfaceome additions for ovarian cancer specificly
curated_additions = {
    "CLDN6", "CLDN3", "CLDN4", "CLDN16", "CLDN18",
    "AMHR2", "FSHR",
    "TACSTD2", "L1CAM", "PTK7",
    "CEACAM6", "CEACAM5",
    "ITGB4", "ITGA6",
    "ICAM1", "NCAM1", "CXCR4",
}

# --- c) Union set --------------------------------------------------------
surface_set: set[str] = hpa_genes | cspa_genes
print(f"   Union surfaceome           : {len(surface_set):,}")

# --- d) Union set with curated additions know for ovarian cancer ---------
surface_set |= curated_additions
print(f"   + curated additions (+{len(curated_additions)}) → {len(surface_set):,} total surface genes")


# --------------------------------------------------------------------------------------
# 4️⃣  Merge | tumour vs normal | compute log2FC | flag surface genes
# --------------------------------------------------------------------------------------
print("\n▶ Merge & compute log2FC …")

# gene-wise means (raw TPM)
tcga_mean = tcga.mean(axis=1)
gtex_mean = gtex.mean(axis=1)

common = tcga_mean.index.intersection(gtex_mean.index)
print("   common protein-coding genes:", len(common))

merged = pd.DataFrame({
    "TCGA_OV_TPM": tcga_mean.loc[common],
    "GTEx_TPM": gtex_mean.loc[common],
})
merged["log2FC"] = (
    np.log2(merged["TCGA_OV_TPM"] + 1) - np.log2(merged["GTEx_TPM"] + 1)
)
merged["is_surface"] = merged.index.str.upper().isin(surface_set)

out_path = PROCESSED / "features_matrix.csv"
merged.to_csv(out_path)
print("   feature matrix saved   :", merged.shape)
print("   surface genes matched  :", merged["is_surface"].sum())

# %% [markdown]
# ## ✅ Sanity check
# * Expect **MUC16, FOLR1, MSLN, EPCAM, CD276 (B7-H3)** among the top‐ranked surface genes.
# * `GTEx_TPM` is ≥ 0 and `log2FC` is positive for tumour-specific antigens.

# %%
print(merged[merged["is_surface"]].sort_values("log2FC", ascending=False).head(10))
# %%
# Look up known CAR-T targets in the final matrix
known_targets = ["MUC16", "FOLR1", "MSLN", "EPCAM", "CD276"]
print(merged.loc[merged.index.intersection(known_targets)].sort_values("log2FC", ascending=False))

# %% [markdown]
## 📖 How to Read `features_matrix.csv`

# At a glance, each row in `features_matrix.csv` is a gene and the key columns mean:
#<br>
#(TPM = transcripts per million)
#<br>
#- **TCGA_OV_TPM**  (We want this high)<br>
#  – *What?* Average expression (in TPM) of that gene across all TCGA ovarian cancer samples (n ≈ 429).  <br>
#  – *Why?* Shows how “on” the gene is, in tumours.<br>
# <br>
#- **GTEx_TPM**  (We want this low)<br>
#  – *What?* Average expression (in TPM) of that gene across all GTEx normal‐tissue samples (n ≈ 17 382).  <br>
#  – *Why?* Tells you how common the gene is in healthy human tissues.<br>
# <br>
#- **log₂FC**  (We want this positive)<br>
#  – *What?* `log₂( TCGA_OV_TPM + 1 ) – log₂( GTEx_TPM + 1 )`  <br>
#  – *Why?*  
#  - **Positive** (↑) → higher in ovarian tumours vs. average normal  <br>
#  - **Negative** (↓) → higher in normal tissues  <br>
#  - The number is in “folds.” E.g.  
#    - log₂FC = 1 → 2× higher in tumour  
#    - log₂FC = –1 → 2× higher in normal

#- **is_surface**  (We want this True)<br>
#  – *What?* `True` if the gene is annotated as a cell‐surface protein (HPA ∪ CSPA).  <br>
#  – *Why?* CAR-T therapies can only target proteins accessible on the cell surface.

# ---

# ### Tips for Interpretation

# 1. **High TCGA_OV_TPM + low GTEx_TPM + log₂FC ≫ 0**  
#    → *Tumour-specific surface antigen*  
# 2. **Negative log₂FC**  
#    → The gene is more abundant in healthy tissues than in the average tumour sample (common for general epithelial or immune markers).  
# 3. **Tissue-filtered fold‐change**  
#    If you want “ovary-specific” comparisons, you can recompute GTEx_TPM using only GTEx “Ovary” (and related) samples—this yields a more focused `log₂FC_ovary`.  

# Use these fields to rank and prioritize your top CAR-T target candidates!


# %%
# --- 5️⃣  Filter candidates ——————————————————————————————————————————————
candidates = (
    merged[merged.is_surface]
      .query("TCGA_OV_TPM > 2 and GTEx_TPM < 1 and log2FC > .5")
      .sort_values("log2FC", ascending=False)
)
candidates.head(20)

# %%
