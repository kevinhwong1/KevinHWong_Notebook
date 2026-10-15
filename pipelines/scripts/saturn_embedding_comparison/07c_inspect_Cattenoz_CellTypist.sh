#!/bin/bash
#BSUB -J inspect_Dmel_Hsap
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 4
#BSUB -W 4:00
#BSUB -R "rusage[mem=96000]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/inspect_Dmel_Hsap.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/inspect_Dmel_Hsap.%J.err
#BSUB -u kxw755@earth.miami.edu
#BSUB -N
# =============================================================================
# 07c_inspect_Cattenoz_CellTypist.sh  --  DIAGNOSTICS ONLY, writes nothing but reports.
# cat-then-python: writes the Python script, then runs it.
# NOTE: the Python below contains no pipe characters (Pegasus heredoc issue).
# Run after 07b finishes.
# Submit:  bsub < /scratch/dark_genes/SATURN_Mnemi/scripts/07c_inspect_Cattenoz_CellTypist.sh
# =============================================================================
set -eo pipefail

BASE="/scratch/dark_genes/SATURN_Mnemi"
OUT="${BASE}/raw_data/inspection_Dmel_Hsap"
mkdir -p "${OUT}"

# Batch jobs don't read ~/.bashrc: activate by full path, as in 05a/05i/05o.
# Falls back to micromamba in case saturn_protein was rebuilt there after the Oct 2026 env wipe.
if [ -f ~/anaconda3/bin/activate ] && ~/anaconda3/bin/conda env list | grep -q "saturn_protein"; then
    source ~/anaconda3/bin/activate
    conda activate saturn_protein
else
    export MAMBA_ROOT_PREFIX=/nethome/kxw755/micromamba
    eval "$(/nethome/kxw755/.local/bin/micromamba shell hook --shell bash)"
    micromamba activate saturn_protein
fi
echo "python: $(which python)"
python -c "import scanpy, anndata; print('scanpy', scanpy.__version__, 'anndata', anndata.__version__)"

cat > "${OUT}/inspect_Dmel_Hsap.py" << 'EOF'
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad

BASE = Path("/scratch/dark_genes/SATURN_Mnemi")
RAW = BASE / "raw_data"
OUT = RAW / "inspection_Dmel_Hsap"
pd.set_option("display.width", 250)
pd.set_option("display.max_rows", 400)
pd.set_option("display.max_columns", 50)

def hdr(t):
    print("\n" + "=" * 90 + "\n" + t + "\n" + "=" * 90, flush=True)

def is_integer_counts(M, n=200000):
    if sp.issparse(M):
        vals = M.data[:n]
    else:
        vals = np.asarray(M).ravel()[:n]
    vals = np.asarray(vals)
    if vals.size == 0:
        return "EMPTY"
    frac_int = float(np.mean(np.isclose(vals, np.round(vals))))
    return "frac_integer=%.4f  min=%.3f  max=%.1f" % (frac_int, vals.min(), vals.max())

# ============================================================================
# FLY: Cattenoz 2020 NI (RRCZ22)
# ============================================================================
hdr("FLY  Cattenoz 2020 NI  (E-MTAB-8698, RRCZ22)")
fdir = RAW / "Dmel_Cattenoz_EMTAB8698" / "RRCZ22_filtered"
mtx_dirs = sorted(set(p.parent for p in fdir.rglob("matrix.mtx*")))
print("matrix dirs found:", [str(p) for p in mtx_dirs])
mdir = mtx_dirs[0]

a = sc.read_10x_mtx(mdir, var_names="gene_ids", cache=False)
print("shape (cells x genes):", a.shape)
print("var columns:", list(a.var.columns))
print(a.var.head(5).to_string())
print("var_names starting with FBgn: %d / %d" % (sum(v.startswith("FBgn") for v in a.var_names), a.n_vars))
if "feature_types" in a.var:
    print("feature_types:", a.var["feature_types"].value_counts().to_dict())
print("barcode examples:", list(a.obs_names[:3]))
print("X integer check:", is_integer_counts(a.X))

tot = np.asarray(a.X.sum(1)).ravel()
ngenes = np.asarray((a.X > 0).sum(1)).ravel()
sym = a.var["gene_symbols"].astype(str) if "gene_symbols" in a.var else pd.Series(a.var_names, index=a.var_names)
mito = sym.str.startswith("mt:").values
pct_mt = 100 * np.asarray(a.X[:, mito].sum(1)).ravel() / np.maximum(tot, 1)
print("n mito genes (symbol starts mt:):", int(mito.sum()))
qc = pd.DataFrame({"total_counts": tot, "n_genes": ngenes, "pct_mt": pct_mt})
print(qc.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95]).round(1).to_string())

lab = pd.read_csv(RAW / "Dmel_Cattenoz_EMTAB8698" / "NI_cell_cluster_ID.txt", sep="\t")
print("\nlabel file:", lab.shape, list(lab.columns))
bc_core = pd.Index([b.split("-")[0] for b in a.obs_names])
lab_ids = pd.Index(lab["Cell_ID"].astype(str))
print("barcodes in matrix: %d   labelled barcodes: %d" % (len(bc_core), len(lab_ids)))
print("matrix barcodes with a label: %d" % bc_core.isin(lab_ids).sum())
print("labels without a matrix barcode: %d" % (~lab_ids.isin(bc_core)).sum())
print("duplicated core barcodes in matrix: %d" % bc_core.duplicated().sum())
print("\nauthor cluster sizes:")
print(lab["Cluster_ID"].value_counts().to_string())

qc["Cluster_ID"] = pd.Series(lab.set_index("Cell_ID")["Cluster_ID"]).reindex(bc_core).values
qc["Cluster_ID"] = qc["Cluster_ID"].fillna("UNLABELLED")
print("\nQC of unlabelled vs labelled cells:")
print(qc.groupby(qc["Cluster_ID"].eq("UNLABELLED"))[["total_counts", "n_genes", "pct_mt"]].median().round(1).to_string())
qc.to_csv(OUT / "Dmel_RRCZ22_cell_qc.csv")

# ============================================================================
# HUMAN: CellTypist Myeloid compartment
# ============================================================================
def describe_celltypist(path, tag):
    hdr("HUMAN  CellTypist " + tag + "  " + path.name)
    h = ad.read_h5ad(path, backed="r")
    print("shape:", h.shape)
    print("obs columns:", list(h.obs.columns))
    print("has raw:", h.raw is not None)
    if h.raw is not None:
        print("raw shape:", h.raw.shape)
    print("var index examples:", list(h.var_names[:5]))
    print("var columns:", list(h.var.columns))
    obs = h.obs.copy()
    label_cols = [c for c in obs.columns if ("celltype" in c.lower()) or ("cell_type" in c.lower())]
    print("candidate label columns:", label_cols)
    for c in label_cols:
        print("  %s : %d levels" % (c, obs[c].nunique()))
    donor_col = "donor_id" if "donor_id" in obs else None
    tissue_col = "tissue" if "tissue" in obs else None
    print("\nassay counts:")
    print(obs["assay"].value_counts().to_string())
    if donor_col:
        print("\ndonor x assay:")
        print(pd.crosstab(obs[donor_col], obs["assay"]).to_string())
    if tissue_col:
        print("\ntissue x assay:")
        print(pd.crosstab(obs[tissue_col], obs["assay"]).to_string())
    return h, obs, label_cols, donor_col, tissue_col

def raw_integer_check(h, n_cells=2000):
    idx = np.sort(np.random.default_rng(0).choice(h.n_obs, size=min(n_cells, h.n_obs), replace=False))
    if h.raw is not None:
        M = h.raw.X[idx]
        print("raw.X integer check (%d cells): %s" % (len(idx), is_integer_counts(M)))
    M = h.X[idx]
    print("X integer check (%d cells): %s" % (len(idx), is_integer_counts(M)))

mye_path = RAW / "Hsap_CellTypist_CXG" / "CellTypist_Myeloid.h5ad"
h, obs, label_cols, donor_col, tissue_col = describe_celltypist(mye_path, "MYELOID")
raw_integer_check(h)
v3 = obs[obs["assay"].astype(str).str.contains("3' v3")]
print("\n>>> MYELOID cells on 10x 3' v3: %d (donors: %s)" % (len(v3), sorted(v3[donor_col].astype(str).unique()) if donor_col else "NA"))
for c in label_cols:
    print("\n--- 3' v3 myeloid: %s counts" % c)
    print(v3[c].value_counts().to_string())
best = [c for c in label_cols if "manual" in c.lower() or "curated" in c.lower()]
lc = best[0] if best else label_cols[0]
print("\nusing label column for tissue breakdown:", lc)
if tissue_col:
    tab = pd.crosstab(v3[lc], v3[tissue_col])
    print(tab.to_string())
    tab.to_csv(OUT / "Hsap_myeloid_3v3_celltype_x_tissue.csv")
obs.to_csv(OUT / "Hsap_myeloid_obs.csv")
h.file.close()

# ============================================================================
# HUMAN: CellTypist Global (obs only) -> lymphocyte slice candidates
# ============================================================================
gl_path = RAW / "Hsap_CellTypist_CXG" / "CellTypist_Global.h5ad"
g, gobs, glabel_cols, gdonor, gtissue = describe_celltypist(gl_path, "GLOBAL")
gv3 = gobs[gobs["assay"].astype(str).str.contains("3' v3")]
glc_best = [c for c in glabel_cols if "manual" in c.lower() or "curated" in c.lower()]
glc = glc_best[0] if glc_best else glabel_cols[0]
print("\n>>> GLOBAL cells on 10x 3' v3: %d" % len(gv3))
print("\n--- 3' v3 global: %s counts" % glc)
print(gv3[glc].value_counts().to_string())
in_mye = gv3.index.isin(obs.index)
print("\n3' v3 global cells also in Myeloid file: %d" % in_mye.sum())
print("3' v3 global cells NOT in Myeloid file (lymphoid + other): %d" % (~in_mye).sum())
print("\n--- 3' v3 non-myeloid cells: %s counts" % glc)
print(gv3.loc[~in_mye, glc].value_counts().to_string())
if gtissue:
    t2 = pd.crosstab(gv3.loc[~in_mye, glc], gv3.loc[~in_mye, gtissue])
    t2.to_csv(OUT / "Hsap_nonmyeloid_3v3_celltype_x_tissue.csv")
shared = obs.index.intersection(gobs.index)
if len(shared) and lc in obs and glc in gobs:
    agree = (obs.loc[shared, lc].astype(str) == gobs.loc[shared, glc].astype(str)).mean()
    print("\nlabel agreement Myeloid-file vs Global-file on shared cells (%s vs %s): %.3f" % (lc, glc, agree))
raw_integer_check(g)
g.file.close()

hdr("DONE  tables written to " + str(OUT))
EOF

python "${OUT}/inspect_Dmel_Hsap.py" > "${OUT}/inspect_Dmel_Hsap_report.txt" 2>&1
echo "Report: ${OUT}/inspect_Dmel_Hsap_report.txt"
