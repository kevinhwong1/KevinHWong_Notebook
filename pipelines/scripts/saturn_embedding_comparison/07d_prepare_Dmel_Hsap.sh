#!/usr/bin/env bash
# ============================================================
# 07d_prepare_Dmel_Hsap.sh
#
# Builds SATURN-ready h5ads for the two new species, same format as
# 01_prepare_h5ads.sh (raw counts in .X as CSR; obs cell_type + species):
#
#   Dmel  Cattenoz 2020 NI (E-MTAB-8698, RRCZ22), 10x 3'
#         - keep the 7,517 cells with author labels (drop 89 unlabelled, low-gene cells)
#         - cell_type = author Cluster_ID (PL-*, CC, LM-1, LM-2)
#         - var_names = FlyBase FBgn IDs
#
#   Hsap  CellTypist cross-tissue atlas (CELLxGENE), 10x 3' v3 only (donors D496, D503)
#         - counts from raw.X; var_names = Ensembl ENSG IDs
#         - cell_type = Manually_curated_celltype
#         - composition (decided 2026-10-15; revised same day to mirror the zebrafish cell types
#           in the earlier Drer/Crob/Mlei SATURN run: progenitors and erythroid cells added, and the
#           lymphocyte slice is now a matched T/B/NK counterpart rather than a negative control):
#             keep all : Nonclassical monocytes, Intestinal macrophages, Intermediate
#                        macrophages, DC1, DC2, migDC, Erythrophagocytic macrophages, Cycling
#             subsample: Classical monocytes 2,400 (600 each from blood, spleen, bone marrow, lung)
#                        Alveolar macrophages 2,400
#                        Mast cells 500 (from Global file)
#                        Lymphocyte slice 2,000 = 7 groups x ~286 (CD4 T, CD8 T, gd T, MAIT,
#                        NK, B, plasma), proportional to fine labels within each group
#                        (counterparts of zebrafish T and B cells)
#                        Progenitor: all 3' v3 cells, capped at 1,000 (EXTRA_TYPES; zebrafish HPC /
#                        myeloid progenitor counterpart)
#                        Erythroid: all 3' v3 cells, capped at 1,000 (EXTRA_TYPES; zebrafish erythroid
#                        counterpart; also a check on erythroid x Ciona cLRP-1 mixing)
#             dropped  : all doublets, Megakaryocytes, pDC, ILC3, GC/Pro-B/Pre-B, ABCs,
#                        T_CD4/CD8, Cycling T&NK, activated CD4
#             Changing the cells does not change the gene list (var = all genes in raw.X), so
#             the ESM embeddings from 02_embed.sh stay valid. The new cells are sampled after all
#             earlier picks, so the original 10,541 cells are unchanged (same seed).
#
# obs_names are prefixed with the species code (Dmel_, Hsap_) so barcodes can never
# collide across species. Extra obs columns (tissue, donor_id, sample_group) are kept
# for downstream interpretation; SATURN only reads cell_type.
#
# Writes: 02_processed_h5ad/Dmel_saturn.h5ad, 02_processed_h5ad/Hsap_saturn.h5ad
#         02_processed_h5ad/Dmel_Hsap_prep_composition.csv
#
# Submit: bsub < /scratch/dark_genes/SATURN_Mnemi/scripts/07d_prepare_Dmel_Hsap.sh
# ============================================================
#BSUB -J prep_Dmel_Hsap
#BSUB -q normal
#BSUB -P dark_genes
#BSUB -n 1
#BSUB -W 02:00
#BSUB -R "rusage[mem=32000] span[hosts=1]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/prep_Dmel_Hsap.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/prep_Dmel_Hsap.%J.err
#BSUB -u kxw755@earth.miami.edu
#BSUB -N

set -eo pipefail

BASE="/scratch/dark_genes/SATURN_Mnemi"
OUT_DIR="${BASE}/02_processed_h5ad"
mkdir -p "${OUT_DIR}"

# Batch jobs don't read ~/.bashrc; same activation as 05a (relax -e/-u for conda scripts)
set +eu
source ~/anaconda3/bin/activate
conda activate saturn_protein
set -eo pipefail
echo "python: $(which python)"

echo "=== 07d_prepare_Dmel_Hsap.sh ==="; date; hostname

cat > "${OUT_DIR}/07d_prepare_Dmel_Hsap.py" << 'EOF'
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

SEED = 42
BASE = Path("/scratch/dark_genes/SATURN_Mnemi")
RAW = BASE / "raw_data"
OUT = BASE / "02_processed_h5ad"
rng = np.random.default_rng(SEED)

def is_integer(mat, n=100000):
    data = mat.data[:n] if sp.issparse(mat) else np.asarray(mat).ravel()[:n]
    data = data[data > 0]
    return len(data) > 0 and float(np.mean(np.abs(data - np.round(data)) < 1e-3)) > 0.95

def write_saturn(a, path):
    a.X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    a.X = a.X.astype(np.float32)
    a.write_h5ad(path, compression="gzip")
    print("  wrote %s  (%d cells x %d genes, %d cell_type labels)" % (path.name, a.n_obs, a.n_vars, a.obs["cell_type"].nunique()))

def sample_idx(index, n):
    index = pd.Index(index)
    if n >= len(index):
        return list(index)
    return list(rng.choice(index, size=n, replace=False))

# ============================================================
# DMEL  Cattenoz 2020 NI
# ============================================================
print("\n[Dmel] loading Cattenoz RRCZ22 filtered matrix")
mdir = sorted(set(p.parent for p in (RAW / "Dmel_Cattenoz_EMTAB8698" / "RRCZ22_filtered").rglob("matrix.mtx*")))[0]
fl = sc.read_10x_mtx(mdir, var_names="gene_ids", cache=False)
assert is_integer(fl.X), "Dmel X not integer"
lab = pd.read_csv(RAW / "Dmel_Cattenoz_EMTAB8698" / "NI_cell_cluster_ID.txt", sep="\t")
lab = lab.set_index("Cell_ID")["Cluster_ID"].astype(str)
core = pd.Index([b.split("-")[0] for b in fl.obs_names])
keep = core.isin(lab.index)
print("  cells: %d total, %d labelled kept, %d unlabelled dropped" % (fl.n_obs, keep.sum(), (~keep).sum()))
fl = fl[keep].copy()
core = core[keep]
fl.obs = pd.DataFrame({
    "cell_type": lab.reindex(core).values,
    "species": "Dmel",
    "tissue": "larval hemolymph",
    "donor_id": "OregonR_pool20F",
    "sample_group": "Cattenoz2020_NI",
}, index=["Dmel_" + b for b in fl.obs_names])
fl.var = pd.DataFrame({"gene_symbol": fl.var["gene_symbols"].astype(str).values}, index=fl.var_names)
write_saturn(fl, OUT / "Dmel_saturn.h5ad")

# ============================================================
# HSAP  CellTypist, 10x 3' v3 only
# ============================================================
LC = "Manually_curated_celltype"
V3 = "10x 3' v3"
KEEP_ALL = ["Nonclassical monocytes", "Intestinal macrophages", "Intermediate macrophages",
            "DC1", "DC2", "migDC", "Erythrophagocytic macrophages", "Cycling"]
CLASSICAL_PER_TISSUE = 600
CLASSICAL_TISSUES = ["blood", "spleen", "bone marrow", "lung"]
N_ALVEOLAR = 2400
N_MAST = 500
LYMPH_GROUPS = {
    "CD4 T":  ["Tnaive/CM_CD4", "Teffector/EM_CD4", "Trm_Th1/Th17", "Tfh", "Tregs"],
    "CD8 T":  ["Tnaive/CM_CD8", "Tem/emra_CD8", "Trm_gut_CD8", "Trm/em_CD8"],
    "gd T":   ["Trm_Tgd", "Tgd_CRTAM+"],
    "MAIT":   ["MAIT"],
    "NK":     ["NK_CD16+", "NK_CD56bright_CD16-"],
    "B":      ["Naive B cells", "Memory B cells"],
    "Plasma": ["Plasma cells", "Plasmablasts"],
}
N_LYMPH_TOTAL = 2000
# Added 2026-10-15: counterparts of zebrafish HPC / myeloid progenitor and erythroid cells.
# Taken from whichever file holds them (Myeloid or Global non-myeloid), 3' v3 only.
EXTRA_TYPES = {"Progenitor": 1000, "Erythroid": 1000}   # label -> cap (all cells if fewer)

print("\n[Hsap] reading obs tables")
mye = ad.read_h5ad(RAW / "Hsap_CellTypist_CXG" / "CellTypist_Myeloid.h5ad", backed="r")
glo = ad.read_h5ad(RAW / "Hsap_CellTypist_CXG" / "CellTypist_Global.h5ad", backed="r")
mo = mye.obs[mye.obs["assay"].astype(str) == V3]
go = glo.obs[glo.obs["assay"].astype(str) == V3]
go_nonmye = go[~go.index.isin(mye.obs.index)]
print("  3' v3 myeloid: %d   3' v3 non-myeloid: %d" % (len(mo), len(go_nonmye)))

picked = []   # (obs_name, sample_group)
for ct in KEEP_ALL:
    ids = mo.index[mo[LC].astype(str) == ct]
    picked += [(i, "myeloid_keepall") for i in ids]
cm = mo[mo[LC].astype(str) == "Classical monocytes"]
for t in CLASSICAL_TISSUES:
    ids = cm.index[cm["tissue"].astype(str) == t]
    picked += [(i, "myeloid_classical_" + t.replace(" ", "_")) for i in sample_idx(ids, CLASSICAL_PER_TISSUE)]
am = mo.index[mo[LC].astype(str) == "Alveolar macrophages"]
picked += [(i, "myeloid_alveolar") for i in sample_idx(am, N_ALVEOLAR)]
mye_ids = [p[0] for p in picked]

glo_picked = []
mast = go_nonmye.index[go_nonmye[LC].astype(str) == "Mast cells"]
glo_picked += [(i, "mast") for i in sample_idx(mast, N_MAST)]
per_group = int(round(N_LYMPH_TOTAL / len(LYMPH_GROUPS)))
for g, labels in LYMPH_GROUPS.items():
    sub = go_nonmye[go_nonmye[LC].astype(str).isin(labels)]
    glo_picked += [(i, "lymph_" + g.replace(" ", "_")) for i in sample_idx(sub.index, per_group)]
for lab_name, cap in EXTRA_TYPES.items():
    in_mye = mo.index[mo[LC].astype(str) == lab_name]
    in_glo = go_nonmye.index[go_nonmye[LC].astype(str) == lab_name]
    print("  %s on 3' v3: %d in Myeloid file, %d in Global (non-myeloid); cap %d"
          % (lab_name, len(in_mye), len(in_glo), cap))
    pool = [(i, "mye") for i in in_mye] + [(i, "glo") for i in in_glo]
    if len(pool) > cap:
        sel = rng.choice(len(pool), size=cap, replace=False)
        pool = [pool[j] for j in sorted(sel)]
    grp_name = lab_name.lower()
    picked += [(i, grp_name) for i, f in pool if f == "mye"]
    glo_picked += [(i, grp_name) for i, f in pool if f == "glo"]
mye_ids = [p[0] for p in picked]
glo_ids = [p[0] for p in glo_picked]
assert not set(mye_ids).intersection(glo_ids), "overlap between myeloid and global picks"

def extract(backed, ids):
    pos = np.sort(backed.obs_names.get_indexer(ids))
    assert (pos >= 0).all()
    X = backed.raw.X[pos]
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    a = ad.AnnData(X=X, obs=backed.obs.iloc[pos].copy(), var=backed.raw.var.copy())
    return a

print("  extracting raw counts")
a1 = extract(mye, mye_ids)
a2 = extract(glo, glo_ids)
assert list(a1.var_names) == list(a2.var_names), "gene order differs between files"
hs = ad.concat([a1, a2], join="outer", merge="first")
hs.var = a1.var.copy()
assert is_integer(hs.X), "Hsap raw.X not integer"
grp = dict(picked + glo_picked)
hs.obs = pd.DataFrame({
    "cell_type": hs.obs[LC].astype(str).values,
    "species": "Hsap",
    "tissue": hs.obs["tissue"].astype(str).values,
    "donor_id": hs.obs["donor_id"].astype(str).values,
    "sample_group": [grp[i] for i in hs.obs_names],
}, index=["Hsap_" + i for i in hs.obs_names])
sym_col = "feature_name" if "feature_name" in hs.var else "gene_symbols"
hs.var = pd.DataFrame({"gene_symbol": hs.var[sym_col].astype(str).values,
                       "feature_biotype": hs.var["feature_biotype"].astype(str).values if "feature_biotype" in hs.var else "NA"},
                      index=hs.var_names)
write_saturn(hs, OUT / "Hsap_saturn.h5ad")
mye.file.close(); glo.file.close()

# ============================================================
# COMPOSITION SUMMARY
# ============================================================
print("\n=== composition ===")
rows = []
for a in [fl, hs]:
    t = a.obs.groupby(["species", "sample_group", "cell_type"]).size().rename("n_cells").reset_index()
    rows.append(t)
comp = pd.concat(rows)
comp = comp[comp["n_cells"] > 0]
comp.to_csv(OUT / "Dmel_Hsap_prep_composition.csv", index=False)
print(comp.to_string(index=False))
print("\nHsap tissue x sample_group:")
print(pd.crosstab(hs.obs["tissue"], hs.obs["sample_group"]).to_string())
print("\nHsap donor counts:", hs.obs["donor_id"].value_counts().to_dict())
print("\nHsap gene biotypes:", hs.var["feature_biotype"].value_counts().to_dict())
print("\nTotals: Dmel %d cells x %d genes; Hsap %d cells x %d genes" % (fl.n_obs, fl.n_vars, hs.n_obs, hs.n_vars))
EOF

python "${OUT_DIR}/07d_prepare_Dmel_Hsap.py"

echo; echo "=== 07d DONE ==="; date
