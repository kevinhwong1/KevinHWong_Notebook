#!/usr/bin/env bash
# 07g_relabel_Dmel_Hsap.sh   (2026-10-19)
# Applies the label decisions agreed after the 07f QC (notebook 2026-10-15 human-and-fly entry, section 5):
#   Hsap: the 18 lymphocyte subtypes are collapsed to the 7 groups used for sampling in 07d (LYMPH_GROUPS):
#         CD4 T, CD8 T, gd T, MAIT, NK, B, Plasma. Myeloid, mast, progenitor and erythroid labels unchanged.
#         31 -> 20 labels.
#   Dmel: LM-1 (6 cells) + LM-2 (2 cells) -> LM. PL-ImpL2 (21 cells) KEPT (author plasmatocyte subtype;
#         QC doubts noted; do not interpret what it clusters with). 16 -> 15 labels.
#   Not changed: no cells removed, no depth downsampling (human depth: optional sensitivity run later).
# Only obs['cell_type'] changes (the column SATURN reads); the author label is kept in obs['cell_type_author'].
# Counts, genes and cells are untouched, so the embeddings and the gene lists are unaffected.
# Re-runnable: if cell_type_author already exists, cell_type is rebuilt from it.
# The first run keeps a backup of each input as <SP>_saturn.pre07g.h5ad.
# Submit:  bsub < /scratch/dark_genes/SATURN_Mnemi/scripts/07g_relabel_Dmel_Hsap.sh
#BSUB -J relabel_DmHs
#BSUB -q normal
#BSUB -P dark_genes
#BSUB -n 2
#BSUB -W 1:00
#BSUB -R "rusage[mem=16000] span[hosts=1]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/relabel_DmHs.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/relabel_DmHs.%J.err

set -euo pipefail
BASE="/scratch/dark_genes/SATURN_Mnemi"
H5AD_DIR="${BASE}/02_processed_h5ad"
PY_DIR="${BASE}/scripts/py"
mkdir -p "${PY_DIR}"

set +u
source ~/anaconda3/bin/activate
conda activate saturn_protein
set -u

cat > "${PY_DIR}/07g_relabel_Dmel_Hsap.py" << 'PYEOF'
import shutil
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

H5AD_DIR = Path(sys.argv[1])

# same groups as LYMPH_GROUPS in 07d_prepare_Dmel_Hsap.sh
LYMPH_GROUPS = {
    "CD4 T":  ["Tnaive/CM_CD4", "Teffector/EM_CD4", "Trm_Th1/Th17", "Tfh", "Tregs"],
    "CD8 T":  ["Tnaive/CM_CD8", "Tem/emra_CD8", "Trm_gut_CD8", "Trm/em_CD8"],
    "gd T":   ["Trm_Tgd", "Tgd_CRTAM+"],
    "MAIT":   ["MAIT"],
    "NK":     ["NK_CD16+", "NK_CD56bright_CD16-"],
    "B":      ["Naive B cells", "Memory B cells"],
    "Plasma": ["Plasma cells", "Plasmablasts"],
}
MAPS = {
    "Hsap": {sub: grp for grp, subs in LYMPH_GROUPS.items() for sub in subs},
    "Dmel": {"LM-1": "LM", "LM-2": "LM"},
}
EXPECTED_N_LABELS = {"Hsap": 20, "Dmel": 15}

rows = []
for sp, mapping in MAPS.items():
    path = H5AD_DIR / f"{sp}_saturn.h5ad"
    backup = H5AD_DIR / f"{sp}_saturn.pre07g.h5ad"
    if not backup.exists():
        shutil.copy2(path, backup)
        print(f"{sp}: backup -> {backup.name}")
    a = ad.read_h5ad(path)
    if "cell_type_author" not in a.obs.columns:
        a.obs["cell_type_author"] = a.obs["cell_type"].astype(str)
    author = a.obs["cell_type_author"].astype(str)
    unused = sorted(set(mapping) - set(author.unique()))
    if unused:
        raise SystemExit(f"{sp}: mapping labels not found in the data: {unused}")
    new = author.map(lambda x: mapping.get(x, x))
    a.obs["cell_type"] = pd.Categorical(new)
    n_before, n_after = author.nunique(), a.obs["cell_type"].nunique()
    print(f"{sp}: {a.n_obs} cells x {a.n_vars} genes; labels {n_before} -> {n_after}")
    if n_after != EXPECTED_N_LABELS[sp]:
        raise SystemExit(f"{sp}: expected {EXPECTED_N_LABELS[sp]} labels, got {n_after}; nothing written")
    tab = (a.obs.groupby(["cell_type", "cell_type_author"], observed=True).size()
           .rename("n_cells").reset_index())
    tab.insert(0, "species", sp)
    rows.append(tab)
    print(tab.to_string(index=False))
    a.write_h5ad(path)
    print(f"{sp}: wrote {path.name}\n")

out = H5AD_DIR / "Dmel_Hsap_relabel_07g.csv"
pd.concat(rows).to_csv(out, index=False)
print(f"label table -> {out}")
PYEOF

echo "=== 07g start $(date) on $(hostname) ==="
python "${PY_DIR}/07g_relabel_Dmel_Hsap.py" "${H5AD_DIR}"
echo "=== 07g done $(date) ==="
