#!/usr/bin/env bash
# 05a_macrogene_only_ESM1b_4sp.sh
# Macrogene-only run: 4 species (Mlei, Cgig, Crob, Drer), ESM-1b embeddings.
# Same params as canonical 5-species run (K=3000, HV=8000, seed=42),
# --epochs 1 since macrogene assignment is fixed by pretraining alone
# (pretrain_epochs default 200 left untouched).
#BSUB -J mg_esm1b_4sp
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 24:00
#BSUB -R "rusage[mem=64000]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/mg_esm1b_4sp.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/mg_esm1b_4sp.%J.err

set -euo pipefail
SCRIPT_DIR="/scratch/dark_genes/SATURN_Mnemi/scripts"
BASE="/scratch/dark_genes/SATURN_Mnemi"
source "${SCRIPT_DIR}/run_config.sh"
# conda's activation scripts reference $PS1, which is unset in a
# non-interactive batch job -- relax nounset just for these two lines.
set +u
source ~/anaconda3/bin/activate
conda activate "${CONDA_ENV}"
set -u
cd "${SATURN_DIR}"

RUN_LABEL="20261002_Mlei_Cgig_Crob_Drer_K3000_ESM1b_macrogene"
RUN_DIR="${BASE}/04_saturn_runs/${RUN_LABEL}"
mkdir -p "${RUN_DIR}/logs"

IN_DATA_CSV="${RUN_DIR}/in_data.csv"
cat > "${IN_DATA_CSV}" <<EOF
path,species,embedding_path
${BASE}/02_processed_h5ad/Mlei_saturn.h5ad,Mlei,${BASE}/03_embeddings/Mlei/Mlei_gene_embeddings.pt
${BASE}/02_processed_h5ad/Cgig_saturn.h5ad,Cgig,${BASE}/03_embeddings/Cgig/Cgig_gene_embeddings.pt
${BASE}/02_processed_h5ad/Crob_saturn.h5ad,Crob,${BASE}/03_embeddings/Crob/Crob_gene_embeddings.pt
${BASE}/02_processed_h5ad/Drer_saturn.h5ad,Drer,${BASE}/03_embeddings/Drer/Drer_gene_embeddings.pt
EOF

# Sanity check: fail fast if any h5ad/embedding is missing, before launching
while IFS=, read -r H5AD SP EMB; do
    [[ "${H5AD}" == "path" ]] && continue
    if [[ ! -s "${H5AD}" ]]; then echo "MISSING h5ad for ${SP}: ${H5AD}" >&2; exit 1; fi
    if [[ ! -s "${EMB}" ]]; then echo "MISSING embedding for ${SP}: ${EMB}" >&2; exit 1; fi
done < "${IN_DATA_CSV}"

CENTROIDS_PKL="${RUN_DIR}/centroids_init.pkl"
python train-saturn.py \
    --in_data              "${IN_DATA_CSV}" \
    --in_label_col         "cell_type" \
    --ref_label_col        "cell_type" \
    --num_macrogenes       "${K_MACROGENES}" \
    --hv_genes              "${HV_GENES}" \
    --centroids_init_path  "${CENTROIDS_PKL}" \
    --embedding_model      "ESM1b" \
    --seed                 "${MACROGENE_SEED}" \
    --epochs               1 \
    --work_dir             "${RUN_DIR}" \
    --device               "cpu" \
    --device_num           0

echo "${RUN_LABEL}" > "${BASE}/04_saturn_runs/LATEST_MACROGENE_ESM1b.txt"
