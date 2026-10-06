#!/usr/bin/env bash
# 05b2_macrogene_only_ESMC600M_rescaled_4sp.sh
# Macrogene-only run: 4 species (Mlei, Cgig, Crob, Drer), ESMC-600M embeddings, LENGTH-NORMALIZED (05e_rescale_esmc_embeddings.py)
# (converted from collaborator's FANTASIA v4 HDF5 output).
# Same params as canonical 5-species run (K=3000, HV=8000, seed=42),
# --epochs 1 since macrogene assignment is fixed by pretraining alone.
# NOTE: --embedding_model is decorative here -- SATURN's loader fully
# overrides it with the per-species embedding_path column below (no
# "ESMC" choice exists in the CLI's fixed --embedding_model list).
# NOTE: ESMC output folder for Mnemiopsis leidyi is named "Mnemi" (FANTASIA's
# label) not "Mlei" -- confirmed same species as the Mlei row in the h5ad.
#BSUB -J mg_esmcR_4sp
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 72:00
#BSUB -R "rusage[mem=64000]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/mg_esmcR_4sp.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/mg_esmcR_4sp.%J.err

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

RUN_LABEL="20261005_Mlei_Cgig_Crob_Drer_K3000_ESMC600M_rescaled_macrogene"
RUN_DIR="${BASE}/04_saturn_runs/${RUN_LABEL}"
mkdir -p "${RUN_DIR}/logs"

IN_DATA_CSV="${RUN_DIR}/in_data.csv"
cat > "${IN_DATA_CSV}" <<EOF
path,species,embedding_path
${BASE}/02_processed_h5ad/Mlei_saturn.h5ad,Mlei,${BASE}/03_embeddings_ESMC600M_rescaled/Mnemi/Mnemi_gene_embeddings.pt
${BASE}/02_processed_h5ad/Cgig_saturn.h5ad,Cgig,${BASE}/03_embeddings_ESMC600M_rescaled/Cgig/Cgig_gene_embeddings.pt
${BASE}/02_processed_h5ad/Crob_saturn.h5ad,Crob,${BASE}/03_embeddings_ESMC600M_rescaled/Crob/Crob_gene_embeddings.pt
${BASE}/02_processed_h5ad/Drer_saturn.h5ad,Drer,${BASE}/03_embeddings_ESMC600M_rescaled/Drer/Drer_gene_embeddings.pt
EOF

while IFS=, read -r H5AD SP EMB; do
    [[ "${H5AD}" == "path" ]] && continue
    if [[ ! -s "${H5AD}" ]]; then echo "MISSING h5ad for ${SP}: ${H5AD}" >&2; exit 1; fi
    if [[ ! -s "${EMB}" ]]; then echo "MISSING embedding for ${SP}: ${EMB}" >&2; exit 1; fi
done < "${IN_DATA_CSV}"

CENTROIDS_PKL="${RUN_DIR}/centroids_init.pkl"
# SATURN silently LOADS an existing centroids_init.pkl and skips KMeans --
# refuse to run on a stale one so the rescaled embeddings actually take effect.
if [[ -e "${CENTROIDS_PKL}" ]]; then
    echo "ERROR: ${CENTROIDS_PKL} already exists; delete it (or the run dir) first." >&2
    exit 1
fi
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

echo "${RUN_LABEL}" > "${BASE}/04_saturn_runs/LATEST_MACROGENE_ESMC600M_rescaled.txt"
