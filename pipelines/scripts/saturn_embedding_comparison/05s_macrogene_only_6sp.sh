#!/usr/bin/env bash
# 05s_macrogene_only_6sp.sh   (2026-10-19)
# Phase 3 of PLAN_6species_OrthoFinder_ESM1b_vs_ESMCpostnorm_20261019: macrogene-only SATURN runs on the
# 6-species panel (Mlei, Cgig, Crob, Drer, Dmel, Hsap), one per protein embedding model:
#     MODEL=ESM1b         03_embeddings/{SP}/{SP}_gene_embeddings.pt
#     MODEL=ESMCpostnorm  03_embeddings_ESMC600M_postnorm/{F}/{F}_gene_embeddings.pt   (05q; F = Mnemi for Mlei)
# One script for both so the only difference between the runs is the embedding file.
# Settings identical to the 4-species runs 05a/05b (2026-10-02): K=3000 macrogenes, 8000 HV genes, seed 42,
# pretrain_epochs left at SATURN's default (200), --epochs 1 (macrogene membership is fixed at pretraining,
# which uses expression only; labels enter only in metric learning), CPU.
# Inputs: Dmel/Hsap h5ads relabelled by 07g (2026-10-19); cells, counts and genes unchanged from 07d.
# Dmel/Hsap ESMC gene sets are identical to ESM-1b (05q restricted to the ESM-1b protein set).
# Speed: the 4-sp runs took ~102 s/epoch on gpu1 (~6 h of pretraining) vs ~630 s/epoch on n051-type nodes,
# so this is pinned to gpu1/gpu2 (bigmem), as in the October runs. 6 species = ~89k cells (4 sp: ~70k).
# A fresh run directory is required: SATURN silently reuses an existing centroids_init.pkl, so the script
# refuses to start if one is already there.
# Submit (from /scratch/dark_genes/SATURN_Mnemi):
#     MODEL=ESM1b        bsub < scripts/05s_macrogene_only_6sp.sh
#     MODEL=ESMCpostnorm bsub < scripts/05s_macrogene_only_6sp.sh
#BSUB -J mg_6sp
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 72:00
#BSUB -R "rusage[mem=64000] select[hname==gpu1 || hname==gpu2]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/mg_6sp.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/mg_6sp.%J.err

set -euo pipefail
SCRIPT_DIR="/scratch/dark_genes/SATURN_Mnemi/scripts"
BASE="/scratch/dark_genes/SATURN_Mnemi"
source "${SCRIPT_DIR}/run_config.sh"
set +u
source ~/anaconda3/bin/activate
conda activate "${CONDA_ENV}"
set -u
cd "${SATURN_DIR}"

MODEL="${MODEL:-}"
H="${BASE}/02_processed_h5ad"
case "${MODEL}" in
    ESM1b)
        E="${BASE}/03_embeddings"
        EMB_Mlei="${E}/Mlei/Mlei_gene_embeddings.pt" ;;
    ESMCpostnorm)
        E="${BASE}/03_embeddings_ESMC600M_postnorm"
        EMB_Mlei="${E}/Mnemi/Mnemi_gene_embeddings.pt" ;;
    *)
        echo "Set MODEL=ESM1b or MODEL=ESMCpostnorm (e.g. MODEL=ESM1b bsub < scripts/05s_macrogene_only_6sp.sh)" >&2
        exit 1 ;;
esac

RUN_LABEL="20261019_Mlei_Cgig_Crob_Drer_Dmel_Hsap_K${K_MACROGENES}_${MODEL}_macrogene"
RUN_DIR="${BASE}/04_saturn_runs/${RUN_LABEL}"
CENTROIDS_PKL="${RUN_DIR}/centroids_init.pkl"
if [[ -e "${CENTROIDS_PKL}" ]]; then
    echo "${CENTROIDS_PKL} exists -- SATURN would reuse it. Move/delete ${RUN_DIR} first." >&2
    exit 1
fi
mkdir -p "${RUN_DIR}/logs"

IN_DATA_CSV="${RUN_DIR}/in_data.csv"
{
    echo "path,species,embedding_path"
    echo "${H}/Mlei_saturn.h5ad,Mlei,${EMB_Mlei}"
    for SP in Cgig Crob Drer Dmel Hsap; do
        echo "${H}/${SP}_saturn.h5ad,${SP},${E}/${SP}/${SP}_gene_embeddings.pt"
    done
} > "${IN_DATA_CSV}"

echo "=== 05s macrogene-only, 6 species, MODEL=${MODEL} ==="
date; hostname
echo "RUN_DIR=${RUN_DIR}"
echo "K=${K_MACROGENES}  HV=${HV_GENES}  seed=${MACROGENE_SEED}"
cat "${IN_DATA_CSV}"

# fail fast if any input is missing
while IFS=, read -r H5AD SP EMB; do
    [[ "${H5AD}" == "path" ]] && continue
    [[ -s "${H5AD}" ]] || { echo "MISSING h5ad for ${SP}: ${H5AD}" >&2; exit 1; }
    [[ -s "${EMB}" ]]  || { echo "MISSING embedding for ${SP}: ${EMB}" >&2; exit 1; }
done < "${IN_DATA_CSV}"

# --embedding_model is decorative when embedding_path is given per species (see 05b notes)
python train-saturn.py \
    --in_data              "${IN_DATA_CSV}" \
    --in_label_col         "cell_type" \
    --ref_label_col        "cell_type" \
    --num_macrogenes       "${K_MACROGENES}" \
    --hv_genes             "${HV_GENES}" \
    --centroids_init_path  "${CENTROIDS_PKL}" \
    --embedding_model      "ESM1b" \
    --seed                 "${MACROGENE_SEED}" \
    --epochs               1 \
    --work_dir             "${RUN_DIR}" \
    --device               "cpu" \
    --device_num           0

echo "${RUN_LABEL}" > "${BASE}/04_saturn_runs/LATEST_MACROGENE_6sp_${MODEL}.txt"
ls -lh "${RUN_DIR}/saturn_results/" 2>/dev/null | head
echo "=== 05s done $(date) ==="
