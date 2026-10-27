#!/usr/bin/env bash
# 05x_full_saturn_6sp_longest.sh   (2026-10-27)
# FULL SATURN training (pretraining + metric learning) on the 6-species panel, one run per protein embedding
# model, both built from ONE protein per gene (longest isoform; same protein for both models):
#     MODEL=ESM1b         03_embeddings_ESM1b_longest/{F}/{F}_gene_embeddings.pt            (05t, job 614089)
#     MODEL=ESMCpostnorm  03_embeddings_ESMC600M_postnorm_longest/{F}/{F}_gene_embeddings.pt  (05q --longest, 05w)
#     (F = Mnemi for Mlei)
# Everything else identical between the two runs and to the canonical full run (03_saturn_run.sh): K=3000
# macrogenes, 8000 HV genes, seed 42, SATURN default pretraining and metric-learning epochs, label column
# cell_type (Dmel/Hsap relabelled by 07g), CPU. Embeddings as delivered (no normalisation); truncation as
# embedded (ESM-1b 1,022 aa, ESMC 2,000 aa).
# Queue (10-27): gpu1 unavailable and gpu2 closed (running gpu_h100 jobs), so no host pin. Default = bigmem CPU,
# -W 168 h (bigmem allows 240 h; 6-sp pretraining alone took ~63 h on parallel nodes).
# GPU option: DEVICE=cuda with a GPU queue on the command line (needs a CUDA build of torch in saturn_protein), e.g.
#     MODEL=ESM1b DEVICE=cuda bsub -q gpu_arg -gpu "num=1" -W 48:00 ... < scripts/05x_full_saturn_6sp_longest.sh
# Both models must run on the same device type (CPU vs GPU results are not bit-identical).
# The pretrained model is saved (--pretrain_model_path), so later metric-learning variants (e.g. training on
# Leiden clusters instead of author labels) can reuse it instead of repeating ~2-3 days of pretraining.
# Fresh run directory: SATURN silently reuses an existing centroids_init.pkl, so the script refuses to start
# if one is already there.
# Submit (from /scratch/dark_genes/SATURN_Mnemi):
#     MODEL=ESM1b        bsub -J full6_ESM1b -o logs/full6_ESM1b.%J.out -e logs/full6_ESM1b.%J.err < scripts/05x_full_saturn_6sp_longest.sh
#     MODEL=ESMCpostnorm bsub -J full6_ESMCpn -o logs/full6_ESMCpn.%J.out -e logs/full6_ESMCpn.%J.err < scripts/05x_full_saturn_6sp_longest.sh
#BSUB -J full_6sp
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 168:00
#BSUB -R "rusage[mem=64000] span[hosts=1]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/full_6sp.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/full_6sp.%J.err

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
DEVICE="${DEVICE:-cpu}"
H="${BASE}/02_processed_h5ad"
case "${MODEL}" in
    ESM1b)        E="${BASE}/03_embeddings_ESM1b_longest" ;;
    ESMCpostnorm) E="${BASE}/03_embeddings_ESMC600M_postnorm_longest" ;;
    *)
        echo "Set MODEL=ESM1b or MODEL=ESMCpostnorm" >&2
        exit 1 ;;
esac

RUN_LABEL="20261027_Mlei_Cgig_Crob_Drer_Dmel_Hsap_K${K_MACROGENES}_${MODEL}_longest_full"
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
    echo "${H}/Mlei_saturn.h5ad,Mlei,${E}/Mnemi/Mnemi_gene_embeddings.pt"
    for SP in Cgig Crob Drer Dmel Hsap; do
        echo "${H}/${SP}_saturn.h5ad,${SP},${E}/${SP}/${SP}_gene_embeddings.pt"
    done
} > "${IN_DATA_CSV}"

echo "=== 05x FULL SATURN, 6 species, longest isoform, MODEL=${MODEL} ==="
date; hostname
echo "RUN_DIR=${RUN_DIR}"
echo "K=${K_MACROGENES}  HV=${HV_GENES}  seed=${MACROGENE_SEED}  epochs: SATURN defaults  device=${DEVICE}"
if [[ "${DEVICE}" == "cuda" ]]; then
    python -c "import torch, sys; ok = torch.cuda.is_available(); print('torch', torch.__version__, 'CUDA', torch.version.cuda, 'available', ok); sys.exit(0 if ok else 1)" \
        || { echo "DEVICE=cuda but torch sees no GPU" >&2; exit 1; }
fi
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
    --work_dir             "${RUN_DIR}" \
    --pretrain_model_path  "${RUN_DIR}/pretrain_model.model" \
    --device               "${DEVICE}" \
    --device_num           0

echo "${RUN_LABEL}" > "${BASE}/04_saturn_runs/LATEST_FULL_6sp_${MODEL}.txt"
ls -lh "${RUN_DIR}/saturn_results/" 2>/dev/null | head
echo "=== 05x done $(date) ==="
