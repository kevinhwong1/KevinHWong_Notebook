#!/usr/bin/env bash
# 05y_full_saturn_6sp.sh   (2026-10-29)
# FULL SATURN (pretraining + metric learning), 6-species panel, for the 2x2 comparison
#   model (ESM-1b / ESMC-600M post-norm)  x  protein per gene (averaged isoforms / longest isoform)
# Generalises 05x (left untouched: jobs 614411/614412 = ESM1b/ESMC longest, seed 42, are runs 3 and 4).
#
#   MODEL=ESM1b        ISOFORM=avg      03_embeddings/{SP}/                          (Mlei: 03_embeddings/Mlei/)
#   MODEL=ESMCpostnorm ISOFORM=avg      03_embeddings_ESMC600M_postnorm_avgmatched/{F}/   (05q --match_esm1b, 05y0)
#   MODEL=ESM1b        ISOFORM=longest  03_embeddings_ESM1b_longest/{F}/                  (as 05x)
#   MODEL=ESMCpostnorm ISOFORM=longest  03_embeddings_ESMC600M_postnorm_longest/{F}/      (as 05x)
#   SEED   default 42 (seed replicate: SEED=43; changes k-means macrogene init and training)
#   SMOKE=1  pipeline test: --pretrain_epochs 1 --epochs 1, run dir under 04_saturn_runs/_smoke/ (~1-2 h).
#            Checks that the pretrained model is saved and metric learning + output writing work end to end.
#
# Identical to 05x / the canonical full run: K=3000, HV=8000, label column cell_type (Dmel/Hsap relabelled by 07g),
# SATURN default pretraining (200) and metric-learning epochs, embeddings as delivered (no normalisation),
# each model's own truncation (ESM-1b 1,022 aa; ESMC 2,000 aa), CPU, bigmem, -n 8.
# Difference from 05x: thread count pinned to the 8 requested slots (OMP/MKL/OPENBLAS), so up to 4-5 jobs can
# share a bigmem node without oversubscribing. Changes only floating-point summation order, not the model.
# Writes provenance.txt (md5 of train-saturn.py, h5ads, embeddings; torch threads; conda list) to the run dir.
#
# Submit (from /scratch/dark_genes/SATURN_Mnemi), e.g.
#   SMOKE=1 MODEL=ESMCpostnorm ISOFORM=avg bsub -J y_smoke -W 06:00 -o logs/y_smoke.%J.out -e logs/y_smoke.%J.err < scripts/05y_full_saturn_6sp.sh
#   MODEL=ESM1b        ISOFORM=avg              bsub -J full6_ESM1b_avg   -o logs/full6_ESM1b_avg.%J.out   -e logs/full6_ESM1b_avg.%J.err   < scripts/05y_full_saturn_6sp.sh
#   MODEL=ESMCpostnorm ISOFORM=avg              bsub -J full6_ESMCpn_avg  -o logs/full6_ESMCpn_avg.%J.out  -e logs/full6_ESMCpn_avg.%J.err  < scripts/05y_full_saturn_6sp.sh
#   MODEL=ESM1b        ISOFORM=longest SEED=43  bsub -J full6_ESM1b_s43   -o logs/full6_ESM1b_s43.%J.out   -e logs/full6_ESM1b_s43.%J.err   < scripts/05y_full_saturn_6sp.sh
#BSUB -J full_6sp_y
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 168:00
#BSUB -R "rusage[mem=64000] span[hosts=1]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/full_6sp_y.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/full_6sp_y.%J.err

set -euo pipefail
SCRIPT_DIR="/scratch/dark_genes/SATURN_Mnemi/scripts"
BASE="/scratch/dark_genes/SATURN_Mnemi"
source "${SCRIPT_DIR}/run_config.sh"
set +u
source ~/anaconda3/bin/activate
conda activate "${CONDA_ENV}"
set -u
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
cd "${SATURN_DIR}"

MODEL="${MODEL:-}"
ISOFORM="${ISOFORM:-}"
SEED="${SEED:-${MACROGENE_SEED}}"
SMOKE="${SMOKE:-0}"
DEVICE="cpu"
H="${BASE}/02_processed_h5ad"

case "${MODEL}:${ISOFORM}" in
    ESM1b:avg)            E="${BASE}/03_embeddings" ;;
    ESMCpostnorm:avg)     E="${BASE}/03_embeddings_ESMC600M_postnorm_avgmatched" ;;
    ESM1b:longest)        E="${BASE}/03_embeddings_ESM1b_longest" ;;
    ESMCpostnorm:longest) E="${BASE}/03_embeddings_ESMC600M_postnorm_longest" ;;
    *)
        echo "Set MODEL=ESM1b|ESMCpostnorm and ISOFORM=avg|longest (got '${MODEL}:${ISOFORM}')" >&2
        exit 1 ;;
esac
# Mlei folder: the ESM-1b averaged set uses Mlei/Mlei_*, the three others Mnemi/Mnemi_*
if [[ "${MODEL}:${ISOFORM}" == "ESM1b:avg" ]]; then
    EMB_Mlei="${E}/Mlei/Mlei_gene_embeddings.pt"
else
    EMB_Mlei="${E}/Mnemi/Mnemi_gene_embeddings.pt"
fi

RUN_LABEL="$(date +%Y%m%d)_Mlei_Cgig_Crob_Drer_Dmel_Hsap_K${K_MACROGENES}_${MODEL}_${ISOFORM}_seed${SEED}_full"
EPOCH_ARGS=()
EPOCH_DESC="SATURN defaults"
if [[ "${SMOKE}" == "1" ]]; then
    RUN_LABEL="${RUN_LABEL}_SMOKE"
    RUN_DIR="${BASE}/04_saturn_runs/_smoke/${RUN_LABEL}"
    EPOCH_ARGS=(--pretrain_epochs 1 --epochs 1)
    EPOCH_DESC="--pretrain_epochs 1 --epochs 1 (SMOKE)"
    python train-saturn.py --help 2>/dev/null | grep -q -- "--pretrain_epochs" \
        || { echo "train-saturn.py has no --pretrain_epochs flag; smoke test cannot shorten pretraining" >&2; exit 1; }
else
    RUN_DIR="${BASE}/04_saturn_runs/${RUN_LABEL}"
fi
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

echo "=== 05y FULL SATURN, 6 species, MODEL=${MODEL} ISOFORM=${ISOFORM} SEED=${SEED} SMOKE=${SMOKE} ==="
date; hostname
echo "RUN_DIR=${RUN_DIR}"
echo "K=${K_MACROGENES}  HV=${HV_GENES}  seed=${SEED}  epochs: ${EPOCH_DESC}  device=${DEVICE}"
cat "${IN_DATA_CSV}"

# fail fast if any input is missing
while IFS=, read -r H5AD SP EMB; do
    [[ "${H5AD}" == "path" ]] && continue
    [[ -s "${H5AD}" ]] || { echo "MISSING h5ad for ${SP}: ${H5AD}" >&2; exit 1; }
    [[ -s "${EMB}" ]]  || { echo "MISSING embedding for ${SP}: ${EMB}" >&2; exit 1; }
done < "${IN_DATA_CSV}"

# provenance
PROV="${RUN_DIR}/provenance.txt"
{
    echo "run ${RUN_LABEL}   host $(hostname)   start $(date)"
    echo "MODEL=${MODEL} ISOFORM=${ISOFORM} SEED=${SEED} SMOKE=${SMOKE} K=${K_MACROGENES} HV=${HV_GENES}"
    echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}  torch threads: $(python -c 'import torch; print(torch.get_num_threads(), torch.__version__)')"
    echo "--- md5"
    md5sum train-saturn.py
    tail -n +2 "${IN_DATA_CSV}" | while IFS=, read -r H5AD SP EMB; do md5sum "${H5AD}" "${EMB}"; done
    echo "--- conda list"
    conda list
} > "${PROV}" 2>&1
head -12 "${PROV}"

# --embedding_model is decorative when embedding_path is given per species (see 05b notes)
python train-saturn.py \
    --in_data              "${IN_DATA_CSV}" \
    --in_label_col         "cell_type" \
    --ref_label_col        "cell_type" \
    --num_macrogenes       "${K_MACROGENES}" \
    --hv_genes             "${HV_GENES}" \
    --centroids_init_path  "${CENTROIDS_PKL}" \
    --embedding_model      "ESM1b" \
    --seed                 "${SEED}" \
    --work_dir             "${RUN_DIR}" \
    --pretrain_model_path  "${RUN_DIR}/pretrain_model.model" \
    ${EPOCH_ARGS[@]+"${EPOCH_ARGS[@]}"} \
    --device               "${DEVICE}" \
    --device_num           0

if [[ "${SMOKE}" == "1" ]]; then
    echo "--- smoke-test checks"
    ok=1
    [[ -s "${RUN_DIR}/pretrain_model.model" ]] && echo "pretrain_model.model saved: OK" || { echo "pretrain_model.model MISSING"; ok=0; }
    n_h5ad=$(find "${RUN_DIR}/saturn_results" -name "*.h5ad" 2>/dev/null | wc -l)
    [[ "${n_h5ad}" -ge 1 ]] && echo "saturn_results h5ad written (${n_h5ad}): OK" || { echo "no h5ad in saturn_results"; ok=0; }
    ls -lh "${RUN_DIR}/saturn_results/" 2>/dev/null | head
    [[ "${ok}" == "1" ]] || { echo "SMOKE TEST FAILED" >&2; exit 1; }
    echo "SMOKE TEST PASSED"
else
    echo "${RUN_LABEL}" > "${BASE}/04_saturn_runs/LATEST_FULL_6sp_${MODEL}_${ISOFORM}_seed${SEED}.txt"
    ls -lh "${RUN_DIR}/saturn_results/" 2>/dev/null | head
fi
echo "=== 05y done $(date) ==="
