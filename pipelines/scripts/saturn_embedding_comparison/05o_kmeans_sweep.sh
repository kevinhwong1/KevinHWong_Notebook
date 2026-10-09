#!/usr/bin/env bash
# 05o_kmeans_sweep.sh
# K sweep of SATURN's macrogene-seeding KMeans (no SATURN run), ESM-1b vs ESMC-600M (rescaled, as in the 20261005 run).
# Array: 1-12 = 2 models x 6 K values; 13 = validation (init vs final macrogenes; reproduce centroids_init.pkl).
# Submit:  bsub < /scratch/dark_genes/SATURN_Mnemi/scripts/05o_kmeans_sweep.sh
# Re-run single elements:  bsub -J "ksweep[13]" < 05o_kmeans_sweep.sh   (or e.g. "ksweep[6,12]")
#BSUB -J "ksweep[1-13]"
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 24:00
#BSUB -R "rusage[mem=32000]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/ksweep.%J.%I.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/ksweep.%J.%I.err

set -euo pipefail
SCRIPT_DIR="/scratch/dark_genes/SATURN_Mnemi/scripts"
source "${SCRIPT_DIR}/run_config.sh"
set +u
source ~/anaconda3/bin/activate
conda activate "${CONDA_ENV}"
set -u
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8

MODELS=(ESM1b ESMC600M_rescaled)
KS=(500 1000 2000 3000 5000 10000)
SEED=42
I=${LSB_JOBINDEX}

if [[ "${I}" -eq 13 ]]; then
    python -u "${SCRIPT_DIR}/05o_kmeans_sweep.py" --mode validate --seed "${SEED}"
else
    M=${MODELS[$(( (I - 1) / ${#KS[@]} ))]}
    K=${KS[$(( (I - 1) % ${#KS[@]} ))]}
    echo "Array element ${I}: model=${M} K=${K} seed=${SEED}"
    python -u "${SCRIPT_DIR}/05o_kmeans_sweep.py" --mode sweep --model "${M}" --k "${K}" --seed "${SEED}"
fi
