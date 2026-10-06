#!/usr/bin/env bash
# 05h_norm_vs_macrogene.sh
# ESMC vector length vs macrogene assignment across ESM-1b / ESMC raw / ESMC rescaled runs (run after 05b2 finishes).
#BSUB -J norm_mg
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 03:00
#BSUB -R "rusage[mem=32000]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/norm_mg.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/norm_mg.%J.err

set -euo pipefail
SCRIPT_DIR="/scratch/dark_genes/SATURN_Mnemi/scripts"
source "${SCRIPT_DIR}/run_config.sh"
set +u
source ~/anaconda3/bin/activate
conda activate "${CONDA_ENV}"
set -u
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8

python "${SCRIPT_DIR}/05h_norm_vs_macrogene.py"
