#!/usr/bin/env bash
# 05f_embedding_descriptives.sh
# Descriptive stats + figures for ESM-1b / ESMC raw / ESMC rescaled gene embeddings (run after 05e).
#BSUB -J emb_desc
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 03:00
#BSUB -R "rusage[mem=32000]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/emb_desc.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/emb_desc.%J.err

set -euo pipefail
SCRIPT_DIR="/scratch/dark_genes/SATURN_Mnemi/scripts"
source "${SCRIPT_DIR}/run_config.sh"
set +u
source ~/anaconda3/bin/activate
conda activate "${CONDA_ENV}"
set -u
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8

python "${SCRIPT_DIR}/05f_embedding_descriptives.py"
