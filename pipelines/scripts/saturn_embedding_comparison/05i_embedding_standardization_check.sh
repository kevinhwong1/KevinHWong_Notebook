#!/usr/bin/env bash
# 05i_embedding_standardization_check.sh
# Can post-processing (unit length, z-score per dimension, removing top PCs) rescue ESMC? Nearest-neighbour ortholog test, no SATURN run.
#BSUB -J emb_std
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 03:00
#BSUB -R "rusage[mem=48000]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/emb_std.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/emb_std.%J.err

set -euo pipefail
SCRIPT_DIR="/scratch/dark_genes/SATURN_Mnemi/scripts"
source "${SCRIPT_DIR}/run_config.sh"
set +u
source ~/anaconda3/bin/activate
conda activate "${CONDA_ENV}"
set -u
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8

python "${SCRIPT_DIR}/05i_embedding_standardization_check.py"
