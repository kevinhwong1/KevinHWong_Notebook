#!/usr/bin/env bash
# 05d_embedding_species_structure.sh
# Pre-SATURN check of raw ESM-1b vs ESMC-600M gene embeddings for the shared
# 4-species HV genes: integrity, species separation, cross-species kNN mixing,
# and within-species cross-model agreement (conversion-artifact check).
#BSUB -J emb_species
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 02:00
#BSUB -R "rusage[mem=32000]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/emb_species.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/emb_species.%J.err

set -euo pipefail
SCRIPT_DIR="/scratch/dark_genes/SATURN_Mnemi/scripts"
source "${SCRIPT_DIR}/run_config.sh"
set +u
source ~/anaconda3/bin/activate
conda activate "${CONDA_ENV}"
set -u

export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8

python "${SCRIPT_DIR}/05d_embedding_species_structure.py" \
    --species "Mlei,Cgig,Crob,Drer" \
    --knn_k 10
