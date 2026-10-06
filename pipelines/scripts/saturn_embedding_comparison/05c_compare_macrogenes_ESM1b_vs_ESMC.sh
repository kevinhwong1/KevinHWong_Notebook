#!/usr/bin/env bash
# 05c_compare_macrogenes_ESM1b_vs_ESMC.sh
# Label-invariant comparison of macrogene assignments between the 4-species
# ESM-1b and ESMC-600M macrogene-only runs (2026-10-02), anchored against the
# canonical 5-species ESM-1b run (20260617). Reads genes_to_macrogenes.pkl only.
# Python lives in its own file (05c_...py) to avoid heredoc pipe truncation.
#BSUB -J mg_compare
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 04:00
#BSUB -R "rusage[mem=32000]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/mg_compare.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/mg_compare.%J.err

set -euo pipefail
SCRIPT_DIR="/scratch/dark_genes/SATURN_Mnemi/scripts"
source "${SCRIPT_DIR}/run_config.sh"
set +u
source ~/anaconda3/bin/activate
conda activate "${CONDA_ENV}"
set -u

export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8

python "${SCRIPT_DIR}/05c_compare_macrogenes_ESM1b_vs_ESMC.py" \
    --species "Mlei,Cgig,Crob,Drer" \
    --focus_ref_mgs "1089,2003" \
    --knn_k 10
