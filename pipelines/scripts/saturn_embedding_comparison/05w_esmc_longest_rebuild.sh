#!/usr/bin/env bash
# 05w_esmc_longest_rebuild.sh   (2026-10-27)
# ESMC-600M post-norm gene vectors from ONE protein per gene (the protein ESM-1b's longest-isoform rebuild used),
# so the full SATURN runs compare ESM-1b and ESMC on the same isoform per gene.
#   1. 05q --longest  -> 03_embeddings_ESMC600M_postnorm_longest/ (6 species; Mlei under Mnemi/)
#                        + esmc_longest_choice_{SP}.tsv (how each gene's protein was chosen)
#   2. 05i            hit@1, ESM-1b longest vs ESMC longest (label "ESM1b" = 03_embeddings_ESM1b_longest here),
#                     6-sp ESM-1b run HV genes, vs 6-sp OrthoFinder OG. Sanity check: ESMC human pairs should
#                     improve like ESM-1b's did (ESMC averaged, 10-26: 0.276 raw cosine)
#   3. per-pair table
# Submit (from /scratch/dark_genes/SATURN_Mnemi):  bsub < scripts/05w_esmc_longest_rebuild.sh
#BSUB -J esmc_longest
#BSUB -q normal
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 06:00
#BSUB -R "rusage[mem=40000] span[hosts=1]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/esmc_longest.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/esmc_longest.%J.err

set -euo pipefail
SCRIPT_DIR="/scratch/dark_genes/SATURN_Mnemi/scripts"
source "${SCRIPT_DIR}/run_config.sh"
set +u
source ~/anaconda3/bin/activate
conda activate "${CONDA_ENV}"
set -u
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
cd "${BASE}"

RUNS="${BASE}/04_saturn_runs"
OF6="${BASE}/06_orthofinder_6sp"
SP6="Mlei,Cgig,Crob,Drer,Dmel,Hsap"
RUN_E="${RUNS}/20261019_Mlei_Cgig_Crob_Drer_Dmel_Hsap_K3000_ESM1b_macrogene"
OUTI="${RUNS}/20261027_longest_ESM1b_vs_ESMC_OG"
echo "=== 05w start $(date) on $(hostname) ==="

echo; echo "--- 1. 05q --longest: ESMC post-norm, one protein per gene ---"
python "${SCRIPT_DIR}/05q_convert_esmc_postnorm.py" --species all --longest --compare_esm1b

echo; echo "--- 2. 05i hit@1: ESM-1b longest vs ESMC longest ---"
python "${SCRIPT_DIR}/05i_embedding_standardization_check.py" --species "${SP6}" \
    --hv_run_a "${RUN_E}" --hv_run_b "${RUN_E}" \
    --esm1b_dir 03_embeddings_ESM1b_longest \
    --esmc_dir 03_embeddings_ESMC600M_postnorm_longest --esmc_label ESMC600M_postnorm_longest \
    --reference "${OF6}/ALL6_gene_orthofinder_OG.csv" \
    --out_dir "${OUTI}"

echo; echo "--- 3. per-pair table ---"
python "${SCRIPT_DIR}/05t_esm1b_longest_isoform.py" --report_pairs "${OUTI}/I2_standardization_by_species_pair.csv"

echo; echo "=== 05w done $(date) ==="
echo "Read: ${BASE}/03_embeddings_ESMC600M_postnorm_longest/esmc_longest_choice_*.tsv"
echo "      ${OUTI}/summary.txt"
