#!/usr/bin/env bash
# 05u_isoform_longest_test.sh   (2026-10-27)
# Isoform test: ESM-1b gene vectors from the longest isoform vs all isoforms averaged (current), 6 species.
#   1. 05t  rebuild 03_embeddings_ESM1b_longest/ from the existing per-protein .pt files (no new embedding)
#   2. 05i  nearest-neighbour ortholog hit@1, ESM1b (averaged) vs ESM1b_longest, on the 6-sp ESM-1b run's
#           HV genes, vs the 6-sp OrthoFinder OG table (same settings as 06e step 3)
#   3. 05t --report_pairs  per-pair table; human pairs side by side
# Submit (from /scratch/dark_genes/SATURN_Mnemi):  bsub < scripts/05u_isoform_longest_test.sh
#BSUB -J iso_longest
#BSUB -q normal
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 06:00
#BSUB -R "rusage[mem=40000] span[hosts=1]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/iso_longest.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/iso_longest.%J.err

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
OUTI="${RUNS}/20261027_isoform_longest_OG"
echo "=== 05u start $(date) on $(hostname) ==="

echo; echo "--- 1. 05t rebuild longest-isoform gene vectors ---"
python "${SCRIPT_DIR}/05t_esm1b_longest_isoform.py" --species "${SP6}"

echo; echo "--- 2. 05i hit@1: ESM1b averaged vs ESM1b_longest ---"
python "${SCRIPT_DIR}/05i_embedding_standardization_check.py" --species "${SP6}" \
    --hv_run_a "${RUN_E}" --hv_run_b "${RUN_E}" \
    --esmc_dir 03_embeddings_ESM1b_longest --esmc_label ESM1b_longest \
    --reference "${OF6}/ALL6_gene_orthofinder_OG.csv" \
    --out_dir "${OUTI}"

echo; echo "--- 3. per-pair table ---"
python "${SCRIPT_DIR}/05t_esm1b_longest_isoform.py" --report_pairs "${OUTI}/I2_standardization_by_species_pair.csv"

echo; echo "=== 05u done $(date) ==="
echo "Read: ${BASE}/03_embeddings_ESM1b_longest/longest_isoform_summary.tsv"
echo "      ${OUTI}/summary.txt"
