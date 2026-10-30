#!/usr/bin/env bash
# 05y0_prep_avg_matched.sh   (2026-10-29)
# Step 0 of the 2x2 full-SATURN comparison (ESM-1b / ESMC x averaged / longest isoform).
#   1. 05q --match_esm1b : averaged ESMC post-norm gene vectors over exactly the isoforms ESM-1b averaged, ESM-1b
#                          gene set only, all 6 species -> 03_embeddings_ESMC600M_postnorm_avgmatched/
#                          (+ esmc_avgmatch_{SP}.tsv). The old 03_embeddings_ESMC600M_postnorm/ is left untouched
#                          (10-19 runs stay reproducible).
#   2. 05y_check_gene_sets.py : all 4 embedding sets give identical usable gene sets per species (-> identical HV
#                          genes); finite values; md5s. Exits 1 on FAIL, so a dependent job (-w "done(...)") won't start.
#   3. 05i hit@1, ESM-1b averaged vs ESMC averaged-matched (sanity: ESMC averaged was 0.276 raw cosine on 10-26;
#                          matching should change it only slightly, mostly Drer pairs)
# Submit (from /scratch/dark_genes/SATURN_Mnemi):  bsub < scripts/05y0_prep_avg_matched.sh
#BSUB -J y0_prep
#BSUB -q normal
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 06:00
#BSUB -R "rusage[mem=40000] span[hosts=1]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/y0_prep.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/y0_prep.%J.err

set -euo pipefail
SCRIPT_DIR="/scratch/dark_genes/SATURN_Mnemi/scripts"
source "${SCRIPT_DIR}/run_config.sh"
set +u
source ~/anaconda3/bin/activate
conda activate "${CONDA_ENV}"
set -u
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
cd "${BASE}"

OUTD="${BASE}/03_embeddings_ESMC600M_postnorm_avgmatched"
if [[ -e "${OUTD}" ]]; then
    echo "${OUTD} exists -- move it aside first (this script does not overwrite)." >&2
    exit 1
fi
RUNS="${BASE}/04_saturn_runs"
OF6="${BASE}/06_orthofinder_6sp"
SP6="Mlei,Cgig,Crob,Drer,Dmel,Hsap"
RUN_E="${RUNS}/20261019_Mlei_Cgig_Crob_Drer_Dmel_Hsap_K3000_ESM1b_macrogene"
OUTI="${RUNS}/20261029_avg_ESM1b_vs_ESMCmatched_OG"
echo "=== 05y0 start $(date) on $(hostname) ==="

echo; echo "--- 1. 05q --match_esm1b: ESMC post-norm, averaged over ESM-1b's isoforms ---"
python "${SCRIPT_DIR}/05q_convert_esmc_postnorm.py" --species all --match_esm1b --compare_esm1b

echo; echo "--- 2. gene-set check across the 4 embedding sets ---"
python "${SCRIPT_DIR}/05y_check_gene_sets.py"

echo; echo "--- 3. 05i hit@1: ESM-1b averaged vs ESMC averaged-matched ---"
python "${SCRIPT_DIR}/05i_embedding_standardization_check.py" --species "${SP6}" \
    --hv_run_a "${RUN_E}" --hv_run_b "${RUN_E}" \
    --esm1b_dir 03_embeddings \
    --esmc_dir 03_embeddings_ESMC600M_postnorm_avgmatched --esmc_label ESMC600M_postnorm_avgmatched \
    --reference "${OF6}/ALL6_gene_orthofinder_OG.csv" \
    --out_dir "${OUTI}"
python "${SCRIPT_DIR}/05t_esm1b_longest_isoform.py" --report_pairs "${OUTI}/I2_standardization_by_species_pair.csv" || true

echo; echo "=== 05y0 done $(date) ==="
echo "Read: ${OUTD}/05y_gene_set_check.tsv, ${OUTD}/esmc_avgmatch_*.tsv, ${OUTI}/summary.txt"
