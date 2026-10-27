#!/usr/bin/env bash
# 06e_benchmark_6sp.sh   (2026-10-19)
# Phase 4 of PLAN_6species_OrthoFinder_ESM1b_vs_ESMCpostnorm_20261019: score the 6-species macrogene-only runs
# (05s) and embeddings against the 6-species OrthoFinder reference (06b6). Run after both have finished.
# Steps
#   1. 06c  OrthoFinder -> 06_orthofinder_6sp/ALL6_gene_orthofinder_OG.csv (+ ALL6_..._<Metazoa node>HOG.csv);
#           outgroup Sros dropped; Metazoa node found from the species tree (expected N1)
#   2. 05g  run descriptives + ortholog benchmark (recall / precision / OG recovery per species pair; pooled,
#           Mlei pairs, Drer-Hsap positive control, invertebrate pairs) vs the OG table, and vs the HOG table
#   3. 05i  nearest-neighbour ortholog hit@1 (raw / unit / zscore / abtt1-3) for ESM-1b vs ESMC postnorm on the
#           6-species HV genes, vs the OG table
# The ESMC postnorm SATURN run is optional: if its run folder has no genes_to_macrogenes pkl, step 2 scores the
# ESM-1b run alone and step 3 uses the ESM-1b run's HV genes for both models.
# Submit (from /scratch/dark_genes/SATURN_Mnemi):  bsub < scripts/06e_benchmark_6sp.sh
#BSUB -J bench_6sp
#BSUB -q normal
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 08:00
#BSUB -R "rusage[mem=32000] span[hosts=1]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/bench_6sp.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/bench_6sp.%J.err

set -euo pipefail
SCRIPT_DIR="/scratch/dark_genes/SATURN_Mnemi/scripts"
source "${SCRIPT_DIR}/run_config.sh"
set +u
source ~/anaconda3/bin/activate
conda activate "${CONDA_ENV}"
set -u
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8

RUNS="${BASE}/04_saturn_runs"
OF6="${BASE}/06_orthofinder_6sp"
SP6="Mlei,Cgig,Crob,Drer,Dmel,Hsap"
RUN_E="${RUNS}/20261019_Mlei_Cgig_Crob_Drer_Dmel_Hsap_K3000_ESM1b_macrogene"
RUN_C="${RUNS}/20261019_Mlei_Cgig_Crob_Drer_Dmel_Hsap_K3000_ESMCpostnorm_macrogene"
TAG="$(date +%Y%m%d)"
has_pkl() { compgen -G "$1/saturn_results/*genes_to_macrogenes*.pkl" > /dev/null; }

has_pkl "${RUN_E}" || { echo "ESM-1b 6-sp run not finished: no pkl in ${RUN_E}/saturn_results" >&2; exit 1; }
compgen -G "${OF6}/Results_*/*/Orthogroups/Orthogroups.tsv" > /dev/null || \
compgen -G "${OF6}/Results_*/Orthogroups/Orthogroups.tsv" > /dev/null || \
    { echo "6-sp OrthoFinder not finished: no Orthogroups.tsv under ${OF6}" >&2; exit 1; }

RUNLIST="ESM1b_6sp=${RUN_E}"
HV_B="${RUN_E}"
if has_pkl "${RUN_C}"; then
    RUNLIST="${RUNLIST},ESMCpostnorm_6sp=${RUN_C}"
    HV_B="${RUN_C}"
else
    echo "NOTE: no ESMC postnorm 6-sp run -- scoring ESM-1b alone in 05g"
fi
echo "=== 06e start $(date) on $(hostname) ==="

echo; echo "--- 1. 06c OrthoFinder -> reference tables ---"
python "${SCRIPT_DIR}/06c_orthofinder_to_reference.py" --of_dir "${OF6}" --saturn_run "${RUN_E}"

echo; echo "--- 2a. 05g vs OrthoFinder OGs ---"
python "${SCRIPT_DIR}/05g_compare_runs_orthology.py" --species "${SP6}" --runs "${RUNLIST}" \
    --reference "${OF6}/ALL6_gene_orthofinder_OG.csv" --out_dir "${RUNS}/${TAG}_6sp_benchmark_OG"

HOG=$(ls "${OF6}"/ALL6_gene_orthofinder_N*HOG.csv 2>/dev/null | head -n 1 || true)
if [[ -n "${HOG}" ]]; then
    echo; echo "--- 2b. 05g vs OrthoFinder Metazoa HOGs ($(basename "${HOG}")) ---"
    python "${SCRIPT_DIR}/05g_compare_runs_orthology.py" --species "${SP6}" --runs "${RUNLIST}" \
        --reference "${HOG}" --out_dir "${RUNS}/${TAG}_6sp_benchmark_HOG"
else
    echo "No HOG table written by 06c; skipping 2b"
fi

echo; echo "--- 3. 05i nearest-neighbour ortholog test, ESM-1b vs ESMC postnorm ---"
python "${SCRIPT_DIR}/05i_embedding_standardization_check.py" --species "${SP6}" \
    --hv_run_a "${RUN_E}" --hv_run_b "${HV_B}" \
    --esmc_dir 03_embeddings_ESMC600M_postnorm --esmc_label ESMC600M_postnorm \
    --reference "${OF6}/ALL6_gene_orthofinder_OG.csv" \
    --out_dir "${RUNS}/${TAG}_6sp_standardization_OG"

echo; echo "=== 06e done $(date) ==="
echo "Read: ${OF6}/06c_summary.txt"
echo "      ${RUNS}/${TAG}_6sp_benchmark_OG/summary.txt  (and _HOG)"
echo "      ${RUNS}/${TAG}_6sp_standardization_OG/summary.txt"
