#!/usr/bin/env bash
# 08c_macrogene_2x2.sh   (2026-10-30)
# Macrogene composition across the 2x2 (ESM-1b / ESMC post-norm x averaged / longest isoform):
#   1. 05g  run descriptives (R1), pairwise ARI/AMI (R2), ortholog benchmark (R3) on the 4 runs, vs OrthoFinder OG
#   2. 05g  same vs the Metazoa HOG table
#   3. 08c  isoform / model effects + interaction (C1), per-gene partner stability (C2), single-species MGs (C3)
# SET=prelim (default, 10-30): averaged conditions = the 10-19 macrogene-only runs. ESM1b_avg has the same inputs
#   as full run 1; ESMC_avg uses the OLD unmatched averaged ESMC files (Drer gene set differs; 05g/08c use the genes
#   common to all runs). Longest = full runs 3 and 4 (614411/614412).
# SET=final: averaged conditions = full runs 1 and 2 (05y, matched ESMC). Use once 622420/622421 have finished.
# Submit (from /scratch/dark_genes/SATURN_Mnemi):  bsub < scripts/08c_macrogene_2x2.sh      (or SET=final bsub < ...)
#BSUB -J mg_2x2
#BSUB -q normal
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 06:00
#BSUB -R "rusage[mem=24000] span[hosts=1]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/mg_2x2.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/mg_2x2.%J.err

set -euo pipefail
SCRIPT_DIR="/scratch/dark_genes/SATURN_Mnemi/scripts"
source "${SCRIPT_DIR}/run_config.sh"
set +u
source ~/anaconda3/bin/activate
conda activate "${CONDA_ENV}"
set -u
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
cd "${BASE}"

R="${BASE}/04_saturn_runs"
OF6="${BASE}/06_orthofinder_6sp"
SP6="Mlei,Cgig,Crob,Drer,Dmel,Hsap"
P6="Mlei_Cgig_Crob_Drer_Dmel_Hsap_K3000"
SET="${SET:-prelim}"
case "${SET}" in
    prelim)
        RUN_E_AVG="${R}/20261019_${P6}_ESM1b_macrogene"
        RUN_C_AVG="${R}/20261019_${P6}_ESMCpostnorm_macrogene"
        EMB_C_AVG="${BASE}/03_embeddings_ESMC600M_postnorm" ;;
    final)
        RUN_E_AVG="${R}/$(cat "${R}/LATEST_FULL_6sp_ESM1b_avg_seed42.txt")"
        RUN_C_AVG="${R}/$(cat "${R}/LATEST_FULL_6sp_ESMCpostnorm_avg_seed42.txt")"
        EMB_C_AVG="${BASE}/03_embeddings_ESMC600M_postnorm_avgmatched" ;;
    *) echo "SET=prelim|final" >&2; exit 1 ;;
esac
RUN_E_LONG="${R}/20261027_${P6}_ESM1b_longest_full"
RUN_C_LONG="${R}/20261027_${P6}_ESMCpostnorm_longest_full"
for d in "${RUN_E_AVG}" "${RUN_C_AVG}" "${RUN_E_LONG}" "${RUN_C_LONG}"; do
    compgen -G "${d}/saturn_results/*genes_to_macrogenes*.pkl" > /dev/null || { echo "no macrogene pkl in ${d}" >&2; exit 1; }
done
RUNLIST="ESM1b_avg=${RUN_E_AVG},ESMC_avg=${RUN_C_AVG},ESM1b_longest=${RUN_E_LONG},ESMC_longest=${RUN_C_LONG}"
EMBLIST="ESM1b_avg=${BASE}/03_embeddings,ESMC_avg=${EMB_C_AVG},ESM1b_longest=${BASE}/03_embeddings_ESM1b_longest,ESMC_longest=${BASE}/03_embeddings_ESMC600M_postnorm_longest"
OUT="${R}/$(date +%Y%m%d)_macrogene_2x2_${SET}"
mkdir -p "${OUT}"
echo "=== 08c start $(date) on $(hostname)  SET=${SET} ==="
echo "${RUNLIST}" | tr ',' '\n'

echo; echo "--- 1. 05g vs OrthoFinder OG ---"
python "${SCRIPT_DIR}/05g_compare_runs_orthology.py" --species "${SP6}" --runs "${RUNLIST}" \
    --reference "${OF6}/ALL6_gene_orthofinder_OG.csv" --out_dir "${OUT}/05g_OG"

HOG=$(ls "${OF6}"/ALL6_gene_orthofinder_N*HOG.csv 2>/dev/null | head -n 1 || true)
if [[ -n "${HOG}" ]]; then
    echo; echo "--- 2. 05g vs Metazoa HOG ($(basename "${HOG}")) ---"
    python "${SCRIPT_DIR}/05g_compare_runs_orthology.py" --species "${SP6}" --runs "${RUNLIST}" \
        --reference "${HOG}" --out_dir "${OUT}/05g_HOG"
fi

echo; echo "--- 3. 08c 2x2 effects, partner stability, single-species macrogenes ---"
python "${SCRIPT_DIR}/08c_macrogene_2x2.py" --runs "${RUNLIST}" --emb "${EMBLIST}" \
    --og_dir "${OUT}/05g_OG" --hog_dir "${OUT}/05g_HOG" --out_dir "${OUT}/08c"

echo; echo "=== 08c done $(date) ==="
echo "Download: ${OUT}/   (05g_OG/, 05g_HOG/, 08c/)"
