#!/usr/bin/env bash
# 05r_esmc_postnorm_check_4sp.sh   (2026-10-19)
# Phase 1B of PLAN_6species_OrthoFinder_ESM1b_vs_ESMCpostnorm_20261019: is Ingrid's second ESMC-600M export
# ("postnorm_maxlen2000") a real contender? Embedding-only test, no SATURN run. 4 species (Mlei Cgig Crob Drer),
# same HV gene set as the Oct 2-9 comparison, so the numbers sit next to the earlier ones:
#   ESM-1b hit@1 (euclidean, mean of 12 pairs) 0.53 raw / 0.61 abtt2-3;  old ESMC v4 0.34 raw / 0.49 abtt3.
# Steps
#   1. 05q   convert postnorm .h5 -> 03_embeddings_ESMC600M_postnorm/{F}/{F}_gene_embeddings.pt
#   2. 05f   descriptives (vector length, geometry, neighbour mixing, eggNOG hit@1) for
#            ESM1b + ESMC v4 raw + ESMC v4 rescaled + ESMC postnorm
#   3. 05i   standardization check (raw/unit/center/zscore/abtt1-3) ESM1b vs ESMC postnorm,
#            a) against eggNOG Metazoa OGs, b) against the 4-species OrthoFinder OGs (06c)
# Prerequisite: Ingrid's folder copied to ${BASE}/raw_embeddings/ESMC600M_postnorm/ (see plan doc for the rsync).
# Submit:  bsub < /scratch/dark_genes/SATURN_Mnemi/scripts/05r_esmc_postnorm_check_4sp.sh
#BSUB -J esmc_pn4
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 06:00
#BSUB -R "rusage[mem=48000]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/esmc_pn4.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/esmc_pn4.%J.err

set -euo pipefail
SCRIPT_DIR="/scratch/dark_genes/SATURN_Mnemi/scripts"
source "${SCRIPT_DIR}/run_config.sh"
set +u
source ~/anaconda3/bin/activate
conda activate "${CONDA_ENV}"
set -u
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8

RUNS="${BASE}/04_saturn_runs"
RAW="${BASE}/raw_embeddings/ESMC600M_postnorm"
OF_OG="${BASE}/06_orthofinder/ALL4_gene_orthofinder_OG.csv"
TAG="20261019"

for f in "${RAW}/proteome_h5_id_mapping.csv" \
         "${RAW}/Mnemi/Mnemi-longestisoform_esmc-600m_postnorm_maxlen2000.h5" \
         "${RAW}/Cgig/Cgig_esmc-600m_postnorm_maxlen2000.h5" \
         "${RAW}/Crob/Crob_esmc-600m_postnorm_maxlen2000.h5" \
         "${RAW}/Drer/Drer_esmc-600m_postnorm_maxlen2000.h5" "${OF_OG}"; do
    [[ -s "$f" ]] || { echo "Missing $f" >&2; exit 1; }
done
echo "=== 05r start $(date) on $(hostname) ==="

echo; echo "--- 1. 05q convert (Mlei Cgig Crob Drer) ---"
python "${SCRIPT_DIR}/05q_convert_esmc_postnorm.py" --species Mlei Cgig Crob Drer --compare_esm1b

echo; echo "--- 2. 05f descriptives ---"
python "${SCRIPT_DIR}/05f_embedding_descriptives.py" \
    --extra_model ESMC600M_postnorm=03_embeddings_ESMC600M_postnorm \
    --out_dir "${RUNS}/${TAG}_embedding_descriptives_ESMCpostnorm"

echo; echo "--- 3a. 05i standardization, eggNOG reference ---"
python "${SCRIPT_DIR}/05i_embedding_standardization_check.py" \
    --esmc_dir 03_embeddings_ESMC600M_postnorm --esmc_label ESMC600M_postnorm \
    --out_dir "${RUNS}/${TAG}_standardization_ESMCpostnorm_eggnog"

echo; echo "--- 3b. 05i standardization, 4-species OrthoFinder reference ---"
python "${SCRIPT_DIR}/05i_embedding_standardization_check.py" \
    --esmc_dir 03_embeddings_ESMC600M_postnorm --esmc_label ESMC600M_postnorm \
    --reference "${OF_OG}" \
    --out_dir "${RUNS}/${TAG}_standardization_ESMCpostnorm_orthofinder"

echo; echo "=== 05r done $(date) ==="
echo "Results: ${RUNS}/${TAG}_embedding_descriptives_ESMCpostnorm/summary.txt"
echo "         ${RUNS}/${TAG}_standardization_ESMCpostnorm_{eggnog,orthofinder}/summary.txt"
