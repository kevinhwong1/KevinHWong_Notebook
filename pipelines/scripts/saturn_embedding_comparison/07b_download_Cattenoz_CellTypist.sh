#!/usr/bin/env bash
# 07b_download_Cattenoz_CellTypist.sh
# Downloads the FINAL two new species for the pLM benchmark (decision 2026-10-15, see
# claude/DATASET_SCOUTING_human_fly_20261014.md). Supersedes the data part of 07a
# (Tattikota + Fly Cell Atlas are no longer used).
#   Dmel  Cattenoz et al. 2020 EMBO J, ArrayExpress E-MTAB-8698, "NI dataset" (not-infested
#         wandering L3, OregonR, 10x 3'), library RRCZ22: Cell Ranger filtered + raw matrices,
#         author cluster labels (NI_cell_cluster_ID.txt), md5, SDRF/IDF.
#   Hsap  CellTypist cross-tissue atlas (Dominguez Conde et al. 2022), CELLxGENE collection
#         62ef75e4-cbea-454e-a0ce-998ec40223d3: "Myeloid compartment" (51,552 cells, 0.7 GB) and
#         "Global" (329,762 cells, 3 GB; for the lymphocyte slice). Raw counts in raw.X.
# Proteomes/GTFs are NOT downloaded here (Ensembl release still to be decided: 113 to match Drer, or 115 as in 07a).
# Runs on the 'transfers' queue (login-node jobs over 30 min are killed). Resumable: wget -c.
# Submit:  bsub < /scratch/dark_genes/SATURN_Mnemi/scripts/07b_download_Cattenoz_CellTypist.sh
#BSUB -J dl_cattenoz_celltypist
#BSUB -q transfers
#BSUB -P dark_genes
#BSUB -n 1
#BSUB -W 6:00
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/dl_cattenoz_celltypist.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/dl_cattenoz_celltypist.%J.err

set -uo pipefail   # no -e: one failed download should not stop the others; failures are listed at the end

BASE="/scratch/dark_genes/SATURN_Mnemi"
DMEL="${BASE}/raw_data/Dmel_Cattenoz_EMTAB8698"
HSAP="${BASE}/raw_data/Hsap_CellTypist_CXG"
mkdir -p "${DMEL}" "${HSAP}" "${BASE}/logs"
FAILED=()

get() {  # get <url> <dest_dir> [output_name]
    local url="$1" dir="$2" name="${3:-}"
    echo "--- $(date '+%F %T')  ${url}"
    if [[ -n "${name}" ]]; then
        wget -c -nv --tries=5 --timeout=60 -O "${dir}/${name}" "${url}" || FAILED+=("${url}")
    else
        wget -c -nv --tries=5 --timeout=60 -P "${dir}" "${url}" || FAILED+=("${url}")
    fi
}

echo "=== start $(date) on $(hostname) ==="

# ---- Fly: Cattenoz NI (RRCZ22)
AE="https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/698/E-MTAB-8698/Files"
for f in filtered_feature_bc_matrix_RRCZ22.tgz raw_feature_bc_matrix_RRCZ22.tgz \
         NI_cell_cluster_ID.txt md5_NI.txt E-MTAB-8698.sdrf.txt E-MTAB-8698.idf.txt; do
    get "${AE}/${f}" "${DMEL}"
done
mkdir -p "${DMEL}/RRCZ22_filtered" "${DMEL}/RRCZ22_raw"
tar -xzf "${DMEL}/filtered_feature_bc_matrix_RRCZ22.tgz" -C "${DMEL}/RRCZ22_filtered" || FAILED+=("untar filtered")
tar -xzf "${DMEL}/raw_feature_bc_matrix_RRCZ22.tgz"      -C "${DMEL}/RRCZ22_raw"      || FAILED+=("untar raw")

# ---- Human: CellTypist via CELLxGENE
CXG="https://datasets.cellxgene.cziscience.com"
get "${CXG}/96de0195-492e-4180-b0b8-49cc04b7293e.h5ad" "${HSAP}" "CellTypist_Myeloid.h5ad"
get "${CXG}/09a381ce-feaf-48c8-a928-58f6768f4265.h5ad" "${HSAP}" "CellTypist_Global.h5ad"

# ---- Check what arrived
echo; echo "=== Files ==="
find "${DMEL}" "${HSAP}" -type f -exec ls -lh {} \;
echo; echo "=== md5 (compare by eye with md5_NI.txt) ==="
cat "${DMEL}/md5_NI.txt"
md5sum "${DMEL}/filtered_feature_bc_matrix_RRCZ22.tgz" "${DMEL}/raw_feature_bc_matrix_RRCZ22.tgz"
echo; echo "=== h5ad headers (should start with HDF5 signature) ==="
for f in "${HSAP}"/*.h5ad; do
    head -c 8 "$f" | od -c | head -1 | grep -q "H   D   F" && echo "OK   $f" || { echo "BAD  $f (not HDF5)"; FAILED+=("not hdf5: $f"); }
done

echo
if (( ${#FAILED[@]} )); then
    echo "=== ${#FAILED[@]} problem(s) ==="; printf '  %s\n' "${FAILED[@]}"
else
    echo "=== all downloads OK ==="
fi
echo "=== end $(date) ==="
