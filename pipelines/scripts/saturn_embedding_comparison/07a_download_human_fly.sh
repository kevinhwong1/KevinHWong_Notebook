#!/usr/bin/env bash
# 07a_download_human_fly.sh
# Downloads the two new species for the pLM benchmark (see DATASET_SCOUTING_human_fly_20261014.md):
#   Hsap  CellTypist cross-tissue immune atlas (Dominguez Conde et al. 2022, Science), Sanger "Download Raw"
#         global object (329,762 cells; raw counts), + Ensembl 115 GRCh38 proteome and GTF
#   Dmel  larval hemocytes: Tattikota et al. 2020 eLife, GEO GSE146596 (whole supplementary folder + SOFT metadata)
#         adult hemocytes:  Fly Cell Atlas "10x Hemocyte (Stringent)" h5ad (snRNA-seq) -- set FCA_HEMO_URL below
#         + Ensembl 115 BDGP6 proteome and GTF (gene IDs = FlyBase FBgn)
# Runs on the 'transfers' queue (login-node jobs over 30 min are killed). Resumable: wget -c, re-submit if cut off.
# Data are public and land on /scratch (purged after 21 days): re-run this script to restore them.
# Submit:  bsub < /scratch/dark_genes/SATURN_Mnemi/scripts/07a_download_human_fly.sh
#BSUB -J dl_hsap_dmel
#BSUB -q transfers
#BSUB -P dark_genes
#BSUB -n 1
#BSUB -W 12:00
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/dl_hsap_dmel.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/dl_hsap_dmel.%J.err

set -uo pipefail   # no -e: one failed download should not stop the others; failures are listed at the end

# Fly Cell Atlas hemocyte h5ad: on https://flycellatlas.org, row "10x Hemocyte (Stringent)", right-click
# "H5AD" -> Copy link. Looks like https://cloud.flycellatlas.org/index.php/s/<token>/download
FCA_HEMO_URL="${FCA_HEMO_URL:-}"

BASE="/scratch/dark_genes/SATURN_Mnemi"
RAW="${BASE}/raw_data"
PROT="${BASE}/01_proteomes"
ENS="ftp://ftp.ensembl.org/pub/release-115"
mkdir -p "${RAW}/Hsap" "${RAW}/Dmel_larva_GSE146596" "${RAW}/Dmel_adult_FCA" "${PROT}/Hsap" "${PROT}/Dmel" "${BASE}/logs"
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
get_pattern() {  # get_pattern <ftp_dir> <pattern> <dest_dir>   (Ensembl file names include the assembly version)
    echo "--- $(date '+%F %T')  ${1}  [${2}]"
    wget -c -nv -r -l1 -nd --no-parent --tries=5 --timeout=60 -A "$2" -P "$3" "$1" || FAILED+=("$1 $2")
}

echo "=== start $(date) on $(hostname) ==="

# ---- Human
get "https://cellgeni.cog.sanger.ac.uk/pan-immune/CountAdded_PIP_global_object_for_cellxgene.h5ad" "${RAW}/Hsap"
get_pattern "${ENS}/fasta/homo_sapiens/pep/" "*.pep.all.fa.gz" "${PROT}/Hsap"
get_pattern "${ENS}/gtf/homo_sapiens/" "Homo_sapiens.GRCh38.115.gtf.gz" "${PROT}/Hsap"

# ---- Fly, larval (GEO): list the supplementary folder first, then mirror it
GEO="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146596"
echo "--- GSE146596 supplementary files:"
curl -s --max-time 120 "${GEO}/suppl/" || echo "(listing failed)"
wget -c -nv -r -nd --no-parent --tries=5 --timeout=60 -P "${RAW}/Dmel_larva_GSE146596" "${GEO}/suppl/" \
    || FAILED+=("${GEO}/suppl/")
get "${GEO}/soft/GSE146596_family.soft.gz" "${RAW}/Dmel_larva_GSE146596"

# ---- Fly, adult (Fly Cell Atlas)
if [[ -n "${FCA_HEMO_URL}" ]]; then
    get "${FCA_HEMO_URL}" "${RAW}/Dmel_adult_FCA" "FCA_hemocyte_stringent_10x.h5ad"
else
    echo "--- FCA_HEMO_URL not set: skipping the Fly Cell Atlas hemocyte h5ad"; FAILED+=("FCA_HEMO_URL not set")
fi

# ---- Fly proteome + GTF (Ensembl; gene IDs are FBgn)
get_pattern "${ENS}/fasta/drosophila_melanogaster/pep/" "*.pep.all.fa.gz" "${PROT}/Dmel"
get_pattern "${ENS}/gtf/drosophila_melanogaster/" "Drosophila_melanogaster.*.115.gtf.gz" "${PROT}/Dmel"

# ---- Check what arrived
echo; echo "=== Files ==="
find "${RAW}/Hsap" "${RAW}/Dmel_larva_GSE146596" "${RAW}/Dmel_adult_FCA" "${PROT}/Hsap" "${PROT}/Dmel" -type f -exec ls -lh {} \;
echo; echo "=== gzip integrity ==="
for f in $(find "${RAW}" "${PROT}/Hsap" "${PROT}/Dmel" -name "*.gz" | sort); do
    if gzip -t "$f" 2>/dev/null; then echo "OK   $f"; else echo "BAD  $f"; FAILED+=("corrupt: $f"); fi
done
echo; echo "=== h5ad headers (should start with HDF5 signature) ==="
for f in $(find "${RAW}" -name "*.h5ad"); do
    head -c 8 "$f" | od -c | head -1 | grep -q "H   D   F" && echo "OK   $f" || { echo "BAD  $f (not HDF5; check the URL)"; FAILED+=("not hdf5: $f"); }
done

echo
if (( ${#FAILED[@]} )); then
    echo "=== ${#FAILED[@]} problem(s) ==="; printf '  %s\n' "${FAILED[@]}"
else
    echo "=== all downloads OK ==="
fi
echo "=== end $(date) ==="
