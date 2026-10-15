#!/usr/bin/env bash
# 07e_download_proteomes_Ens113.sh
# Ensembl release 113 proteomes + GTFs for the two new species (matches Drer, Ensembl r113 GRCz11).
#   Hsap  Homo_sapiens.GRCh38.pep.all.fa.gz + Homo_sapiens.GRCh38.113.gtf.gz   (gene IDs = ENSG, as in Hsap_saturn.h5ad)
#   Dmel  Drosophila_melanogaster.BDGP6.46.pep.all.fa.gz + Drosophila_melanogaster.BDGP6.46.113.gtf.gz (gene IDs = FBgn, as in Dmel_saturn.h5ad)
# Writes to 01_proteomes/{Hsap,Dmel}/ (same layout as the other species), plus unzipped copies.
# Runs on the 'transfers' queue. Resumable (wget -c).
# Submit:  bsub < /scratch/dark_genes/SATURN_Mnemi/scripts/07e_download_proteomes_Ens113.sh
#BSUB -J dl_ens113
#BSUB -q transfers
#BSUB -P dark_genes
#BSUB -n 1
#BSUB -W 3:00
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/dl_ens113.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/dl_ens113.%J.err

set -uo pipefail

BASE="/scratch/dark_genes/SATURN_Mnemi"
PROT="${BASE}/01_proteomes"
ENS="https://ftp.ensembl.org/pub/release-113"
mkdir -p "${PROT}/Hsap" "${PROT}/Dmel" "${BASE}/logs"
FAILED=()

get() {  # get <url> <dest_dir>   (explicit file URLs: Ensembl robots.txt blocks recursive wget)
    echo "--- $(date '+%F %T')  ${1}"
    wget -c -nv --tries=5 --timeout=60 -P "$2" "$1" || FAILED+=("$1")
}

echo "=== start $(date) on $(hostname) ==="
get "${ENS}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"                          "${PROT}/Hsap"
get "${ENS}/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz"                                 "${PROT}/Hsap"
get "${ENS}/fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.46.pep.all.fa.gz" "${PROT}/Dmel"
get "${ENS}/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.113.gtf.gz"         "${PROT}/Dmel"

echo; echo "=== gzip integrity + unzip copies ==="
for f in "${PROT}"/Hsap/*.gz "${PROT}"/Dmel/*.gz; do
    [ -e "$f" ] || continue
    if gzip -t "$f"; then
        echo "OK   $f"
        gunzip -c "$f" > "${f%.gz}"
    else
        echo "BAD  $f"; FAILED+=("corrupt: $f")
    fi
done

echo; echo "=== Files ==="
ls -lh "${PROT}/Hsap" "${PROT}/Dmel"
echo; echo "=== protein / gene counts in FASTAs ==="
for f in "${PROT}"/Hsap/*.pep.all.fa "${PROT}"/Dmel/*.pep.all.fa; do
    [ -e "$f" ] || continue
    nprot=$(grep -c ">" "$f")
    ngene=$(grep ">" "$f" | grep -o "gene:[^ ]*" | sort -u | wc -l)
    echo "$(basename "$f"): ${nprot} proteins, ${ngene} genes"
done

for f in "${PROT}/Hsap/Homo_sapiens.GRCh38.pep.all.fa" "${PROT}/Dmel/Drosophila_melanogaster.BDGP6.46.pep.all.fa"; do
    [ -s "$f" ] || FAILED+=("missing or empty: $f")
done

echo
if (( ${#FAILED[@]} )); then
    echo "=== ${#FAILED[@]} problem(s) ==="; printf '  %s\n' "${FAILED[@]}"
else
    echo "=== all downloads OK ==="
fi
echo "=== end $(date) ==="
