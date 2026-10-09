#!/usr/bin/env bash
# 06b_orthofinder_4sp.sh
# OrthoFinder v3.1.0 on the 4 SATURN species (longest protein per gene, from 06a), as a de novo orthology
# reference to compare with eggNOG Metazoa OGs in 05n / 05p.
# Same environment and flags as the working coral pipeline (submit_pipeline.sh, job 02):
#   conda env found automatically (name/path containing "ortho", or any env with orthofinder; override with OF_ENV); -M msa -A mafft -T fasttree; gpu1/gpu2 hosts (FastTree/FAMSA AVX2 issue on older nodes).
# Fixed species tree (-s): (Mlei,(Cgig,(Crob,Drer))); with these 4 species the root (ctenophore vs bilaterians)
# is the same under ctenophore-first and sponge-first hypotheses, so N0 = the Metazoa-level HOGs.
# Submit (after 06a):  bsub < /scratch/dark_genes/SATURN_Mnemi/scripts/06b_orthofinder_4sp.sh
#BSUB -J of_saturn4
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 72:00
#BSUB -R "rusage[mem=15000] select[hname==gpu1 || hname==gpu2]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/of_saturn4.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/of_saturn4.%J.err

set -euo pipefail
OF_DIR="/scratch/dark_genes/SATURN_Mnemi/06_orthofinder"
OUT="${OF_DIR}/Results_20261009"

for sp in Mlei Cgig Crob Drer; do
    [[ -s "${OF_DIR}/input/${sp}.faa" ]] || { echo "Missing ${OF_DIR}/input/${sp}.faa -- run 06a first" >&2; exit 1; }
done
[[ -s "${OF_DIR}/species_tree.nwk" ]] || { echo "Missing species_tree.nwk -- run 06a first" >&2; exit 1; }
if [[ -e "${OUT}" ]]; then echo "${OUT} exists -- move it or change OUT" >&2; exit 1; fi

set +u
source ~/anaconda3/etc/profile.d/conda.sh
# Environment: OF_ENV if set (bsub -env "all, OF_ENV=name"), else the first conda env whose name or path
# contains "ortho" (case-insensitive), else any env that has an orthofinder executable.
if [[ -z "${OF_ENV:-}" ]]; then
    OF_ENV=$(conda env list | awk 'tolower($0) ~ /ortho/ && $1 !~ /^#/ {print $NF; exit}')
fi
if [[ -z "${OF_ENV:-}" ]]; then
    for p in $(conda env list | awk '$1 !~ /^#/ {print $NF}'); do
        if [[ -x "${p}/bin/orthofinder" || -x "${p}/bin/orthofinder.py" ]]; then OF_ENV="${p}"; break; fi
    done
fi
[[ -n "${OF_ENV:-}" ]] || { echo "No conda env with OrthoFinder found; set OF_ENV" >&2; conda env list >&2; exit 1; }
echo "Using conda env: ${OF_ENV}"
conda activate "${OF_ENV}"
set -u
command -v orthofinder >/dev/null || { echo "orthofinder not on PATH in ${OF_ENV}" >&2; exit 1; }

echo "=== OrthoFinder, 4 SATURN species ===" && date && hostname
orthofinder --version || true
orthofinder \
    -f "${OF_DIR}/input" \
    -o "${OUT}" \
    -s "${OF_DIR}/species_tree.nwk" \
    -t 16 -a 16 -M msa -A mafft -T fasttree
echo "=== OrthoFinder complete ===" && date
find "${OUT}" -name "Orthogroups.tsv" -o -name "N0.tsv" | head
