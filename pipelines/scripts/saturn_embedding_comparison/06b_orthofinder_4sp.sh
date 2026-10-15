#!/usr/bin/env bash
# 06b_orthofinder_4sp.sh
# OrthoFinder v3.1.0 on the 4 SATURN species (longest protein per gene, from 06a), as a de novo orthology
# reference to compare with eggNOG Metazoa OGs in 05n / 05p.
# Flags as in the coral pipeline (submit_pipeline.sh, job 02): -M msa -A mafft -T fasttree.
# 2026-10-15: conda envs on Pegasus were wiped; orthofinder_env was rebuilt with micromamba (lean install,
#   see notebook 2026-10-15). Activation now uses micromamba by full path (batch jobs don't read .bashrc).
#   Queue: general is closed; bigmem has a long backlog, so this runs on normal. The old gpu1/gpu2 host pin
#   (FastTree/FAMSA AVX2 crash on older nodes) is replaced by a FastTree check at the start of the job:
#   if the node can't run FastTree the job stops in seconds instead of failing hours in. If that happens,
#   resubmit with:  -q bigmem -R "rusage[mem=15000] select[hname==gpu1 || hname==gpu2]".
# Fixed species tree (-s): (Mlei,(Cgig,(Crob,Drer))); with these 4 species the root (ctenophore vs bilaterians)
# is the same under ctenophore-first and sponge-first hypotheses, so N0 = the Metazoa-level HOGs.
# Submit (after 06a):  bsub < /scratch/dark_genes/SATURN_Mnemi/scripts/06b_orthofinder_4sp.sh
#BSUB -J of_saturn4
#BSUB -q normal
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 72:00
#BSUB -R "rusage[mem=15000] span[hosts=1]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/of_saturn4.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/of_saturn4.%J.err

set -euo pipefail
OF_DIR="/scratch/dark_genes/SATURN_Mnemi/06_orthofinder"
OUT="${OF_DIR}/Results_20261015"

for sp in Mlei Cgig Crob Drer; do
    [[ -s "${OF_DIR}/input/${sp}.faa" ]] || { echo "Missing ${OF_DIR}/input/${sp}.faa -- run 06a first" >&2; exit 1; }
done
[[ -s "${OF_DIR}/species_tree.nwk" ]] || { echo "Missing species_tree.nwk -- run 06a first" >&2; exit 1; }
if [[ -e "${OUT}" ]]; then echo "${OUT} exists -- move it or change OUT" >&2; exit 1; fi

set +u
# micromamba env (full paths: batch jobs don't read ~/.bashrc). Override env name with OF_ENV.
OF_ENV="${OF_ENV:-orthofinder_env}"
export MAMBA_ROOT_PREFIX=/nethome/kxw755/micromamba
eval "$(/nethome/kxw755/.local/bin/micromamba shell hook --shell bash)"
micromamba activate "${OF_ENV}"
set -u
echo "Using micromamba env: ${OF_ENV}"
which orthofinder diamond mcl mafft FastTree

# Preflight: FastTree on this node (old nodes without AVX2 crash with 'Illegal instruction')
PRE=$(mktemp -d)
printf ">a\nACDEFGHIKLMNPQ\n>b\nACDEFGHIKLMNPR\n>c\nACDQFGHIKLMNPQ\n>d\nACDEFGHIRLMNPQ\n" > "${PRE}/t.fa"
if ! FastTree -quiet "${PRE}/t.fa" > "${PRE}/t.nwk" 2> "${PRE}/t.err"; then
    echo "FastTree preflight FAILED on $(hostname); see message below. Resubmit pinned to gpu1/gpu2 (see header)." >&2
    cat "${PRE}/t.err" >&2; exit 1
fi
echo "FastTree preflight OK on $(hostname)"; rm -rf "${PRE}"
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
