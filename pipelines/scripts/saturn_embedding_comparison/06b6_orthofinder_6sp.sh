#!/usr/bin/env bash
# 06b6_orthofinder_6sp.sh  (2026-10-19; copy of 06b_orthofinder_4sp.sh for the 6-species panel)
# OrthoFinder v3.1.0 on the 6 SATURN species (Mlei, Cgig, Crob, Drer, Dmel, Hsap) + outgroup Sros
# (choanoflagellate Salpingoeca rosetta; 06a0_download_outgroup.sh), longest protein per gene, from:
#   python scripts/06a_orthofinder_prep.py --species Mlei,Cgig,Crob,Drer,Dmel,Hsap,Sros \
#          --out_dir /scratch/dark_genes/SATURN_Mnemi/06_orthofinder_6sp
# De novo orthology reference for the 6-species ESM-1b vs ESMC benchmarks (06c -> 05g / 05p / 06d).
# Differences from 06b: input/output folder 06_orthofinder_6sp, 7 proteomes, mem 32 GB (4sp peaked at 15.9 GB
# with 91k proteins; this run has ~1.5x the proteins), tree from 06a with the outgroup at the root.
# --- original 06b notes below ---
# Flags as in the coral pipeline (submit_pipeline.sh, job 02): -M msa -A mafft -T fasttree.
# 2026-10-15: conda envs on Pegasus were wiped; orthofinder_env was rebuilt with micromamba (lean install,
#   see notebook 2026-10-15). Activation now uses micromamba by full path (batch jobs don't read .bashrc).
#   Queue: general is closed; bigmem has a long backlog, so this runs on normal. The old gpu1/gpu2 host pin
#   (FastTree/FAMSA AVX2 crash on older nodes) is replaced by a FastTree check at the start of the job:
#   if the node can't run FastTree the job stops in seconds instead of failing hours in. If that happens,
#   resubmit with:  -q bigmem -R "rusage[mem=15000] select[hname==gpu1 || hname==gpu2]".
# Fixed species tree (-s): (Sros,(Mlei,((Cgig,Dmel),(Crob,(Drer,Hsap))))); N0 = Choanozoa (root);
# Metazoa-level HOGs = the node joining Mlei and the bilaterians (expected N1; check SpeciesTree_rooted_node_labels.txt).
# Submit (after 06a0 + 06a):  bsub < /scratch/dark_genes/SATURN_Mnemi/scripts/06b6_orthofinder_6sp.sh
#BSUB -J of_saturn6
#BSUB -q normal
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 72:00
#BSUB -R "rusage[mem=32000] span[hosts=1]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/of_saturn6.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/of_saturn6.%J.err

set -euo pipefail
OF_DIR="/scratch/dark_genes/SATURN_Mnemi/06_orthofinder_6sp"
OUT="${OF_DIR}/Results_$(date +%Y%m%d)"

for sp in Mlei Cgig Crob Drer Dmel Hsap Sros; do
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

echo "=== OrthoFinder, 6 SATURN species + Sros outgroup ===" && cat "${OF_DIR}/species_tree.nwk" && date && hostname
orthofinder --version || true
orthofinder \
    -f "${OF_DIR}/input" \
    -o "${OUT}" \
    -s "${OF_DIR}/species_tree.nwk" \
    -t 16 -a 16 -M msa -A mafft -T fasttree
echo "=== OrthoFinder complete ===" && date
find "${OUT}" -name "Orthogroups.tsv" -o -name "N*.tsv" -o -name "SpeciesTree_rooted_node_labels.txt" | sort | head -20
find "${OUT}" -name "SpeciesTree_rooted_node_labels.txt" -exec cat {} \;
