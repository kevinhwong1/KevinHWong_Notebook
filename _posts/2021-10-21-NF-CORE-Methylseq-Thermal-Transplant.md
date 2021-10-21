---
layout: post
title: NF-CORE Methylseq Thermal Transplant
date: '2021-10-21'
categories: Analysis
tags: Porites astreoides, WGBS, methylseq, bioinformatics
---

## 20211021

Initial analysis of the Porites astreoides WGBS data after first round of sequencing.
- Using scripts from the [Coral Methylome project](https://github.com/hputnam/Coral_Methylomes/blob/main/BioInf.md)


### Making a new directory

`mkdir Thermal_Transplant_WGBS`

### Location of WGBS files:

`/data/putnamlab/KITT/hputnam/20211008_Past_ThermalTransplant_WGBS`

### Location of reference Porites astreoides genome:

`/data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta`


### Creating methylseq script

`nano methylseq.sh`

```
#!/bin/bash
#SBATCH --job-name="methylseq"
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

# load modules needed

module load Nextflow

# run nextflow methylseq

nextflow run nf-core/methylseq -profile singularity \
--aligner bismark \
--fasta /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--save_reference \
--reads '/data/putnamlab/KITT/hputnam/20211008_Past_ThermalTransplant_WGBS/*_R{1,2}_001.fastq.gz' \
--clip_r1 10 \
--clip_r2 10 \
--three_prime_clip_r1 10 --three_prime_clip_r2 10 \
--non_directional \
--cytosine_report \
--relax_mismatches \
--unmapped \
--outdir WGBS_methylseq \
-name WGBS_methylseq
```

### Running methylseq script

`sbatch /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq.sh`

- Started October 21 2021 @ 13:50
