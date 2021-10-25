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
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

# load modules needed

module load Nextflow

# run nextflow methylseq

nextflow run nf-core/methylseq -profile singularity \
--aligner bismark \
--igenomes_ignore \
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


-------

I am having errors with versions of Nextflow and the [igenome conflicts](https://nf-co.re/methylseq/1.5/usage) so I am changing the following:

- adding the `--igenomes_ignore` flag
- updating Nextflow and using methylseq v1.6.1
- used the `--input` flag instead of `--reads`

```
#!/bin/bash
#SBATCH --job-name="methylseq"
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

# load modules needed

module load Nextflow/21.03.0

# run nextflow methylseq

nextflow run nf-core/methylseq -profile singularity \
--aligner bismark \
--igenomes_ignore \
--fasta /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--save_reference \
--input '/data/putnamlab/KITT/hputnam/20211008_Past_ThermalTransplant_WGBS/*_R{1,2}_001.fastq.gz' \
--clip_r1 10 \
--clip_r2 10 \
--three_prime_clip_r1 10 --three_prime_clip_r2 10 \
--non_directional \
--cytosine_report \
--relax_mismatches \
--unmapped \
--outdir WGBS_methylseq \
-name WGBS_methylseq_past_TT_1
```

### Running methylseq script

`sbatch /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq.sh`

- Started October 21 2021 @ 15:00
- Job ID 93817

## Resuming methylseq script because it timed out

`nano methylseq_resume.sh`

```
#!/bin/bash
#SBATCH --job-name="methylseq"
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS
#SBATCH -p putnamlab
#SBATCH --exclusive
#SBATCH --cpus-per-task=3

# load modules needed

module load Nextflow/21.03.0

# run nextflow methylseq

nextflow run nf-core/methylseq -resume \
-profile singularity \
--aligner bismark \
--igenomes_ignore \
--fasta /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--save_reference \
--input '/data/putnamlab/KITT/hputnam/20211008_Past_ThermalTransplant_WGBS/*_R{1,2}_001.fastq.gz' \
--clip_r1 10 \
--clip_r2 10 \
--three_prime_clip_r1 10 --three_prime_clip_r2 10 \
--non_directional \
--cytosine_report \
--relax_mismatches \
--unmapped \
--outdir WGBS_methylseq
```

### Running methylseq script

`sbatch /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_resume.sh`

- Started October 25 2021 @ 10:54
- Job ID 94131
