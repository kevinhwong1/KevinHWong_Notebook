---
layout: post
title: BisSNP Analysis
date: '2023-08-04'
categories: Analysis
tags: DNA methylation, porites astreoides
---

# Bis-SNP analysis 

## Introduction

Bis-SNP is used to extract SNPs from whole genome bisulfite seq data. This markdown file will follow the [user guide](https://people.csail.mit.edu/dnaase/bissnp2011/BisSNP-UserGuide-latest.pdf) and [website](https://people.csail.mit.edu/dnaase/bissnp2011/). 

Bis-SNP is the public available free software (GPL v3 license) for genotyping in bisulfite treated massively parallel
sequencing (whole genome Bisulfite-seq(BS-seq), NOMe-seq and RRBS) on Illumina platform. It works for both of
single-end and paired-end reads in Illumina directional Bisulfite-Seq library. It is implemented in Java and based on
GATK map-reduce framework for the parallel computation. Copyright belongs to USC Epigenome Center.

![](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/BisSNP_diagram.png?raw=true)

## Data

I will be using whole genome bisulfte data from the Thermal Transplant Moleclar project looking at different adult/larval thermal histories. 

* 47 samples in total

input data path: `/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated/`

## Modules

Test to see if bis-snp is installed on Andromeda using `module av -t |& grep -i Bis` :

```bash
bis-snp/1.0.1.3-Java-13
```

It should be available to use. 

## Mark duplicated reads

Use Picard tools to mark the duplicated reads which are mostly come from PCR duplication. This step could be
done before indel realignment if there are too many duplicated reads which would mislead indel alignment:


java -Xmx10g -jar MarkDuplicates.jar I=sample.withRG.realigned.bam O=sample.withRG.realigned.mdups.bam
METRICS_FILE=sample.withRG.realigned.metric.txt CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT

`nano episnp.sh`

```bash
#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=18
#SBATCH --export=NONE
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --error=output_messages/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output_messages/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed (specific need for my computer)
source /usr/share/Modules/init/sh # load the module function

# load modules needed
echo "START" $(date)
module load Nextflow/20.07.1 #this pipeline requires this version 
module load SAMtools/1.9-foss-2018b 
Pysam/0.15.1-foss-2018b-Python-3.6.6

# define location for fasta_generate_regions.py
fasta_generate_regions.py = /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/fasta_generate_regions.py

# only need to direct to input folder not *bam files 
NXF_VER=20.07.1 nextflow run epidiverse/snp -resume \
--input /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated/ \
--reference /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--output /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/ \
--clusters \
--variants \
--coverage 5 \
--take 47 # Number of samples 

echo "STOP" $(date) # this will output the time it takes to run within the output message

```