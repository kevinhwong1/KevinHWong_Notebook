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

## Sort files

We are using this [pipeline](https://github.com/lyijin/pdae_dna_meth/tree/master/genetic_contribution/bissnp). 


`nano sort_bam.sh`

```bash
#!/bin/bash
#SBATCH --job-name="sort_bam"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated

# load modules needed

module load SAMtools/1.9-foss-2018b
module load picard/2.25.1-Java-11

for f in *.deduplicated.bam
do
  STEM=$(basename "${f}" _L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam)
  samtools sort "${f}" \
  -o "${STEM}".deduplicated_sorted.bam 
  java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I="${STEM}".deduplicated_sorted.bam O="${STEM}".deduplicated_sorted_rg.bam LB=lib1 PL=illumina PU=unit1 SM="${f}"
done

samtools index *rg.bam
```


```bash
#!/bin/bash
#SBATCH --job-name="bssnp_geno"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated

# load modules needed

module load bis-snp/1.0.1.3-Java-13

for a in *rg.bam; 
do 
b=`echo ${a} | sed 's/.deduplicated_sorted_rg.bam/vfn/'` && java -Xmx10g -jar BisSNP-1.0.1.jar -T BisulfiteGenotyper -C CG,1 -I ${a} -R ../../../../Past_Genome/past_filtered_assembly.fasta -vfn1 ${b}1.vcf -vfn2 snp_vcfs/${b}2.vcf 
done

```

usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java

https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-

https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups