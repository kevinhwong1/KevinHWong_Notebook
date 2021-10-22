---
layout: post
title: HIAST2 alignment with Mansour transcriptome
date: '2021-10-22'
categories: Analysis
tags: Porites astreoides, tagseq, mansour
---

Goal: aligning TagSeq data to the Mansour Porites astreoides reference transcriptome

## script

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH -D /data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/hisat2/Mansour_REF
#SBATCH --cpus-per-task=3


#load packages
module load HISAT2/2.1.0-foss-2018b #Alignment to reference genome: HISAT2
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

# symbolically link 'clean' reads to hisat2 dir
ln -s /data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/clean/clean*.fastq.gz ./

# index the reference transcriptome for Porites astreoides Mansour Transcriptome output index to working directory
hisat2-build -f /data/putnamlab/kevin_wong1/REFS/Past_Mansour/p_ast2016.fasta ./Past_trans_Mansour_ref # called the reference genome (scaffolds)
echo "Reference genome indexed. Starting alignment" $(date)

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed
array=($(ls *.fastq.gz)) # call the symbolically linked sequences - make an array to align
for i in ${array[@]}; do
        sample_name=`echo $i| awk -F [.] '{print $2}'`
        hisat2 -p 8 --dta -x Past_trans_Mansour_ref -U ${i} -S ${sample_name}.sam
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
                echo "${i} bam-ified!"
        rm ${sample_name}.sam
done
```

## outputs

Very low alignment rates (3-5%)
