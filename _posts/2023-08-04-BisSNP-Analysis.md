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


Resources: 

- https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-
- https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups


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


# Create reference sequence directory

`cp ../../../../Past_Genome/past_filtered_assembly.fasta`
`cp ../../../../Past_Genome/past_filtered_assembly.fasta.fai`

`mv past_filtered_assembly.fasta.fai past_filtered_assembly.fa.fai`

```
interactive
module load SAMtools/1.9-foss-2018b

samtools faidx past_filtered_assembly.fa
```

`nano seq.dir.sh`

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
module load picard/2.25.1-Java-11

java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary -R past_filtered_assembly.fasta -O past_filtered_assembly.dict
```

# Run Bisulfite genotyper

`mkdir vfn`

`nano bssnp_geno.sh`

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

module load bis-snp/1.0.1.3-Java-8

for a in *rg.bam; 
do 
b=`echo ${a} | sed 's/.deduplicated_sorted_rg.bam/vfn/'` && bis-snp -Xmx10g -T BisulfiteGenotyper -C CG,1 -I ${a} -R past_filtered_assembly.fa -vfn1 ${b}1.vcf -vfn2 snp_vcfs/${b}2.vcf 
done
```

Error message:
```bash
##### ERROR MESSAGE: SAM/BAM/CRAM file 18-118_S162.deduplicated_sorted_rg.bam is malformed. Please see https://software.broadinstitute.org/gatk/documentation/article?id=1317for more information. Error details: SAM
 file doesn't have any read groups defined in the header.  The GATK no longer supports SAM files without read groups
```

Lets check if the read groups were assigned:

### Checking if read groups were assigned and sorted
```
for f in *_sorted_rg.bam
do
  samtools view -H "${f}" | grep '^@RG'
done
```

There is an error with one file: 18-118.
```bash
samtools view: failed to open "18-118_S162.deduplicated_sorted_rg.bam" for reading: Exec format error
@RG     ID:1    LB:lib1 PL:illumina     SM:18-130_S172_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-142_S189_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-167_S166_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-178_S191_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-190_S186_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-202_S188_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-20_S202_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-227_S170_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-239_S185_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-250_S195_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-262_S179_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-311_S187_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-322_S180_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-32_S178_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-334_S164_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-346_S193_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-358_S201_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-370_S171_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-394_S192_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-406_S177_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-418_S196_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-442_S165_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-44_S198_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-454_S197_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-466_S199_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-55_S190_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-67_S176_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-79_S181_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-91_S160_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-9_S159_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam  PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-1029_S183_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-1038_S184_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-1053_S167_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-1059_S175_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-1093_S168_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-1257_S205_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-1263_S173_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-562_S174_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-571_S194_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-661_S182_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-704_S169_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-728_S161_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-862_S200_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-924_S204_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-933_S203_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
```


### Validating BAM files

Lets see if this BAM file is vaild

https://gatk.broadinstitute.org/hc/en-us/articles/360035891231-Errors-in-SAM-or-BAM-files-can-be-diagnosed-with-ValidateSamFile

```
java -jar $EBROOTPICARD/picard.jar ValidateSamFile \
        I=18-118_S162_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam \
        IGNORE_WARNINGS=true \
        MODE=VERBOSE
```

ERROR MESSAGE

```bash
(base) [kevin_wong1@n094 bismark_deduplicated]$ java -jar $EBROOTPICARD/picard.jar ValidateSamFile \
>         I=18-118_S162_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam \
>         IGNORE_WARNINGS=true \
>         MODE=VERBOSE
INFO    2023-10-24 10:56:57     ValidateSamFile

********** NOTE: Picard's command line syntax is changing.
**********
********** For more information, please see:
********** https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
**********
********** The command line looks like this in the new syntax:
**********
**********    ValidateSamFile -I 18-118_S162_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam -IGNORE_WARNINGS true -MODE VERBOSE
**********


10:56:58.093 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/software/picard/2.25.1-Java-11/picard.jar!/com/intel/gkl/native/libgkl_compression.so
[Tue Oct 24 10:56:58 EDT 2023] ValidateSamFile INPUT=18-118_S162_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam MODE=VERBOSE IGNORE_WARNINGS=true    MAX_OUTPUT=100 VALIDATE_INDEX=true INDEX_VALIDATION_STRINGENCY=EXHAUSTIVE IS_BISULFITE_SEQUENCED=false MAX_OPEN_TEMP_FILES=8000 SKIP_MATE_VALIDATION=false VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
[Tue Oct 24 10:56:58 EDT 2023] Executing as kevin_wong1@n094.cluster.com on Linux 3.10.0-1160.102.1.el7.x86_64 amd64; OpenJDK 64-Bit Server VM 11.0.2+9; Deflater: Intel; Inflater: Intel; Provider GCS is not available; Picard version: 2.25.1
WARNING 2023-10-24 10:56:58     ValidateSamFile NM validation cannot be performed without the reference. All other validations will still occur.
ERROR::MISSING_READ_GROUP:Read groups is empty
INFO    2023-10-24 10:56:58     SamFileValidator        Seen many non-increasing record positions. Printing Read-names as well.
 
INFO    2023-10-24 10:57:35     SamFileValidator        Validated Read    10,000,000 records.  Elapsed time: 00:00:36s.  Time for last 10,000,000:   36s.  Last read position: 000177F:273,416.  Last read name: A00547:195:HG2MYDSX2:4:2106:28438:13135_1:N:0:CAAGCTAG+CGCTATGT
INFO    2023-10-24 10:58:10     SamFileValidator        Validated Read    20,000,000 records.  Elapsed time: 00:01:12s.  Time for last 10,000,000:   35s.  Last read position: 000046F:343,128.  Last read name: A00547:195:HG2MYDSX2:4:1172:14660:10692_1:N:0:CAAGCTAG+CGCTATGT
INFO    2023-10-24 10:58:46     SamFileValidator        Validated Read    30,000,000 records.  Elapsed time: 00:01:47s.  Time for last 10,000,000:   35s.  Last read position: 000977F:43,335.  Last read name: A00547:195:HG2MYDSX2:4:2353:22417:4523_1:N:0:CAAGCTAG+CGCTATGT
[Tue Oct 24 10:59:00 EDT 2023] picard.sam.ValidateSamFile done. Elapsed time: 2.04 minutes.
Runtime.totalMemory()=2035417088
To get help, see http://broadinstitute.github.io/picard/index.html#GettingHelp
```

## Lets try indexing this file


`samtools index 18-118_S162.deduplicated_sorted_rg.bam`

This didn't work suggesting that the file is not sorted. Lets backtrack and repeat the previous sorting steps.


## Running this on the one problematic sample
```bash
samtools sort 18-118_S162_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam \
-o 18-118_S162.deduplicated_sorted.bam 
```

There was an error because temporary files were formed. I removed the tmp files and it ran... this was likely the source of our error. 

Running picard now: 

```bash
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=18-118_S162.deduplicated_sorted.bam O=18-118_S162.deduplicated_sorted_rg.bam LB=lib1 PL=illumina PU=unit1 SM=18-118_S162_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam
```

### Checking again if read groups were assigned and sorted
```
for f in *_sorted_rg.bam
do
  samtools view -H "${f}" | grep '^@RG'
done
```

YAYYY IT WORKED
```
@RG     ID:1    LB:lib1 PL:illumina     SM:18-106_S163_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-118_S162_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-130_S172_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-142_S189_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-167_S166_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-178_S191_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-190_S186_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-202_S188_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-20_S202_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-227_S170_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-239_S185_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-250_S195_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-262_S179_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-311_S187_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-322_S180_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-32_S178_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-334_S164_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-346_S193_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-358_S201_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-370_S171_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-394_S192_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-406_S177_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-418_S196_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-442_S165_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-44_S198_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-454_S197_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-466_S199_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-55_S190_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-67_S176_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-79_S181_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-91_S160_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:18-9_S159_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam  PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-1029_S183_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-1038_S184_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-1053_S167_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-1059_S175_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-1093_S168_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-1257_S205_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-1263_S173_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam        PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-562_S174_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-571_S194_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-661_S182_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-704_S169_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-728_S161_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-862_S200_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-924_S204_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
@RG     ID:1    LB:lib1 PL:illumina     SM:L-933_S203_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam PU:unit1
```


### indexing files

There were no .bai files, suggesting that the indexing did not run (probably because of the error in one file above). Re-running this now as a loop: 

```bash
for f in *_sorted_rg.bam
do
  samtools index "${f}"
done
```



# Re-run Bisulfite genotyper

`mkdir snp_vcfs`

`nano bssnp_geno.sh`

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

module load bis-snp/1.0.1.3-Java-8

for a in *rg.bam; 
do 
b=`echo ${a} | sed 's/.deduplicated_sorted_rg.bam/vfn/'` && bis-snp -Xmx10g -T BisulfiteGenotyper -C CG,1 -I ${a} -R past_filtered_assembly.fa -vfn1 ${b}1.vcf -vfn2 snp_vcfs/${b}2.vcf 
done
```