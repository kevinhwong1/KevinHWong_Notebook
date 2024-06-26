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

# Sort files

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

This script completed in 9 days. 


# Tabulate vcfs

First we need to get scaffold lengths of our genome:

```bash
interactive

cd /data/putnamlab/kevin_wong1/Past_Genome

module load pyfaidx/0.7.0-GCCcore-11.3.0

faidx past_filtered_assembly.fasta -i chromsizes > past_scaffold_lengths.tsv
```

Second, we need to modify the script to our own paths:

`wget https://github.com/lyijin/pdae_dna_meth/blob/master/genetic_contribution/bissnp/tabulate_snp_vcfs.py`


Copy this file to the same directory because it was struggling trying to find this file in other directories. 

`cp ~/../../data/putnamlab/kevin_wong1/Past_Genome/past_scaffold_lengths.tsv .`


Install [natsort](https://pypi.org/project/natsort/)

`pip install --kevin_wong1 natsort`


Modify the python file 

`nano tabulate_snp_vcfs.py`

```bash
#!/usr/bin/env python3

"""
> tabulate_snp_vcfs.py <

Compiles SNPs called by Bis-SNP from multiple files, then dumps genotype
information from all files into a giant tsv table.
"""
import argparse
import csv
import re

import numpy as np

from natsort import natsort

def convert_gt_to_int(gt_string):
    """
    This function converts GT calls from GATK into ints:
        0/0 (homozygous reference) --> 0
        0/1 (heterozygous reference) --> 1
        1/1, 1/2, ... (homozygous non-reference) --> 2
    
    Most homozygous non-reference bases are 1/1 (i.e. homozygous for first
    gt allele) anyway.
    """
    gt_ints = [int(x) for x in gt_string.split('/')]
    sum_gt_ints = sum(gt_ints)
    
    if sum_gt_ints > 2:
        sum_gt_ints = 2
    
    return sum_gt_ints

parser = argparse.ArgumentParser(description="""
Compiles SNPs called by Bis-SNP from multiple files, then dumps genotype
information from all files into a giant tsv table.""")

parser.add_argument('vcfs', metavar='vcf_files',
                    type=argparse.FileType('r'), nargs='+',
                    help='VCFs containing coverage information.')

args = parser.parse_args()

vcf_filenames = natsort.natsorted([x.name for x in args.vcfs])

# read styl scaffold lengths to create appropriately-sized NumPy arrays
scaf_lens = {}

tsv_reader = csv.reader(open(
    'past_scaffold_lengths.tsv'), delimiter='\t')
for row in tsv_reader:
    if not row: continue
    
    scaf_lens[row[0]] = int(row[1])

# chuck gt info into a dict containing many NumPy arrays
#   gt_info[file][scaf] = convert_gt_to_int(gt_bases)
gt_info = {}

for v in vcf_filenames:
    gt_info[v] = {}
    for s in scaf_lens:
        gt_info[v][s] = np.zeros(scaf_lens[s], dtype='int8')
    
    tsv_reader = csv.reader(open(v), delimiter='\t')
    for row in tsv_reader:
        if not row: continue
        if row[0][0] == '#': continue       # ignore commented lines
        
        scaf = row[0]
        pos = int(row[1]) - 1               # vcf files are 1-based
        gt_string = row[9].split(':')[0]
        gt_info[v][scaf][pos] = convert_gt_to_int(gt_string)

# print stuff out
print ('scaf', 'pos', *natsort.natsorted(gt_info), sep='\t')
for s in natsort.natsorted(scaf_lens):
    for n in range(scaf_lens[s]):
        pos_covs = [gt_info[x][s][n] for x in natsort.natsorted(gt_info)]
        
        # do not print positions that are 0 coverage in all files
        if not any(pos_covs): continue
        
        # remember to convert 0-based positions back to 1-based positions
        results = [s, n + 1] + pos_covs
        print (*results, sep='\t')
```


`nano tabulate_vcfs.sh`

```bash
#!/bin/bash
#SBATCH --job-name="tab_vcfs"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated/snp_vcfs

# load modules needed

module load Python/3.10.8-GCCcore-12.2.0
module load SciPy-bundle/2023.02-gfbf-2022b

# tabulate
python3 tabulate_snp_vcfs.py *.vcf > tabulated_genotypes.tsv

```

## Calculate fst

* [vcf files](http://samtools.github.io/hts-specs/VCFv4.2.pdf)
* Left off - 

https://github.com/lyijin/pdae_dna_meth/blob/master/genetic_contribution/calc_indiv_fst/calc_hudson_fst.py
- download tsv file and see what it looks like
- figure out filtering prior to FST analysis
- choose which analysis we want to perform and how to blacklist snps/genes 



`mkdir fst`

`nano fst.sh`

```bash
#!/bin/bash
#SBATCH --job-name="fst"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated/snp_vcfs


#load modules
module load VCFtools/0.1.16-GCC-11.2.0

#calculate fast
vcftools --gzvcf *.vcf \
--weir-fst-pop \
--out ./fst

```