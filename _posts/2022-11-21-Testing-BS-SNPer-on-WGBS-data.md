---
layout: post
title: Testing BS-SNPer on WGBS data
date: '2022-11-21'
categories: Analysis
tags: WGBS, BS-SNPer
---

Code is inspired by [this post](https://github.com/RobertsLab/project-gigas-oa-meth/blob/master/code/07-BS-SNPer.ipynb)

In this notebook post, I'll use BS-Snper to call SNP variants from the Porites astreoides WGBS data from the Themal Transplant Molecular project. Adult ccolonies originated from 2 thermally disincct reef sites, subjected to ambient (28) or heated (31) temperature conditions, then transplanted to high thermal variability site for 11 months. After the transplantation, adults were sampled along with subsequent larve from each ccolony. There were 4 adult parental histories with a n=4 for each history, with a larval pair for each colony sampled. DNA was extracted from whole tissue for WGBS. Reads were aligned to the *P. astreoides* genome through the nf-ccore methylseq pipeline.


# 1. Set working directory

`cd /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS`

`mkdir BS-SNPer`

# 2. Identify SNP variants

I will identify variants in individual files, as well as SNPs across all samples.

## 2a. Merge SNP variants

First I need to sort all the deduplicated bam files

`cd /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated`

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

for f in *.deduplicated.bam
do
  STEM=$(basename "${f}" _L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam)
  samtools sort "${f}" \
  -o "${STEM}".deduplicated_sorted.bam
done
```
*revisit this to see if this sort actually worked*

To identify SNPs across all samples, I need to merge my samples, then use that as the input file for BS-Snper.

`nano merge_bam.sh`

```bash
#!/bin/bash
#SBATCH --job-name="merge_snp"
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

# Merge Samples with SAMtools

samtools merge \
merged-sorted-deduplicated.bam \
*sorted.bam

```

View output file header


`samtools view merged-sorted-deduplicated.bam | head `

```
A00547:195:HG2MYDSX2:4:2253:28881:9518_1:N:0:CTTGGTAT+GGACTTGG  99      000000F 2       3       106M    =       30      129     GGAGATATAAAATTTTTTTTTGAGTGTTGAAAAATATTTTGCGAGTGAGTGTAGTGAACGAGTGAAATATTTTTTTTAATACGAGAAGAGAAATTTCGTATTTTTA   FFF:FF:FFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF:FFFFF:F:FFFFFFFFFFF      NM:i:18 MD:Z:2T4C0G6C1C2C4C13C0A0T7C1C2C21C2C21C1C0C1   XM:Z:.......z.......h.h..z..................h.........z.x..z...Z.................h..h.Z..............Z....h.hh.      XR:Z:CT XG:Z:CT
A00547:195:HG2MYDSX2:4:2230:2437:27790_1:N:0:GCAATGCA+AACGTTCC  99      000000F 9       2       106M    =       104     201     TGAAATTTTTTTTTGAGTGTTGAAAAATATTTTATGAGTGAGTGTAGTGAATGAGTGAAATATTTTTTTTAATATGAGAAGAGATATTTTGTATTTTTAAGAGATT   FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF      NM:i:21 MD:Z:0C7C1C2C4C13C9C1C2C3C17C2C1C9A4C4C1C0C3C2C0C0      XM:Z:z.......h.h..z..................h.........z.x..z...z.................h..h.z..............z....h.hh......hh      XR:Z:CT XG:Z:CT
A00547:195:HG2MYDSX2:4:1407:26793:24283_1:N:0:CTACGACA+TTGGACTC 163     000000F 16      0       105M    =       165     255     TCTCTTCAAATATTAAAAAATATTTCACAAATAAACACAACAAACAAATAAAATATTTTTTTCAACACAAAAAAAAACATTTCATATCTCCAAAAAACCATATAA    FFFFFFFFFFFFFFFFFFFFFFFFFFFF:,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF       NM:i:25 MD:Z:7G1G1C2G12T0G1G1G1G1G2G1G3G1G1G18G1G2G1G1A5G9G0C0G5G3      XM:Z:.......z.h....h.............h.h.h.h.z..x.z...z.h.h..................z.h..h.h.......z.........h.z.....h...       XR:Z:GA XG:Z:GA
A00547:195:HG2MYDSX2:4:2253:28881:9518_1:N:0:CTTGGTAT+GGACTTGG  147     000000F 30      3       85M1D2M1D8M7I4M =       2       -129    GAAAAATATTTTGCGAGTGAGTGTAGTGAACGAGTGAAATATTTTTTTTAATACGAGAAGAGAAATTTCGTATTTTTAAGCGATTTGATCGTTTTTGTTGGGATGT   FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF      NM:i:26 MD:Z:11C0A0T7C1C2C21C2C21C1C0C6C0C0^A2^T1A0T3C3T1       XM:Z:...........h.........z.x..z...Z.................h..h.Z..............Z....h.hh...Z..hx........u............      XR:Z:GA XG:Z:CT
A00547:195:HG2MYDSX2:4:1105:7735:26318_1:N:0:GGTACCTT+GACGTCTT  163     000000F 31      8       105M    =       82      150     AAAAATATTTCACAAATAAACACAACAAACAAATAAAATATTTTTTTCAACACAAAAAAAAACATTTCATATCTCCAAAAAACCATATAATATTCTATTTATTAT    FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF       NM:i:22 MD:Z:12T0G1G1G1G1G2G1G3G1G1G18G1G2G1G1A5G9G0C0G5G4G13   XM:Z:.............h.h.h.h.z..x.z...z.h.h..................z.h..h.h.......z.........h.z.....h....h.............       XR:Z:GA XG:Z:GA
A00547:195:HG2MYDSX2:4:1522:10303:17206_1:N:0:CGGCGTGA+ACAGGCGC 163     000000F 36      2       5M1I1M3I4M1I91M =       41      109     CATCTACAAAATAATATAAATACAACGAACGAATAAAATATTTTTTTCGACACGAAAAAAAAAATTTCGCATCTCCAAACAACCATATAATATTCTATTTAATATA   F:FFFFFFFF:,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF:FFF:FFFFFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFFF:FFFFFF:      NM:i:26 MD:Z:0T2T4G1G1G1G0C0G2G7G1G13A6G2G1G8T8G1G5G4G9T4       XM:Z:............h..u.h.h.z..x.Z...Z.h.h..................Z.h..h.h.......Z.........h.z.....h....h..............      XR:Z:GA XG:Z:GA
A00547:195:HG2MYDSX2:4:1522:10303:17206_1:N:0:CGGCGTGA+ACAGGCGC 83      000000F 41      2       5M1I99M =       36      -109    AATAATATAAATACACCGAACGAATAAAATATTTTTTTCGACACGAAAAAAAAAATTTCGCATCTCCAAACCACCCTATAATATTCTATTTAATATATAAACCCC    FFFFFFF,FFFFFF:,F:F:FFFF:F:FF,F,,FFFFFFFFFF:FFFFFFFFFFFFFFFFF,FFFFFF::F,FFF,F,,F:F,,FF::F,F:FFFF:FFF,F,FF       NM:i:23 MD:Z:0C2G1G1G1G0C0G2G7G1G13A6G2G1G8T8G1G3A1G4G9T9A2     XM:Z:...h..u.h.h.z....Z...Z.h.h..................Z.h..h.h.......Z.........h.......h....h......................       XR:Z:CT XG:Z:GA
A00547:195:HG2MYDSX2:4:2678:12292:36714_1:N:0:TAAGTGGT+GGCTTAAG 99      000000F 44      2       28M2I68M2I6M    =       55      112     AAGTGAGGGTAGTAAATGAGTGAAATATTATTTTTTTAATATAAGAAGAGAAATTTTGTATTTTTAAGTGATTATGTAATGTTTTATTTATTTTATAAATATAGTA   FFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFF:FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFF      NM:i:26 MD:Z:0G6C1C2C0G2C17C2C1C0G13C4C1C0C3C2C0C10C8A6C1C0C1   XM:Z:.........x..z...z...................h..h.z..............z....h.hh...z..hh..........h.................h..h.      XR:Z:CT XG:Z:CT
A00547:195:HG2MYDSX2:4:1273:17381:7639_1:N:0:AAGTCCAA+TACTCATA  163     000000F 44      0       101M    =       56      118     AAATAAACACAACAAACAAATAAAATATTTTTTTCAAAACAAAAAAAAAAATTCCATATCTCCAAACAACCATATAATATTCTATTTATTATATAAACACC        FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF   NM:i:21 MD:Z:0G1G1G1G1G2G1G3G1G1G15C2G1G2G1G5T1G9G1G5G4G22      XM:Z:h.h.h.h.z..x.z...z.h.h..................z.h..h.h.......z.........h.z.....h....h......................   XR:Z:GA XG:Z:GA
A00547:195:HG2MYDSX2:4:1570:9932:29027_1:N:0:GACCTGAA+CTCACCAA  163     000000F 49      6       90M     =       49      90      AACACAACTAACAAATAAAATATTTTTTTCAACACAAAAAAAAAAATTTCATATCTCCAAACAACCATATAATATTCTATTTATAATATA  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF       NM:i:17 MD:Z:1G1G2G1G3G1G1G18G1G2G1G7G9G1G5G4G10T5      XM:Z:.h.z..x.....z.h.h..................z.h..h.h.......z.........h.z.....h....h................      XR:Z:GA XG:Z:GA
```

Create index file for IGV

`nano index_bam.sh`

```bash
#!/bin/bash
#SBATCH --job-name="index_bam"
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

# Index merged bam file 

samtools index -b merged-sorted-deduplicated.bam 

find merged-sorted-deduplicated.bam.bai
```


## 2b. Identify SNPs 

Options for the script are found [here](https://github.com/hellbelly/BS-Snper/blob/master/README.txt) and below. 

```
--fa: Reference genome file in fasta format
--input: Input bam file (I'm using deduplicated sorted bams)
--output: Temporary file storing SNP candidates
--methcg: CpG methylation information
--methchg: CHG methylation information
--methchh: CHH methylation information
--minhetfreq: Threshold of frequency for calling heterozygous SNP
--minhomfreq: Threshold of frequency for calling homozygous SNP
--minquali: Threshold of base quality
--mincover: Threshold of minimum depth of covered reads
--maxcover: Threshold of maximum depth of covered reads
--minread2: Minimum mutation reads number
--errorate: Minimum mutation rate
--mapvalue: Minimum read mapping value
SNP.out: Final SNP result file
ERR.log: Log file
```

You can run BS-SNPer in Linux or MAC OS, using the command like:

```
perl BS-Snper.pl <sorted_bam_file> --fa <reference_file> --output <snp_result_file> --methcg <meth_cg_result_file> --methchg <meth_chg_result_file> --methchh <meth_chh_result_file> --minhetfreq 0.1 --minhomfreq 0.85 --minquali 15 --mincover 10 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20 >SNP.out 2>ERR.log
```

**Attention**: Both of the input and output file arguments should be passed to BS-SNPer in the form of absolute paths. 

`nano bs_snper_merged.sh`

```bash
#!/bin/bash
#SBATCH --job-name="BS_snper"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/BS-SNPer/merged_SNP_output

module load BS-Snper/1.0-foss-2021b

perl /opt/software/BS-Snper/1.0-foss-2021b/bin/BS-Snper.pl \
--fa /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--input /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated/merged-sorted-deduplicated.bam \
--output /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/BS-SNPer/merged_SNP_output/SNP-candidates.txt \
--methcg /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/BS-SNPer/merged_SNP_output/CpG-meth-info.tab \
--methchg /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/BS-SNPer/merged_SNP_output/CHG-meth-info.tab \
--methchh /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/BS-SNPer/merged_SNP_output/CHH-meth-info.tab \
--mincover 5 \
> /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/BS-SNPer/merged_SNP_output/SNP-results.vcf 2> /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/BS-SNPer/merged_SNP_output/merged.ERR.log
```

Current issue: 

```
Unknown option: input
FLAG: 1
refSeqFile = /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta.
bamFileName = /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated/merged-sorted-deduplicated.bam.
snpFileName = /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/BS-SNPer/merged_SNP_output/SNP-candidates.txt.
methCgFileName = /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/BS-SNPer/merged_SNP_output/CpG-meth-info.tab.
methChgFileName = /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/BS-SNPer/merged_SNP_output/CHG-meth-info.tab.
methChhFileName = /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/BS-SNPer/merged_SNP_output/CHH-meth-info.tab.
vQualMin = 15.
nLayerMax = 1000.
vSnpRate = 0.100000.
vSnpPerBase = 0.020000.
mapqThr = 20.
Too many characters in one row! Try to split the long row into several short rows (fewer than 1000000 characters per row).
Error! at /opt/software/BS-Snper/1.0-foss-2021b/bin/BS-Snper.pl line 110.
```

## 2c. Individual SNP variants

## 2d. Unique SNP variants

`mkdir all_SNP_output`

`nano bs_snper_all.sh`

```bash
#!/bin/bash
#SBATCH --job-name="BS_snper"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/BS-SNPer/all_SNP_output

# load modules
module load BS-Snper/1.0-foss-2021b

# symbolically link files
ln -s /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated/*.deduplicated_sorted.bam .

# Loop BSsnper for each file
FILES=$(ls *.deduplicated_sorted.bam)
echo ${FILES}

for file in ${FILES}
do
    NAME=$(echo ${file} | awk -F "." '{print $1}')
    echo ${NAME}

    perl /opt/software/BS-Snper/1.0-foss-2021b/bin/BS-Snper.pl \
    --fa /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
    --input ${NAME}.deduplicated_sorted.bam \
    --output ${NAME}.SNP-candidates.txt \
    --methcg ${NAME}.CpG-meth-info.tab \
    --methchg ${NAME}.CHG-meth-info.tab \
    --methchh ${NAME}.CHH-meth-info.tab \
    --mincover 5 \
    > ${NAME}.SNP-results.vcf 2> ${NAME}.ERR.log

done
```

I keep on getting the same error:

```
Unknown option: input
FLAG: 1
refSeqFile = /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta.
bamFileName = L-933_S203.deduplicated_sorted.bam.
snpFileName = L-933_S203.SNP-candidates.txt.
methCgFileName = L-933_S203.CpG-meth-info.tab.
methChgFileName = L-933_S203.CHG-meth-info.tab.
methChhFileName = L-933_S203.CHH-meth-info.tab.
vQualMin = 15.
nLayerMax = 1000.
vSnpRate = 0.100000.
vSnpPerBase = 0.020000.
mapqThr = 20.
Too many characters in one row! Try to split the long row into several short rows (fewer than 1000000 characters per row).
Error! at /opt/software/BS-Snper/1.0-foss-2021b/bin/BS-Snper.pl line 110.
```


I am going to try on the original bam files (potentially not sorted)

`mkdir all_notsorted_SNP_output`

`nano bs_snper_all_notsorted.sh`

```bash
#!/bin/bash
#SBATCH --job-name="BS_snper"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/BS-SNPer/all_notsorted_SNP_output

# load modules
module load BS-Snper/1.0-foss-2021b

# symbolically link files
ln -s /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated/*_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam .

# Loop BSsnper for each file
FILES=$(ls *_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam)
echo ${FILES}

for file in ${FILES}
do
    NAME=$(echo ${file} | awk -F "." '{print $1}')
    echo ${NAME}

    perl /opt/software/BS-Snper/1.0-foss-2021b/bin/BS-Snper.pl \
    --fa /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
    --input ${NAME}.deduplicated.sorted.bam \
    --output ${NAME}.SNP-candidates.txt \
    --methcg ${NAME}.CpG-meth-info.tab \
    --methchg ${NAME}.CHG-meth-info.tab \
    --methchh ${NAME}.CHH-meth-info.tab \
    --mincover 5 \
    > ${NAME}.SNP-results.vcf 2> ${NAME}.ERR.log

done
```

Same error: 

```
Unknown option: input
FLAG: 1
refSeqFile = /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta.
bamFileName = L-933_S203_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.sorted.bam.
snpFileName = L-933_S203_L004_R1_001_val_1_bismark_bt2_pe.SNP-candidates.txt.
methCgFileName = L-933_S203_L004_R1_001_val_1_bismark_bt2_pe.CpG-meth-info.tab.
methChgFileName = L-933_S203_L004_R1_001_val_1_bismark_bt2_pe.CHG-meth-info.tab.
methChhFileName = L-933_S203_L004_R1_001_val_1_bismark_bt2_pe.CHH-meth-info.tab.
vQualMin = 15.
nLayerMax = 1000.
vSnpRate = 0.100000.
vSnpPerBase = 0.020000.
mapqThr = 20.
Too many characters in one row! Try to split the long row into several short rows (fewer than 1000000 characters per row).
Error! at /opt/software/BS-Snper/1.0-foss-2021b/bin/BS-Snper.pl line 110.
```


## 20240220 Attempt

we are following Danielle's post now: https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2024-02-12-Testing-BS-SNPer-Molec-Underpinnings-WGBS.md


`nano BS_SNPER.sh`

```bash
#!/bin/bash
#SBATCH --job-name="BS_snper"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/BS-SNPer/merged_SNP_output/

# load modules
module load BS-Snper/1.0-foss-2021b

perl /opt/software/BS-Snper/1.0-foss-2021b/bin/BS-Snper.pl \
--fa /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--input /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated/merged-sorted-deduplicated.bam \
--output SNP-candidates.txt \
--methcg CpG-meth-info.tab \
--methchg CHG-meth-info.tab \
--methchh CHH-meth-info.tab \
--minhetfreq 0.1 \
--minhomfreq 0.85 \
--minquali 15 \
--mincover 5 \
--maxcover 1000 \
--minread2 2 \
--errorate 0.02 \
--mapvalue 20 \

>SNP-results.vcf 2>SNP.log

```

According to Kevin Bryan and Danielle, the number of scaffold/chromosomes and lines must be adjusted in the program. The following code determines how many lines are in our genome file: 

`cd /data/putnamlab/kevin_wong1/Past_Genome`

`nano max_ln.py`

````
import sys
mx = 0
with open(sys.argv[1]) as fd:
 for ln in fd.readlines():
  mx = max(len(ln), mx)
print(mx)
````

`interactive`

`module load Python/3.10.8-GCCcore-12.2.0`

`python max_ln.py past_filtered_assembly.fasta` 

This outputs: 3369716

I asked Kevin Bryan to increase the line input to accomodate 3.3 million lines and re-ran the same script. It is running now!!






