---
layout: post
title: Old Tag-seq analysis for PJB
date: '2022-04-27'
categories: Analysis
tags: TagSeq
---

# Porites bleaching TagSeq Bioinformatics

## Script modified from my awesome lab mates [S. Gurr](https://github.com/SamGurr/Pgenerosa_TagSeq_Metabolomics/blob/main/TagSeq/HPC_work/Geoduck_TagSeq_Bioinf.md#Initial-diagnostics-upon-sequence-upload-to-HPC) and [E. Chille]().


# Initial diagnostics upon sequence upload to HPC
--------------------------------------------

HP uploaded data to Bluewaves and we are waiting for the md5 checksum from UT.

HP ran initial multiQC

# Trimming and post-trim quality check of 'clean' reads
-------------------------------------------

### shell script: <span style="color:green">**fastp_multiqc.sh**<span>

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --output=../../../../kevin_wong1/20210315_Past_Tagseq/output/clean/"%x_out.%j"
#SBATCH --error=../../../../kevin_wong1/20210315_Past_Tagseq/output/clean/"%x_err.%j"
#SBATCH -D /data/putnamlab/KITT/hputnam/20210312_Pastreoides_TagSeq/Pastreoides_TagSeq_JA21015
#SBATCH --cpus-per-task=3

# load modules needed
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.7-foss-2018b-Python-2.7.15

# Make an array of sequences to trim
array1=($(ls *.fastq.gz))

# fastp loop; trim the Read 1 TruSeq adapter sequence; trim poly x default 10 (to trim polyA)
for i in ${array1[@]}; do
	fastp --in1 ${i} --out1 ../../../../kevin_wong1/20210315_Past_Tagseq/output/clean/clean.${i} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 20 -y -Y 50
        fastqc ../../../../kevin_wong1/20210315_Past_Tagseq/output/clean/clean.${i}
done

echo "Read trimming of adapters complete." $(date)

# Quality Assessment of Trimmed Reads

cd ../../../../kevin_wong1/20210315_Past_Tagseq/output/clean/ #The following command will be run in the /clean directory

multiqc ./ #Compile MultiQC report from FastQC files

echo "Cleaned MultiQC report generated." $(date)
```

20210318 job #1883852

### EXPORT MULTIQC REPORT
*exit bluewaves and run from terminal*
- save to gitrepo as multiqc_clean.html
```
scp kevin_wong1@bluewaves.uri.edu:/data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/clean/multiqc_report.html  MyProjects/Porites_Rim_Bleaching_2019/output/TagSeq
```


# HPC Job: HISAT2 Index Reference and Alignment
-----------------------------------------------------------------

This step was not totally necessary, as I am not using these outputs. However, this step did help me determine the mapping percentages using the Kenkel et al. 2013 reference transcriptome.

- create directory output\hisat2

``` mkdir HISAT2 ```


### Mapping #1 - Porites astreodes transcriptome from [Kenkel et al. 2013](https://matzlab.weebly.com/data--code.html)

- index reference and alignment

uploading reference transcriptome from [Kenkel et al. 2013](https://matzlab.weebly.com/data--code.html)

```
scp Desktop/URI_PHD/pastreoides_2014_ref/pastreoides_may2014/past.fasta kevin_wong1@bluewaves.uri.edu:/data/putnamlab/kevin_wong1/REFS/Past/Kenkel2013_past_transcriptome.fasta
```

**input**
- Kenkel2013_past_transcriptome.fasta *= reference transcriptome*
- clean/*.fastq.gz *= all clean TagSeq reads*

**ouput**
- Past_transcriptome_ref *= indexed reference by hisat2-build; stored in the output/hisat2 folder as 1.hy2, 2.ht2... 8.ht2*
- <clean.fasta>.sam *=hisat2 output, readable text file; removed at the end of the script*
- <clean.fasta>.bam *=converted binary file complementary to the hisat sam files*

### shell script: <span style="color:green">**HISAT2.sh**<span>

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
#SBATCH -D /data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/hisat2/Kenkel2013_REF
#SBATCH --cpus-per-task=3

#load packages
module load HISAT2/2.1.0-foss-2018b #Alignment to reference genome: HISAT2
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

# symbolically link 'clean' reads to hisat2 dir
ln -s ../../../clean/clean*.fastq.gz ./

# index the reference transcriptome for Porites astreoides output index to working directory
hisat2-build -f ../../../../REFS/Past/Kenkel2013_past_transcriptome.fasta ./Past_transcriptome_ref # called the reference transcriptome (scaffolds)
echo "Referece transcriptome indexed. Starting alingment" $(date)

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed
array=($(ls *.fastq.gz)) # call the symbolically linked sequences - make an array to align
for i in ${array[@]}; do
        sample_name=`echo $i| awk -F [.] '{print $2}'`
        hisat2 -p 8 --dta -x Past_transcriptome_ref -U ${i} -S ${sample_name}.sam
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
                echo "${i} bam-ified!"
        rm ${sample_name}.sam
done

```

20210319 Submitted job 1883908

Took ~2 days

This job "failed" but I think it was because I had an extra "clean*.fastq.qz" file. All of the samples had mapping percentages and sam files were converted too bam then deleted.

```
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
Extra parameter(s) specified: "clean.A10-KW10_S233_L002_R1_001.fastq.gz", "clean.A11-KW11_S236_L001_R1_001.fastq.gz", "clean.A11-KW11_S236_L002_R1_001.fastq.gz", "clean.A12-KW12_S239_L001_R1_001.fastq.gz", "clean.A12-KW12_S239_L002_R1_001.fastq.gz", "clean.A$
Error: Encountered internal HISAT2 exception (#1)
Command: /opt/software/HISAT2/2.1.0-foss-2018b/bin/hisat2-align-s --wrapper basic-0 -p 8 --dta -x Past_transcriptome_ref -S A10-KW10_S233_L001_R1_001.sam -U /tmp/31840.unp clean.A10-KW10_S233_L002_R1_001.fastq.gz clean.A11-KW11_S236_L001_R1_001.fastq.gz clea$
(ERR): hisat2-align exited with value 1
[E::hts_open_format] Failed to open file A10-KW10_S233_L001_R1_001.sam
samtools sort: can't open "A10-KW10_S233_L001_R1_001.sam": No such file or directory
rm: cannot remove `A10-KW10_S233_L001_R1_001.sam': No such file or directory
```

Exporting to report mapping percentages

```
scp kevin_wong1@bluewaves.uri.edu:/data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/hisat2/Kenkel2013_REF/HISAT2.sh_err.1883908  MyProjects/Porites_Rim_Bleaching_2019/output/TagSeq
```


# Salmon
-----------------------------------------------------------------

Replaces HISAT2 and StringTie because I have no reference genome with a GFF3 file

- Need to also try with SAM/BAM files

### Merging lanes of cleaned reads for Salmon

```
#!/bin/bash
#SBATCH --job-name="cat.clean"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/clean
#SBATCH --mem=100GB

module load SAMtools/1.9-foss-2018b

echo "Concatenate files" $(date)

ls *R1_001.fastq.gz | awk -F '[_]' '{print $1"_"$2}' | sort | uniq > ID

for i in `cat ./ID`;
	do cat $i\_L001_R1_001.fastq.gz $i\_L002_R1_001.fastq.gz > $i\_ALL.fastq.gz;
	done

echo "Mission complete." $(date)
```

### Salmon shell script

https://combine-lab.github.io/salmon/getting_started/

```
#!/bin/bash
#SBATCH --job-name="Salmon"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/salmon
#SBATCH --mem=100GB

module load Salmon/1.1.0-gompi-2019b

#echo "Building index..." $(date)
#salmon index -t ../../../REFS/Past/Kenkel2013_past_transcriptome.fasta -i PAST_Kenkel_ref_index -k 31

# symbolically link 'clean' reads to salmon dir
ln -s ../clean/*_ALL.fastq.gz ./

echo "Running the alignment and abundance estimation." $(date)

array1=($(ls *_ALL.fastq.gz))
for i in ${array1[@]}; do
salmon quant -i PAST_Kenkel_ref_index -l A \
        -r ${i} \
        --validateMappings -p 8 -o quants/${i}_quant
done

echo "Mission complete." $(date)

```

### Quant merge

This takes all the quant.sf files in each individual folder and combined them into a single gene count matrix.

```
cd /data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/salmon/quants
```

```
salmon quantmerge --quants {clean.A10-KW10_S233_ALL.fastq.gz_quant,clean.A11-KW11_S236_ALL.fastq.gz_quant,clean.A12-KW12_S239_ALL.fastq.gz_quant,clean.A1-KW1_S197_ALL.fastq.gz_quant,clean.A2-KW2_S201_ALL.fastq.gz_quant,clean.A3-KW3_S205_ALL.fastq.gz_quant,clean.A4-KW4_S209_ALL.fastq.gz_quant,clean.A5-KW5_S213_ALL.fastq.gz_quant,clean.A6-KW6_S217_ALL.fastq.gz_quant,clean.A7-KW7_S221_ALL.fastq.gz_quant,clean.A8-KW8_S225_ALL.fastq.gz_quant,clean.A9-KW9_S229_ALL.fastq.gz_quant,clean.B10-KW22_S234_ALL.fastq.gz_quant,clean.B11-KW23_S237_ALL.fastq.gz_quant,clean.B12-KW24_S240_ALL.fastq.gz_quant,clean.B1-KW13_S198_ALL.fastq.gz_quant,clean.B2-KW14_S202_ALL.fastq.gz_quant,clean.B3-KW15_S206_ALL.fastq.gz_quant,clean.B4-KW16_S210_ALL.fastq.gz_quant,clean.B5-KW17_S214_ALL.fastq.gz_quant,clean.B6-KW18_S218_ALL.fastq.gz_quant,clean.B7-KW19_S222_ALL.fastq.gz_quant,clean.B8-KW20_S226_ALL.fastq.gz_quant,clean.B9-KW21_S230_ALL.fastq.gz_quant,clean.C10-KW34_S235_ALL.fastq.gz_quant,clean.C11-KW35_S238_ALL.fastq.gz_quant,clean.C12-KW36_S241_ALL.fastq.gz_quant,clean.C1-KW25_S199_ALL.fastq.gz_quant,clean.C2-KW26_S203_ALL.fastq.gz_quant,clean.C3-KW27_S207_ALL.fastq.gz_quant,clean.C4-KW28_S211_ALL.fastq.gz_quant,clean.C5-KW29_S215_ALL.fastq.gz_quant,clean.C6-KW30_S219_ALL.fastq.gz_quant,clean.C7-KW31_S223_ALL.fastq.gz_quant,clean.C8-KW32_S227_ALL.fastq.gz_quant,clean.C9-KW33_S231_ALL.fastq.gz_quant,clean.D1-KW37_S200_ALL.fastq.gz_quant,clean.D2-KW38_S204_ALL.fastq.gz_quant,clean.D3-KW39_S208_ALL.fastq.gz_quant,clean.D4-KW40_S212_ALL.fastq.gz_quant,clean.D5-KW41_S216_ALL.fastq.gz_quant,clean.D6-KW42_S220_ALL.fastq.gz_quant,clean.D7-KW43_S224_ALL.fastq.gz_quant,clean.D8-KW44_S228_ALL.fastq.gz_quant,clean.D9-KW45_S232_ALL.fastq.gz_quant} --column numreads --output /data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/salmon/full_count_matrix.txt

```

```
scp kevin_wong1@bluewaves.uri.edu:/data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/salmon/full_count_matrix.txt  MyProjects/Porites_Rim_Bleaching_2019/output/TagSeq
```
