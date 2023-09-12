---
layout: post
title: Epidiverse WBGS testing
date: '2023-08-21'
categories: Analysis
tags: WGBS, Porites astreoides
---

Testing the WGBS methylation calling pipeline from [Epidiverse](https://github.com/EpiDiverse/wgbs/blob/master/docs/usage.md#running-the-pipeline). 

`mkdir /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/wgbs`

`mkdir /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/wgbs/raw`

`mkdir /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/wgbs/input`

`mkdir /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/wgbs/output`

#symbiolically link genome files
ln s- /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta ./input
ln s- /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta.fai ./input

```bash
#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=18
#SBATCH --export=NONE
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/wgbs/
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed (specific need for my computer)
#source /usr/share/Modules/init/sh # load the module function

# load modules needed
echo "START" $(date)
#module load Anaconda3/2022.05
module load Nextflow/20.07.1 #this pipeline requires this version 

#run epidiverse wgbs pipeline
nextflow run epidiverse/wgbs \
-profile docker \
--input /data/putnamlab/KITT/hputnam/20211008_Past_ThermalTransplant_WGBS/*{1,2}.fastq.gz\
--reference input/past_filtered_assembly.fasta \
--INDEX \
--output output \
--trim \
--clip5 10 \
--clip3 10 \
--fastqc \
--noLambda 

```

# 20230911

So we think we have to copy the files into the raw folder and run it from there. I will make a script to copy the files the re-run wgbs. 

`raw.cp.sh`

```bash
#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=18
#SBATCH --export=NONE
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/wgbs/raw
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file


# load modules needed
echo "START" $(date)

cp /data/putnamlab/KITT/hputnam/20211008_Past_ThermalTransplant_WGBS/*.fastq.gz ./
```

move the error files out of the folder so it is only fastqs

`mv r* ../`

I have to rename the files to remove the _001 since the pipeline requires a *{1,2}.fastq.gz format

```bash
interactive
for file in *fastq.gz; do
    mv "$file" "${file/_001/}"
done
```

Run the wgbs pipeline

`nano epidiverse_wgbs.sh`

```bash
#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=18
#SBATCH --export=NONE
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/wgbs/
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed (specific need for my computer)
#source /usr/share/Modules/init/sh # load the module function

# load modules needed
echo "START" $(date)
#module load Anaconda3/2022.05
module load Nextflow/20.07.1 #this pipeline requires this version 

#run epidiverse wgbs pipeline
nextflow run epidiverse/wgbs \
-profile docker \
--input raw \
--reference input/past_filtered_assembly.fasta \
--INDEX \
--output output \
--trim \
--clip5 10 \
--clip3 10 \
--fastqc \
--noLambda 
```

This runs without terminating the job but I think I still have a Docker issue...

```bash
Error executing process > 'WGBS:read_trimming (18-130_S172_L004_R)'

Caused by:
  Process `WGBS:read_trimming (18-130_S172_L004_R)` terminated with an error exit status (127)

Command executed:

  mkdir fastq fastq/logs
  cutadapt -j 2 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -u 10 -u -10 \
  -q 20 -m 36 -O 3 \
  -o fastq/18-130_S172_L004_R1.fastq.gz \
  -p fastq/18-130_S172_L004_R2.fastq.gz 18-130_S172_L004_R1.fastq.gz 18-130_S172_L004_R2.fastq.gz \
  > fastq/logs/cutadapt.18-130_S172_L004_R.input.log 2>&1

Command exit status:
  127

Command output:
  (empty)

Command error:
  .command.run: line 259: docker: command not found

Work dir:
  /glfs/brick01/gv0/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/wgbs/work/ce/8eae537f967936c0ed58a31d005c09

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```