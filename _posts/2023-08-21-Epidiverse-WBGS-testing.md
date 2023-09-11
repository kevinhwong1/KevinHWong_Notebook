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