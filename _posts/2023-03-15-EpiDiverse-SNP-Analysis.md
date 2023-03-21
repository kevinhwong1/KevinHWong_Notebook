---
layout: post
title: EpiDiverse SNP Analysis
date: '2023-03-15'
categories: Analysis
tags: WGBS
---

#### Documentation 

- https://github.com/EpiDiverse/snp
- https://github.com/EpiDiverse/snp/blob/master/docs/usage.md (usage for options while running)
- [Emma Strand's Notebook Post](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2023-02-06-EpiDiverse-Bleaching-Pairs-Analysis.md#troubleshooting)

#### Set up 

Make new directory: 

```bash
[kevin_wong1@ssh3 Past_WGBS]$ pwd
/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS
[kevin_wong1@ssh3 Past_WGBS]$ mkdir EpiDiverse
[kevin_wong1@ssh3 Past_WGBS]$ pwd
/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS
[kevin_wong1@ssh3 Past_WGBS]$ cd EpiDiverse/
[kevin_wong1@ssh3 EpiDiverse]$ pwd
/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse
```

Download the freebayes `fasta_generate_regions.py` script:

```
[kevin_wong1@ssh3 EpiDiverse]$ wget https://github.com/freebayes/freebayes/blob/master/scripts/fasta_generate_regions.py
```

Clone the EpiDiverse GitHub:

```
[kevin_wong1@ssh3 EpiDiverse]$ git clone https://github.com/EpiDiverse/snp.git
```

Create the conda environment with mamba. This steo takes ~ 2 hours.

```
[kevin_wong1@ssh3 EpiDiverse]$ interactive
[kevin_wong1@n063 EpiDiverse]$ module load Mamba/22.11.1-4
[kevin_wong1@n063 EpiDiverse]$ mamba env create -f /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/snp/env/environment.yml
```

To view any environments: `conda info --envs`

You may be prompted to exit and re-enter Andromeda, that's OK. Log back in and re-enter the interactive node.

`conda activate snps`

Now run sbatch episnp.sh command to run the desired script.

#### Run jobs in interactive node but be able to switch wifi

1. ssh into HPC.
2. Run the following command: `screen -S <name-of-my-session>` (Replace the (including the < and >) with a descriptive name of your job - don't use spaces in the description.
Run your job(s).
3. Close your computer, go home, switch WiFi networks, whatever. The job will continue running!!
4. To get back to that session, ssh back into the HPC.
5. The, resume the screen session: `screen -r <name-of-my-session>`. If you can't remember the name, you can run screen -list. That will list any running screen sessions.

I named my session epi_TTM.

#### Run EpiDiverse

Run the following script:

```bash
[kevin_wong1@n063 EpiDiverse]$ screen -S epi_TTM
[kevin_wong1@n063 EpiDiverse]$ interactive
[kevin_wong1@n063 EpiDiverse]$ module load Mamba/22.11.1-4
[kevin_wong1@n063 EpiDiverse]$ conda info --envs #double check `snps` is there
[kevin_wong1@n063 EpiDiverse]$ conda activate snps
```

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


trying to resolve container(?) issue: 

```bash
(base) [kevin_wong1@n063 EpiDiverse]$ conda init
no change     /opt/software/Mamba/22.11.1-4/condabin/conda
no change     /opt/software/Mamba/22.11.1-4/bin/conda
no change     /opt/software/Mamba/22.11.1-4/bin/conda-env
no change     /opt/software/Mamba/22.11.1-4/bin/activate
no change     /opt/software/Mamba/22.11.1-4/bin/deactivate
no change     /opt/software/Mamba/22.11.1-4/etc/profile.d/conda.sh
no change     /opt/software/Mamba/22.11.1-4/etc/fish/conf.d/conda.fish
no change     /opt/software/Mamba/22.11.1-4/shell/condabin/Conda.psm1
no change     /opt/software/Mamba/22.11.1-4/shell/condabin/conda-hook.ps1
no change     /opt/software/Mamba/22.11.1-4/lib/python3.10/site-packages/xontrib/conda.xsh
no change     /opt/software/Mamba/22.11.1-4/etc/profile.d/conda.csh
no change     /home/kevin_wong1/.bashrc
No action taken.
```