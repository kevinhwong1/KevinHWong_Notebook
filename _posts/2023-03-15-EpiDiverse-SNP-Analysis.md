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

# 20230315 Attempt

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
Did not work... trying to resolve container(?) issue:

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

# 20230808 Attempt

We are going to try this without activiating conda before running the script. 

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
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed (specific need for my computer)
source /usr/share/Modules/init/sh # load the module function

# load modules needed
echo "START" $(date)
module load Nextflow/20.07.1 #this pipeline requires this version 
module load SAMtools/1.9-foss-2018b 
module load Pysam/0.15.1-foss-2018b-Python-3.6.6

# define location for fasta_generate_regions.py
#fasta_generate_regions.py = ./fasta_generate_regions.py

# only need to direct to input folder not *bam files 
NXF_VER=20.07.1 nextflow run epidiverse/snp -resume \
-profile conda \
--input /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated/ \
--reference /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--output /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/ \
--clusters \
--variants \
--coverage 5 \
--take 47 # Number of samples 

echo "STOP" $(date) # this will output the time it takes to run within the output message

```

```bash
(base) [kevin_wong1@n063 EpiDiverse]$ less episnp.sh_error.273737
Exception in thread "Thread-3" groovy.lang.GroovyRuntimeException: exception while reading process stream
        at org.codehaus.groovy.runtime.ProcessGroovyMethods$TextDumper.run(ProcessGroovyMethods.java:496)
        at java.base/java.lang.Thread.run(Thread.java:834)
Caused by: java.io.IOException: Stream closed
        at java.base/java.io.BufferedInputStream.getBufIfOpen(BufferedInputStream.java:176)
        at java.base/java.io.BufferedInputStream.read1(BufferedInputStream.java:289)
        at java.base/java.io.BufferedInputStream.read(BufferedInputStream.java:351)
        at java.base/sun.nio.cs.StreamDecoder.readBytes(StreamDecoder.java:284)
        at java.base/sun.nio.cs.StreamDecoder.implRead(StreamDecoder.java:326)
        at java.base/sun.nio.cs.StreamDecoder.read(StreamDecoder.java:178)
        at java.base/java.io.InputStreamReader.read(InputStreamReader.java:185)
        at java.base/java.io.BufferedReader.fill(BufferedReader.java:161)
        at java.base/java.io.BufferedReader.readLine(BufferedReader.java:326)
        at java.base/java.io.BufferedReader.readLine(BufferedReader.java:392)
        at org.codehaus.groovy.runtime.ProcessGroovyMethods$TextDumper.run(ProcessGroovyMethods.java:489)
        ... 1 more
```

```bash
[kevin_wong1@ssh3 EpiDiverse]$ less episnp.sh_output.273737
Error executing process > 'SNPS:preprocessing (18-227_S170_L004_R1_001_val_1_bismark_bt2_pe.deduplicated)'

Caused by:
  Failed to create Conda environment
  command: conda env create --prefix /glfs/brick01/gv0/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/work/conda/snps-8882ee7ea1a0aa0094bec65b6ca3edc3 --file /ho
me/kevin_wong1/.nextflow/assets/epidiverse/snp/env/environment.yml
  status : 120
  message:
    ==> WARNING: A newer version of conda exists. <==
      current version: 22.11.1
      latest version: 23.7.2
    
    Please update conda by running
    
        $ conda update -n base -c conda-forge conda
    
    Or to minimize the number of packages updated during conda update use
    
         conda install conda=23.7.2

```

It looks like we need to update or reference a different conda version. Check which conda versions are available on Andromeda

```bash
kevin_wong1@ssh3 EpiDiverse]$ module av -t |& grep -i conda
all/Anaconda3/2020.11
all/Anaconda3/2021.11
all/Anaconda3/2022.05
all/Anaconda3/4.2.0
all/Anaconda3/5.3.0
all/Anaconda3/default
all/Miniconda3/22.11.1-1 # I think it is using this one
all/Miniconda3/4.6.14
all/Miniconda3/4.7.10
all/Miniconda3/4.9.2
lang/Anaconda3/2020.11
lang/Anaconda3/2021.11
lang/Anaconda3/2022.05
lang/Anaconda3/4.2.0
lang/Anaconda3/5.3.0
lang/Miniconda3/22.11.1-1 
lang/Miniconda3/4.6.14
lang/Miniconda3/4.7.10
lang/Miniconda3/4.9.2
Anaconda3/2020.11
Anaconda3/2021.11
Anaconda3/2022.05
Anaconda3/4.2.0
Anaconda3/5.3.0
Anaconda3/default
Miniconda3/22.11.1-1
Miniconda3/4.6.14
Miniconda3/4.7.10
Miniconda3/4.9.2
```

# 20230811 Attempt

Running the same script again but removing all previous auto-generated files. 

Still got the same issue.
```bash
Error executing process > 'SNPS:preprocessing (18-227_S170_L004_R1_001_val_1_bismark_bt2_pe.deduplicated)'

Caused by:
  Failed to create Conda environment
  command: conda env create --prefix /glfs/brick01/gv0/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/work/conda/snps-8882ee7ea1a0aa0094bec65b6ca3edc3 --fil
e /home/kevin_wong1/.nextflow/assets/epidiverse/snp/env/environment.yml
  status : 120
  message:
    ==> WARNING: A newer version of conda exists. <==
      current version: 22.11.1
      latest version: 23.7.2
    
    Please update conda by running
    
        $ conda update -n base -c conda-forge conda
    
    Or to minimize the number of packages updated during conda update use
    
         conda install conda=23.7.2
```

I am going to try this with profile singularity instead


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
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed (specific need for my computer)
source /usr/share/Modules/init/sh # load the module function

# load modules needed
echo "START" $(date)
module load Nextflow/20.07.1 #this pipeline requires this version 
module load SAMtools/1.9-foss-2018b 
module load Pysam/0.15.1-foss-2018b-Python-3.6.6

# define location for fasta_generate_regions.py
#fasta_generate_regions.py = ./fasta_generate_regions.py

# only need to direct to input folder not *bam files 
NXF_VER=20.07.1 nextflow run epidiverse/snp -resume \
-profile singularity \
--input /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated/ \
--reference /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--output /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/ \
--clusters \
--variants \
--coverage 5 \
--take 47 # Number of samples 

echo "STOP" $(date) # this will output the time it takes to run within the output message

```

Error in output file: 

```bash
Error executing process > 'SNPS:preprocessing (18-346_S193_L004_R1_001_val_1_bismark_bt2_pe.deduplicated)'

Caused by:
  Process `SNPS:preprocessing (18-346_S193_L004_R1_001_val_1_bismark_bt2_pe.deduplicated)` terminated with an error exit status (1)

Command executed:

  samtools sort -T deleteme -m 966367642 -@ 4 \
  -o sorted.bam 18-346_S193_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam || exit $?
  samtools calmd -b sorted.bam past_filtered_assembly.fasta 1> calmd.bam 2> /dev/null && rm sorted.bam
  samtools index calmd.bam

Command exit status:
  1

Command output:
  (empty)

Command error:
  INFO:    Environment variable SINGULARITYENV_TMP is set, but APPTAINERENV_TMP is preferred
  INFO:    Environment variable SINGULARITYENV_TMPDIR is set, but APPTAINERENV_TMPDIR is preferred
  INFO:    Environment variable SINGULARITYENV_NXF_DEBUG is set, but APPTAINERENV_NXF_DEBUG is preferred
  [E::hts_open_format] Failed to open file "18-346_S193_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam" : No such file or directory
  samtools sort: can't open "18-346_S193_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bam": No such file or directory

Work dir:
  /glfs/brick01/gv0/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/work/bd/1bc1be4895ef173a94c9bbad215226

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

# 20230814 Attempt

As per the epidivere instructions, I need to create a folder for each sample.bam 

First I will move the sorted bam files:

```bash
mkdir test_epi
cd ../bismark_deduplicated/
cp *.deduplicated_sorted.bam ../test_epi
```

Second, I will loop to make folders and move the correct sample into the folder:

```bash
for x in ./*.bam; do
  mkdir "${x%.*}" && mv "$x" "${x%.*}"
done
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
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed (specific need for my computer)
source /usr/share/Modules/init/sh # load the module function

# load modules needed
echo "START" $(date)
module load Nextflow/20.07.1 #this pipeline requires this version 
module load SAMtools/1.9-foss-2018b 
module load Pysam/0.15.1-foss-2018b-Python-3.6.6

# define location for fasta_generate_regions.py
#fasta_generate_regions.py = ./fasta_generate_regions.py

# only need to direct to input folder not *bam files 
NXF_VER=20.07.1 nextflow run epidiverse/snp -resume \
-profile singularity \
--input /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/test_epidiverse/ \
--reference /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--output /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/ \
--clusters \
--variants \
--coverage 5 \
--take 47 # Number of samples 

echo "STOP" $(date) # this will output the time it takes to run within the output message

```

This produced an immediate error. I will try this on the non-sorted bam files. 


# 20230815 Attempt

```bash
cd /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq
mkdir test_epidiverse2
cd bismark_deduplicated
```

`nano epidiverse_prep`

```bash
#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=18
#SBATCH --export=NONE
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# copy files 
cp *R1_001_val_1_bismark_bt2_pe.deduplicated.bam ../test_epidiverse2
```

Now make folders for each file: 

```bash
cd ../test_epidiverse2

for x in ./*.bam; do
  mkdir "${x%.*}" && mv "$x" "${x%.*}"
done
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
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed (specific need for my computer)
source /usr/share/Modules/init/sh # load the module function

# load modules needed
echo "START" $(date)
module load Nextflow/20.07.1 #this pipeline requires this version 
module load SAMtools/1.9-foss-2018b 
module load Pysam/0.15.1-foss-2018b-Python-3.6.6

# define location for fasta_generate_regions.py
#fasta_generate_regions.py = ./fasta_generate_regions.py

# only need to direct to input folder not *bam files 
NXF_VER=20.07.1 nextflow run epidiverse/snp -resume \
-profile singularity \
--input /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/test_epidiverse2/ \
--reference /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--output /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/ \
--clusters \
--variants \
--coverage 5 \
--take 47 # Number of samples 

echo "STOP" $(date) # this will output the time it takes to run within the output message

```

Immediately stops with this error: 

```
ERROR: cannot find valid *.bam files in dir: /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/test_epidiverse2/
```

Maybe they don't have to be in the different folders? Will try again with all the sample files in one folder

# 20230815 Attempt 2


```bash
cd /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq
mkdir test_epidiverse
cd bismark_deduplicated
```

`nano epidiverse_prep.sh`

```bash
#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=18
#SBATCH --export=NONE
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# copy files 
cp *.deduplicated_sorted.bam ../test_epidiverse
```

`cd `

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
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed (specific need for my computer)
source /usr/share/Modules/init/sh # load the module function

# load modules needed
echo "START" $(date)
module load Nextflow/20.07.1 #this pipeline requires this version 
module load SAMtools/1.9-foss-2018b 
module load Pysam/0.15.1-foss-2018b-Python-3.6.6

# define location for fasta_generate_regions.py
#fasta_generate_regions.py = ./fasta_generate_regions.py

# only need to direct to input folder not *bam files 
NXF_VER=20.07.1 nextflow run epidiverse/snp -resume \
-profile singularity \
--input /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/test_epidiverse/ \
--reference /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--output /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/ \
--clusters \
--variants \
--coverage 5 \
--take 47 # Number of samples 

echo "STOP" $(date) # this will output the time it takes to run within the output message

```

Got this error again: 

```bash
Error executing process > 'SNPS:preprocessing (18-67_S176.deduplicated_sorted)'

Caused by:
  Process `SNPS:preprocessing (18-67_S176.deduplicated_sorted)` terminated with an error exit status (1)

Command executed:

  samtools sort -T deleteme -m 966367642 -@ 4 \
  -o sorted.bam 18-67_S176.deduplicated_sorted.bam || exit $?
  samtools calmd -b sorted.bam past_filtered_assembly.fasta 1> calmd.bam 2> /dev/null && rm sorted.bam
  samtools index calmd.bam

Command exit status:
  1

Command output:
  (empty)

Command error:
  INFO:    Environment variable SINGULARITYENV_TMP is set, but APPTAINERENV_TMP is preferred
  INFO:    Environment variable SINGULARITYENV_TMPDIR is set, but APPTAINERENV_TMPDIR is preferred
  INFO:    Environment variable SINGULARITYENV_NXF_DEBUG is set, but APPTAINERENV_NXF_DEBUG is preferred
  [E::hts_open_format] Failed to open file "18-67_S176.deduplicated_sorted.bam" : No such file or directory
  samtools sort: can't open "18-67_S176.deduplicated_sorted.bam": No such file or directory

Work dir:
  /glfs/brick01/gv0/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/work/5d/e6bdcee831ab4508170d1317d6025f
```

# 20230816 Attempt 1

Kevin Bryan suggested this: 

```
Can you try putting this in your batch file before calling nextflow?

`export APPTAINER_BINDPATH=/data,/glfs`

If that doesn’t work, try referring to all of your paths as starting with /glfs/brick01/gv0 instead of /data. You can see that the “files” that it said it couldn’t find are symlinks to the /glfs path, but because the singularity container only binds what it thinks the working directory is, it binds it with /data instead of /glfs, but nextflow has created the symlinks using the canonical path (/glfs before launching the container), so the links aren’t valid.
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
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed (specific need for my computer)
source /usr/share/Modules/init/sh # load the module function

#Forcing the glfs path
export APPTAINER_BINDPATH=/data,/glfs

# load modules needed
echo "START" $(date)
module load Nextflow/20.07.1 #this pipeline requires this version 
module load SAMtools/1.9-foss-2018b 
module load Pysam/0.15.1-foss-2018b-Python-3.6.6

# define location for fasta_generate_regions.py
#fasta_generate_regions.py = ./fasta_generate_regions.py

# only need to direct to input folder not *bam files 
NXF_VER=20.07.1 nextflow run epidiverse/snp -resume \
-profile singularity \
--input /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/test_epidiverse/ \
--reference /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--output /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/ \
--clusters \
--variants \
--coverage 5 \
--take 47 # Number of samples 

echo "STOP" $(date) # this will output the time it takes to run within the output message

```

This produced the same error. Let me try putting the full glfs path

# 20230816 Attempt 2

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
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed (specific need for my computer)
source /usr/share/Modules/init/sh # load the module function

#Forcing the glfs path
export APPTAINER_BINDPATH=/data,/glfs/brick01/gv0

# load modules needed
echo "START" $(date)
module load Nextflow/20.07.1 #this pipeline requires this version 
module load SAMtools/1.9-foss-2018b 
module load Pysam/0.15.1-foss-2018b-Python-3.6.6

# define location for fasta_generate_regions.py
#fasta_generate_regions.py = ./fasta_generate_regions.py

# only need to direct to input folder not *bam files 
NXF_VER=20.07.1 nextflow run epidiverse/snp -resume \
-profile singularity \
--input /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/test_epidiverse/ \
--reference /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--output /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/ \
--clusters \
--variants \
--coverage 5 \
--take 47 # Number of samples 

echo "STOP" $(date) # this will output the time it takes to run within the output message

```

This also did not work. Let me try replacing all of the paths as Kevin suggested: 

# 20230816 Attempt 3

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
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed (specific need for my computer)
source /usr/share/Modules/init/sh # load the module function

# load modules needed
echo "START" $(date)
module load Nextflow/20.07.1 #this pipeline requires this version 
module load SAMtools/1.9-foss-2018b 
module load Pysam/0.15.1-foss-2018b-Python-3.6.6

# define location for fasta_generate_regions.py
#fasta_generate_regions.py = ./fasta_generate_regions.py

# only need to direct to input folder not *bam files 
NXF_VER=20.07.1 nextflow run epidiverse/snp -resume \
-profile singularity \
--input /glfs/brick01/gv0/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/test_epidiverse/ \
--reference /glfs/brick01/gv0/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--output /glfs/brick01/gv0/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/ \
--clusters \
--variants \
--coverage 5 \
--take 47 # Number of samples 

echo "STOP" $(date) # this will output the time it takes to run within the output message

```

Still getting the same error.... 

```bash
Error executing process > 'SNPS:preprocessing (18-67_S176.deduplicated_sorted)'

Caused by:
  Process `SNPS:preprocessing (18-67_S176.deduplicated_sorted)` terminated with an error exit status (1)

Command executed:

  samtools sort -T deleteme -m 966367642 -@ 4 \
  -o sorted.bam 18-67_S176.deduplicated_sorted.bam || exit $?
  samtools calmd -b sorted.bam past_filtered_assembly.fasta 1> calmd.bam 2> /dev/null && rm sorted.bam
  samtools index calmd.bam

Command exit status:
  1

Command output:
  (empty)

Command error:
  INFO:    Environment variable SINGULARITYENV_TMP is set, but APPTAINERENV_TMP is preferred
  INFO:    Environment variable SINGULARITYENV_TMPDIR is set, but APPTAINERENV_TMPDIR is preferred
  INFO:    Environment variable SINGULARITYENV_NXF_DEBUG is set, but APPTAINERENV_NXF_DEBUG is preferred
  [E::hts_open_format] Failed to open file "18-67_S176.deduplicated_sorted.bam" : No such file or directory
  samtools sort: can't open "18-67_S176.deduplicated_sorted.bam": No such file or directory

Work dir:
  /glfs/brick01/gv0/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/work/cf/d934a34af1e0ccf6bf52ac672daffa
```

# 20230816 Attempt 4

Next I will try Kevin's other suggestion of modifying the config file: 

`nano data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/snp/assets/custom.config`

```bash
process {

        executor = 'pbspro'

        // with conda
        module = ['Miniconda3']
        conda = "${baseDir}/env/environment.yml"

        // with docker/singularity
        container = "epidiverse/dmr"
        containerOptions '--volume /glfs:/glfs' #adding this line here
```

# 20230816 

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
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed (specific need for my computer)
source /usr/share/Modules/init/sh # load the module function

# load modules needed
echo "START" $(date)
module load Nextflow/20.07.1 #this pipeline requires this version 
module load SAMtools/1.9-foss-2018b 
module load Pysam/0.15.1-foss-2018b-Python-3.6.6

# define location for fasta_generate_regions.py
#fasta_generate_regions.py = ./fasta_generate_regions.py

# only need to direct to input folder not *bam files 
NXF_VER=20.07.1 nextflow run epidiverse/snp -resume \
-profile singularity \
--input /glfs/brick01/gv0/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/test_epidiverse/ \
--reference /glfs/brick01/gv0/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--output /glfs/brick01/gv0/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/ \
--clusters \
--variants \
--coverage 5 \
--take 47 # Number of samples 

echo "STOP" $(date) # this will output the time it takes to run within the output message

```

Got the same error :(

```bash
Error executing process > 'SNPS:preprocessing (18-202_S188.deduplicated_sorted)'

Caused by:
  Process `SNPS:preprocessing (18-202_S188.deduplicated_sorted)` terminated with an error exit status (1)

Command executed:

  samtools sort -T deleteme -m 966367642 -@ 4 \
  -o sorted.bam 18-202_S188.deduplicated_sorted.bam || exit $?
  samtools calmd -b sorted.bam past_filtered_assembly.fasta 1> calmd.bam 2> /dev/null && rm sorted.bam
  samtools index calmd.bam

Command exit status:
  1

Command output:
  (empty)

Command error:
  INFO:    Environment variable SINGULARITYENV_TMP is set, but APPTAINERENV_TMP is preferred
  INFO:    Environment variable SINGULARITYENV_TMPDIR is set, but APPTAINERENV_TMPDIR is preferred
  INFO:    Environment variable SINGULARITYENV_NXF_DEBUG is set, but APPTAINERENV_NXF_DEBUG is preferred
  [E::hts_open_format] Failed to open file "18-202_S188.deduplicated_sorted.bam" : No such file or directory
  samtools sort: can't open "18-202_S188.deduplicated_sorted.bam": No such file or directory

Work dir:
  /glfs/brick01/gv0/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/work/c8/e2a9691fd0b2b35e2b03bbf311a5dd

Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line
```

# 20230821 Attempt

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
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed (specific need for my computer)
#source /usr/share/Modules/init/sh # load the module function

# load modules needed
echo "START" $(date)
module load Anaconda3/2022.05
module load Nextflow/20.07.1 #this pipeline requires this version 
#module load SAMtools/1.9-foss-2018b 
#module load Pysam/0.15.1-foss-2018b-Python-3.6.6

# define location for fasta_generate_regions.py
#fasta_generate_regions.py = ./fasta_generate_regions.py

#make conda env

conda env create --prefix /glfs/brick01/gv0/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/work/conda/snps-8882ee7ea1a0aa0094bec65b6ca3edc2 --file /home/kevin_wong1/.nextflow/assets/epidiverse/snp/env/environment.yml --force

conda activate snps

# only need to direct to input folder not *bam files 
NXF_VER=20.07.1 nextflow run epidiverse/snp -resume \
-profile conda \
--input /glfs/brick01/gv0/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/test_epidiverse/ \
--reference /glfs/brick01/gv0/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--output /glfs/brick01/gv0/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/EpiDiverse/ \
--clusters \
--variants \
--coverage 5 \
--take 47 # Number of samples 

echo "STOP" $(date) # this will output the time it takes to run within the output message

```

