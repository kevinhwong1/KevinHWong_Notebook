---
layout: post
title: Pegasus conda environments
date: '2024-10-30'
categories: Analysis
tags: Pegasus, HPC, UM
---

# Anaconda

```bash
# Download the Anaconda installer
cd /nethome/kxw755/
wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh

# Execute the anaconda installer
bash Anaconda3-2021.05-Linux-x86_64.sh

# Activate conda with the new Anaconda3 folder in your home directory (Depending on your download <path to conda> might be named ‘ENTER’ or 'Anaconda3')
source <path to conda>/bin/activate
conda init

# 20241030 I updated by conda to resolve channel issues:
conda install conda=24.5.0

# When logging into Pegasus, I have to activate my conda each time 
cd nethome/kxw755/
source anaconda3/bin/activate 
```

# SAMTools

```bash
cd /nethome/kxw755/
conda create -n samtools_env
conda activate samtools_env
conda install -c bioconda samtools
samtools --version
```

output
```
samtools 1.3.1
Using htslib 1.3.1
Copyright (C) 2016 Genome Research Ltd.
```

# STAR and MultiQC

```bash
cd /nethome/kxw755/
conda create -n star_env
conda activate star_env
conda install -c bioconda star=2.7.11b
STAR --version #2.7.11b

conda install -c bioconda multiqc
multiqc --version #1.25.1
```

# AGAT

```bash
cd /nethome/kxw755/
conda create -n agat
conda activate agat
conda install agat
```

# NCBI datasets

```bash
cd /nethome/kxw755/
conda create -n ncbi_datasets python=3.8
conda activate ncbi_datasets
conda install -c conda-forge -c bioconda ncbi-datasets-cli
```

# GeneExt

[GeneExt](https://github.com/sebepedroslab/GeneExt/tree/main) is a program from the Sebé-Pedrós Lab which has a great manual that explains this probem with 3' biased sequencing approaches in non-model organism in great detail and with helpful nuance. The program allows the user to create a modified GTF/GFF file based on a BAM file of mapping data. It appends 3' UTRs to the genes in a way that is informed by the mapping of reads, so the length of the 3' UTRs added can change dynamically for each gene based on where reads are acutally mapping.

See manual [here](https://github.com/sebepedroslab/GeneExt/blob/main/Manual.md)

The code from this post was inspired by Zoe Dellaert. 

### Installing GeneExt on Pegasus

```bash
cd nethome/kxw755/
source anaconda3/bin/activate 

# Copy repo
git clone https://github.com/sebepedroslab/GeneExt.git

# create environment
conda env create --prefix geneext -f GeneExt/environment.yaml

#activate conda
conda activate /projectnb/pegasus/nethome/kxw755/geneext

# Test run 
cd GeneExt
python geneext.py -g test_data/annotation.gtf -b test_data/alignments.bam -o result.gtf --peak_perc 0
```

output:

```


      ____                 _____      _
     / ___| ___ _ __   ___| ____|_  _| |_
    | |  _ / _ \ '_ \ / _ \  _| \ \/ / __|
    | |_| |  __/ | | |  __/ |___ >  <| |_
     \____|\___|_| |_|\___|_____/_/\_\__|

          ______    ___    ______
    -----[______]==[___]==[______]===>----

    Gene model adjustment for improved single-cell RNA-seq data counting


╭──────────────────╮
│ Preflight checks │
╰──────────────────╯
Genome annotation warning: Could not find "gene" features in test_data/annotation.gtf! Trying to fix ...
Running macs2 ... 
########## macs2 FAILED ##############
return code:  1 
Output:  Traceback (most recent call last):
  File "/nethome/kxw755/conda/envs/GeneExt/geneext/bin/macs2", line 653, in <module>
    main()
  File "/nethome/kxw755/conda/envs/GeneExt/geneext/bin/macs2", line 49, in main
    from MACS2.callpeak_cmd import run
  File "/nethome/kxw755/conda/envs/GeneExt/geneext/lib/python3.9/site-packages/MACS2/callpeak_cmd.py", line 23, in <module>
    from MACS2.OptValidator import opt_validate
  File "/nethome/kxw755/conda/envs/GeneExt/geneext/lib/python3.9/site-packages/MACS2/OptValidator.py", line 20, in <module>
    from MACS2.IO.Parser import BEDParser, ELANDResultParser, ELANDMultiParser, \
  File "__init__.pxd", line 206, in init MACS2.IO.Parser
ValueError: numpy.dtype size changed, may indicate binary incompatibility. Expected 96 from C header, got 88 from PyObject

Traceback (most recent call last):
  File "/projectnb/pegasus/nethome/kxw755/conda/envs/GeneExt/geneext.py", line 703, in <module>
    helper.run_macs2(tempdir+'/' + 'plus.bam','plus',tempdir,verbose = verbose)
  File "/projectnb/pegasus/nethome/kxw755/conda/envs/GeneExt/geneext/helper.py", line 79, in run_macs2
    ps.check_returncode()
  File "/nethome/kxw755/conda/envs/GeneExt/geneext/lib/python3.9/subprocess.py", line 460, in check_returncode
    raise CalledProcessError(self.returncode, self.args, self.stdout,
subprocess.CalledProcessError: Command '('macs2', 'callpeak', '-t', 'tmp/plus.bam', '-f', 'BAM', '--keep-dup', '20', '-q', '0.01', '--shift', '1', '--extsize', '100', '--broad', '--nomodel', '--min-length', '30', '-n', 'plus', '--outdir', 'tmp')' returned non-zero exit status 1.
```

Okay so Jill was getting this error too... lets try to fix this error. 

```bash
# Reinstall numpy with a Compatible Version:
conda update numpy

# Verify version of numpy
python -c "import numpy; print(numpy.__version__)" #1.26.4

#Reinstall macs2 to Ensure Compatibility:
conda install -c bioconda macs2=2.2.7.1

# Update Pandas
conda update pandas

# check to see if they are compatable
python -c "import pandas as pd; import numpy as np; print(pd.__version__, np.__version__)" # 2.2.2 1.26.4
```

Lets try running it again

```
rm -r tmp
python geneext.py -g test_data/annotation.gtf -b test_data/alignments.bam -o result.gtf --peak_perc 0
```

Output
```

      ____                 _____      _   
     / ___| ___ _ __   ___| ____|_  _| |_ 
    | |  _ / _ \ '_ \ / _ \  _| \ \/ / __|
    | |_| |  __/ | | |  __/ |___ >  <| |_ 
     \____|\___|_| |_|\___|_____/_/\_\__|
     
          ______    ___    ______    
    -----[______]==[___]==[______]===>----

    Gene model adjustment for improved single-cell RNA-seq data counting


╭──────────────────╮
│ Preflight checks │
╰──────────────────╯
Genome annotation warning: Could not find "gene" features in test_data/annotation.gtf! Trying to fix ...
╭───────────╮
│ Execution │
╰───────────╯
Running macs2 ... done
Filtering macs2 peaks ... done
Extending genes ... done
╭───────────╮
│ All done! │
╰───────────╯
Extended 23/35 genes
Median extension length: 931.0 bp
```
