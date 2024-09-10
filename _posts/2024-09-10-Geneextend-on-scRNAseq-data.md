---
layout: post
title: Geneextend on scRNAseq data
date: '2024-09-10'
categories: Analysis
tags: scRNAseq
---

# GeneExt

[GeneExt](https://github.com/sebepedroslab/GeneExt/tree/main) is a program from the Sebé-Pedrós Lab which has a great manual that explains this probem with 3' biased sequencing approaches in non-model organism in great detail and with helpful nuance. The program allows the user to create a modified GTF/GFF file based on a BAM file of mapping data. It appends 3' UTRs to the genes in a way that is informed by the mapping of reads, so the length of the 3' UTRs added can change dynamically for each gene based on where reads are acutally mapping.

See manual [here](https://github.com/sebepedroslab/GeneExt/blob/main/Manual.md)

### Installing GeneExt on Unity

Following Unity guidance on conda environments: https://docs.unity.rc.umass.edu/documentation/software/conda/

```
cd /work/pi_hputnam_uri_edu
mkdir -p conda/envs
cd conda/envs

salloc -c 6 -p cpu #request space on cluster, non-login node

module load miniconda/22.11.1-1 #load miniconda

git clone https://github.com/sebepedroslab/GeneExt.git #clone GeneExt repo

cd GeneExt

# create environment
conda env create --prefix /work/pi_hputnam_uri_edu/conda/envs/GeneExt/geneext -f environment.yaml

# OPTIONAL, make a symlink (shortcut) to home directory
ln -s /work/pi_hputnam_uri_edu/conda/envs/GeneExt/geneext ~/geneext

# activate environment
conda activate /work/pi_hputnam_uri_edu/conda/envs/GeneExt/geneext  #wow that environment name is very long and may be annoying

# Test run 
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
Median extension length: 965.0 bp
```

### Yay! A successful installation!

### Running GeneExt: Run this on the bam file from cellranger count, with the GTF file you used for cellranger mkref

```
cd /work/pi_hputnam_uri_edu/snRNA_analysis
mkdir GeneExt

cd /work/pi_hputnam_uri_edu/snRNA_analysis/scripts
nano GeneExt_POC_L001.sh
```

```
#!/bin/bash
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=24
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --time 06:00:00
#SBATCH --error="%x_error.%j" 
#SBATCH --output="%x_output.%j" 
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails

module load miniconda/22.11.1-1 #load miniconda

# activate environment
conda activate /work/pi_hputnam_uri_edu/conda/envs/GeneExt/geneext

cd /work/pi_hputnam_uri_edu/snRNA_analysis/GeneExt

#the gtf has to have the gene and transcript IDs be different.

# Append -T to each transcript_id and save to a new file
sed 's/transcript_id "\([^"]*\)"/transcript_id "\1-T"/g' /work/pi_hputnam_uri_edu/snRNA_analysis/references/Pocillopora_acuta_HIv2.gtf > /work/pi_hputnam_uri_edu/snRNA_analysis/references/Pocillopora_acuta_HIv2_modified.gtf

# use --clip_strand both to not allow GeneExt to create overlaps on the same strand

python /work/pi_hputnam_uri_edu/conda/envs/GeneExt/geneext.py \
    -g /work/pi_hputnam_uri_edu/snRNA_analysis/references/Pocillopora_acuta_HIv2_modified.gtf \
    -b /scratch/workspace/zdellaert_uri_edu-shared/count_POC_test/outs/possorted_genome_bam.bam \
    -o Pocillopora_acuta_GeneExt.gtf \
    -j 24 -v 3 --clip_strand both
```

Submitted batch job 23863089

```
      ____                 _____      _
     / ___| ___ _ __   ___| ____|_  _| |_
    | |  _ / _ \ '_ \ / _ \  _| \ \/ / __|
    | |_| |  __/ | | |  __/ |___ >  <| |_
     \____|\___|_| |_|\___|_____/_/\_\__|

          ______    ___    ______
    -----[______]==[___]==[______]===>----

    Gene model adjustment for improved single-cell RNA-seq data counting


python geneext.py -genome /work/pi_hputnam_uri_edu/snRNA_analysis/references/Pocillopora_acuta_HIv2_modified.gtf -bam /scratch/workspace/zdellaert_uri_edu-shared/count_POC_test/outs/possorted_genome_bam.bam -output Pocillopora_acuta_GeneExt.gtf -tag GeneExt -verbose 3 -jobs 24 -output_mode new_transcript -clip_strand both -peak_perc 25
╭──────────────────╮
│ Preflight checks │
╰──────────────────╯
Temporary directory isn't set, setting to tmp/
Temporary directory created: tmp/
Alignment file ... OK
Genome annotation file .... OK
Input: /work/pi_hputnam_uri_edu/snRNA_analysis/references/Pocillopora_acuta_HIv2_modified.gtf, guessed format: gtf
Output: Pocillopora_acuta_GeneExt.gtf, guessed format: gtf
Checking gene exons...
Genome annotation warning: Could not find "gene" features in
/work/pi_hputnam_uri_edu/snRNA_analysis/references/Pocillopora_acuta_HIv2_modifi
ed.gtf! Trying to fix ...
Loading the database ...
```

lots of other output...

DONE!

```
	Extended genes written: Pocillopora_acuta_GeneExt.gtf
done
╭───────────╮
│ All done! │
╰───────────╯
Extended **14284/33730 genes**
Median extension length: 1622.0 bp
Running:
	Rscript geneext/plot_extensions.R tmp/extensions.tsv Pocillopora_acuta_GeneExt.gtf.extension_length.pdf
Running:
	Rscript geneext/peak_density.R tmp/genic_peaks.bed tmp/allpeaks_noov.bed Pocillopora_acuta_GeneExt.gtf.peak_coverage.pdf 25
Removing tmp/_genes_peaks_closest
Removing tmp/_genes_tmp
Removing tmp/_genes_tmp_sorted
Removing tmp/_peaks_tmp
Removing tmp/_peaks_tmp_sorted
Removing tmp/minus.bam
Removing tmp/plus.bam
```

Copy file to references directory:

```
cd /work/pi_hputnam_uri_edu/snRNA_analysis/GeneExt
cp Pocillopora_acuta_GeneExt.gtf ../references/
```

### Rerun cellranger mkref with GeneExt GTF file

```
cd /work/pi_hputnam_uri_edu/snRNA_analysis/scripts
nano mkref_Pacuta_GeneExt.sh
```


```
#!/usr/bin/env bash
#
#
# =============================================================================
# Job Script
# =============================================================================
#
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=24
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=250GB
#SBATCH --error="%x_error.%j" 
#SBATCH --output="%x_output.%j" 
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails

# Usage: cellranger mkref [OPTIONS] --genome <GENOME_NAMES> --fasta <FASTA_FILES> --genes <GTF_FILES>

cd /work/pi_hputnam_uri_edu/snRNA_analysis/references

cellranger mkref --genome=Pocillopora_acuta_GeneExt --fasta=Pocillopora_acuta_HIv2.assembly.fasta  --genes=Pocillopora_acuta_GeneExt.gtf --nthreads=$SLURM_CPUS_ON_NODE --jobmode=local --localcores=24 --localmem=225
```

```
sbatch mkref_Pacuta.sh
```

### Rerun cellranger on POC_L001 sample with modified custom reference

```
cd /work/pi_hputnam_uri_edu/snRNA_analysis/scripts
nano count_POC_test_GeneExt.sh
```

```
#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=24
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=250GB
#SBATCH --time 06:00:00
#SBATCH --error="%x_error.%j" 
#SBATCH --output="%x_output.%j" 
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails

# Usage: cellranger count [OPTIONS] --id <ID> --transcriptome <PATH> --create-bam <true|false>

cd /scratch/workspace/zdellaert_uri_edu-shared/

cellranger count --id=count_POC_test_GeneExt \
           --transcriptome=/work/pi_hputnam_uri_edu/snRNA_analysis/references/Pocillopora_acuta_GeneExt \
           --fastqs=/scratch/workspace/zdellaert_uri_edu-shared/POC_L001 \
           --sample=TD006822-POC-2 \
           --create-bam=true \
           --jobmode=local \
           --localcores=24 \
           --localmem=225
```
