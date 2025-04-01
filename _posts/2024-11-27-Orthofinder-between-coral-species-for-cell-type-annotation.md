---
layout: post
title: Orthofinder between coral species for cell type annotation
date: '2024-11-27'
categories: Analysis
tags: scRNAseq, orthofinder
---

# Orthofinder

| **Comp_Name** | **Comp_Number** | **Speices_1** | **Species_2** | **Speices_3** | **Species_4** | **Speices_5** | **Species_6** | **Species_7** |
|:-------------:|:---------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|
|   Gfas_Spis   |        1        |      spis     |      gfas     |               |               |               |               |               |
|   Pdam_Spis   |        2        |      spis     |      pdam     |               |               |               |               |               |
|   Acer_Spis   |        3        |      spis     |      acer     |               |               |               |               |               |
|  Coral_ortho  |        4        |      spis     |      gfas     |      pdam     |      acer     |      oarb     |               |               |
|   Cnid_ortho  |        5        |      spis     |      gfas     |      pdam     |      acer     |      oarb     |      nvec     |     xenia     |
|    Stemcell   |        6        |      nvec     |      pdam     |      acer     |               |               |               |               |


```bash
source anaconda3/bin/activate 

conda create -n orthofinder_env -y -c bioconda orthofinder diamond fasttree
conda activate orthofinder_env

cd /scratch/projects/dark_genes
mkdir orthofinder_comps
```

Paths:

```bash
Spis: /nethome/kxw755/genomes/GCA_002571385.2_Spis/spis_protein.faa
Gfas: /nethome/kxw755/genomes/GCA_948470475.1_Gfas/data/GCA_948470475.1/gfas_1.0.proteins.fasta
Acer: /nethome/kxw755/genomes/GCA_032359415.1_Acer/acer_protein.faa
Pdam: /nethome/kxw755/genomes/pdam_genome/pdam_proteins.fasta
Apoc: /nethome/kxw755/genomes/Apoc_14110456/apoculata_proteins.fasta.gz
Xenia: /nethome/kxw755/genomes/GCF_021976095.1_Xenia/xenia_protein.faa
Nvec: /nethome/kxw755/genomes/GCF_932526225.1_Nvec/nvec_protein.faa
```

## Gfas_markers

```bash
mkdir Gfas_Spis
cd Gfas_Spis

ln -s /nethome/kxw755/genomes/GCA_948470475.1_Gfas/data/GCA_948470475.1/gfas_1.0.proteins.fasta
ln -s /nethome/kxw755/genomes/GCA_002571385.2_Spis/spis_protein.faa

mkdir results
nano of_gfas_spis.job
```

```bash
#!/bin/bash

#BSUB -J OF_Gfas_Spis
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o OF_Gfas_Spis_%J.out
#BSUB -e OF_Gfas_Spis_%J.err
#BSUB -B
#BSUB -N
###################################################################

# Parameters
GENOME_DIR="/scratch/projects/dark_genes/orthofinder_comps/Gfas_Spis"  # Update this with the path to your genome/proteome files
OUTPUT_DIR="/scratch/projects/dark_genes/orthofinder_comps/Gfas_Spis/results"  # Update this with the desired output directory
THREADS=4  # Adjust based on the number of available CPU cores

# Run OrthoFinder
echo "Running OrthoFinder..."
orthofinder -f "$GENOME_DIR" -o "$OUTPUT_DIR" -t $THREADS -a $THREADS

# Explanation of parameters:
# -f : Path to the folder containing genome/protein files
# -o : Path to the output folder
# -t : Number of threads for OrthoFinder to use
# -a : Number of threads for DIAMOND alignment

echo "OrthoFinder analysis complete. Results saved to $OUTPUT_DIR."
```


## Pdam_markers

```bash
mkdir Pdam_Spis
cd Pdam_Spis

ln -s /nethome/kxw755/genomes/pdam_genome/pdam_proteins.fasta
ln -s /nethome/kxw755/genomes/GCA_002571385.2_Spis/spis_protein.faa

mkdir results
nano of_pdam_spis.job
```

```bash
#!/bin/bash

#BSUB -J OF_Pdam_Spis
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o OF_Pdam_Spis_%J.out
#BSUB -e OF_Pdam_Spis_%J.err
#BSUB -B
#BSUB -N
###################################################################

# Parameters
GENOME_DIR="/scratch/projects/dark_genes/orthofinder_comps/Pdam_Spis"  # Update this with the path to your genome/proteome files
OUTPUT_DIR="/scratch/projects/dark_genes/orthofinder_comps/Pdam_Spis/results"  # Update this with the desired output directory
THREADS=8  # Adjust based on the number of available CPU cores

# Run OrthoFinder
echo "Running OrthoFinder..."
orthofinder -f "$GENOME_DIR" -o "$OUTPUT_DIR" -t $THREADS -a $THREADS

# Explanation of parameters:
# -f : Path to the folder containing genome/protein files
# -o : Path to the output folder
# -t : Number of threads for OrthoFinder to use
# -a : Number of threads for DIAMOND alignment

echo "OrthoFinder analysis complete. Results saved to $OUTPUT_DIR."
```

`bsub < of_pdam_spis.job`


## Acer_markers

```bash
mkdir Acer_Spis
cd Acer_Spis

ln -s /nethome/kxw755/genomes/GCA_032359415.1_Acer/acer_protein.faa
ln -s /nethome/kxw755/genomes/GCA_002571385.2_Spis/spis_protein.faa

nano of_acer_spis.job
```

```bash
#!/bin/bash

#BSUB -J OF_Acer_Spis
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o OF_Acer_Spis_%J.out
#BSUB -e OF_Acer_Spis_%J.err
#BSUB -B
#BSUB -N
###################################################################

# Parameters
GENOME_DIR="/scratch/projects/dark_genes/orthofinder_comps/Acer_Spis"  # Update this with the path to your genome/proteome files
OUTPUT_DIR="/scratch/projects/dark_genes/orthofinder_comps/Acer_Spis/results"  # Update this with the desired output directory
THREADS=8  # Adjust based on the number of available CPU cores

# Run OrthoFinder
echo "Running OrthoFinder..."
orthofinder -f "$GENOME_DIR" -o "$OUTPUT_DIR" -t $THREADS -a $THREADS

# Explanation of parameters:
# -f : Path to the folder containing genome/protein files
# -o : Path to the output folder
# -t : Number of threads for OrthoFinder to use
# -a : Number of threads for DIAMOND alignment

echo "OrthoFinder analysis complete. Results saved to $OUTPUT_DIR."
```

`bsub < of_acer_spis.job`


## Coral_ortho

```bash
mkdir Coral_ortho
cd Coral_ortho

ln -s /nethome/kxw755/genomes/GCA_032359415.1_Acer/acer_protein.faa
ln -s /nethome/kxw755/genomes/GCA_002571385.2_Spis/spis_protein.faa
ln -s /nethome/kxw755/genomes/GCA_948470475.1_Gfas/data/GCA_948470475.1/gfas_1.0.proteins.fasta
ln -s /nethome/kxw755/genomes/pdam_genome/pdam_proteins.fasta
ln -s /nethome/kxw755/genomes/Apoc_14110456/apoculata_proteins.fasta

nano of_coral_ortho.job
```

```bash
#!/bin/bash

#BSUB -J OF_Coral_ortho
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o OF_Coral_ortho_%J.out
#BSUB -e OF_Coral_ortho_%J.err
#BSUB -B
#BSUB -N
###################################################################

# Parameters
GENOME_DIR="/scratch/projects/dark_genes/orthofinder_comps/Coral_ortho"  # Update this with the path to your genome/proteome files
OUTPUT_DIR="/scratch/projects/dark_genes/orthofinder_comps/Coral_ortho/results"  # Update this with the desired output directory
THREADS=8  # Adjust based on the number of available CPU cores

# Run OrthoFinder
echo "Running OrthoFinder..."
orthofinder -f "$GENOME_DIR" -o "$OUTPUT_DIR" -t $THREADS -a $THREADS

# Explanation of parameters:
# -f : Path to the folder containing genome/protein files
# -o : Path to the output folder
# -t : Number of threads for OrthoFinder to use
# -a : Number of threads for DIAMOND alignment

echo "OrthoFinder analysis complete. Results saved to $OUTPUT_DIR."
```

`bsub < of_coral_ortho.job`



## Cnid_ortho

```bash
mkdir Cnid_ortho
cd Cnid_ortho

ln -s /nethome/kxw755/genomes/GCA_032359415.1_Acer/acer_protein.faa
ln -s /nethome/kxw755/genomes/GCA_002571385.2_Spis/spis_protein.faa
ln -s /nethome/kxw755/genomes/GCA_948470475.1_Gfas/data/GCA_948470475.1/gfas_1.0.proteins.fasta
ln -s /nethome/kxw755/genomes/pdam_genome/pdam_proteins.fasta
ln -s /nethome/kxw755/genomes/Apoc_14110456/apoculata_proteins.fasta
ln -s /nethome/kxw755/genomes/GCF_021976095.1_Xenia/xenia_protein.faa
ln -s /nethome/kxw755/genomes/GCF_932526225.1_Nvec/nvec_protein.faa

nano of_cnid_ortho.job
```

```bash
#!/bin/bash

#BSUB -J OF_cnid_ortho
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o OF_cnid_ortho_%J.out
#BSUB -e OF_cnid_ortho_%J.err
#BSUB -B
#BSUB -N
###################################################################

# Parameters
GENOME_DIR="/scratch/projects/dark_genes/orthofinder_comps/Cnid_ortho"  # Update this with the path to your genome/proteome files
OUTPUT_DIR="/scratch/projects/dark_genes/orthofinder_comps/Cnid_ortho/results"  # Update this with the desired output directory
THREADS=16  # Adjust based on the number of available CPU cores

# Run OrthoFinder
echo "Running OrthoFinder..."
orthofinder -f "$GENOME_DIR" -o "$OUTPUT_DIR" -t $THREADS -a $THREADS

# Explanation of parameters:
# -f : Path to the folder containing genome/protein files
# -o : Path to the output folder
# -t : Number of threads for OrthoFinder to use
# -a : Number of threads for DIAMOND alignment

echo "OrthoFinder analysis complete. Results saved to $OUTPUT_DIR."
```

`bsub < of_cnid_ortho.job`


## Stem_ortho

```bash
mkdir Stem_ortho
cd Stem_ortho

ln -s /nethome/kxw755/genomes/GCA_032359415.1_Acer/acer_protein.faa
ln -s /nethome/kxw755/genomes/pdam_genome/pdam_proteins.fasta
ln -s /nethome/kxw755/genomes/GCF_932526225.1_Nvec/nvec_protein.faa

nano of_stem_ortho.job
```

```bash
#!/bin/bash

#BSUB -J OF_stem_ortho
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o OF_stem_ortho_%J.out
#BSUB -e OF_stem_ortho_%J.err
#BSUB -B
#BSUB -N
###################################################################

# Parameters
GENOME_DIR="/scratch/projects/dark_genes/orthofinder_comps/Stem_ortho"  # Update this with the path to your genome/proteome files
OUTPUT_DIR="/scratch/projects/dark_genes/orthofinder_comps/Stem_ortho/results"  # Update this with the desired output directory
THREADS=16  # Adjust based on the number of available CPU cores

# Run OrthoFinder
echo "Running OrthoFinder..."
orthofinder -f "$GENOME_DIR" -o "$OUTPUT_DIR" -t $THREADS -a $THREADS

# Explanation of parameters:
# -f : Path to the folder containing genome/protein files
# -o : Path to the output folder
# -t : Number of threads for OrthoFinder to use
# -a : Number of threads for DIAMOND alignment

echo "OrthoFinder analysis complete. Results saved to $OUTPUT_DIR."
```

`bsub < of_stem_ortho.job`


`scp -r kxw755@pegasus.ccs.miami.edu:/scratch/projects/dark_genes/orthofinder_comps/Pdam_Spis/results/Results_Dec24/Orthogroups/Orthogroups.tsv ./Pdam_Spis_Orthogroups.tsv`






## Gfas_Pdam

```bash
mkdir Gfas_Pdam
cd Gfas_Pdam

ln -s /nethome/kxw755/genomes/GCA_948470475.1_Gfas/data/GCA_948470475.1/gfas_1.0.proteins.fasta
ln -s /nethome/kxw755/genomes/pdam_genome/pdam_proteins.fasta

nano of_Gfas_Pdam.job
```

```bash
#!/bin/bash

#BSUB -J OF_Gfas_Pdam
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o OF_Gfas_Pdam_%J.out
#BSUB -e OF_Gfas_Pdam_%J.err
#BSUB -B
#BSUB -N
###################################################################

# Parameters
GENOME_DIR="/scratch/projects/dark_genes/orthofinder_comps/Gfas_Pdam"  # Update this with the path to your genome/proteome files
OUTPUT_DIR="/scratch/projects/dark_genes/orthofinder_comps/Gfas_Pdam/results"  # Update this with the desired output directory
THREADS=16  # Adjust based on the number of available CPU cores

# Run OrthoFinder
echo "Running OrthoFinder..."
orthofinder -f "$GENOME_DIR" -o "$OUTPUT_DIR" -t $THREADS -a $THREADS

# Explanation of parameters:
# -f : Path to the folder containing genome/protein files
# -o : Path to the output folder
# -t : Number of threads for OrthoFinder to use
# -a : Number of threads for DIAMOND alignment

echo "OrthoFinder analysis complete. Results saved to $OUTPUT_DIR."
```

`bsub < of_Gfas_Pdam.job`


`scp -r kxw755@pegasus.ccs.miami.edu:/scratch/projects/dark_genes/orthofinder_comps/Gfas_Pdam/results/Results_Dec27/Orthogroups/Orthogroups.tsv ./Gfas_Pdam_Orthogroups.tsv`



## Nvec_Nvec

```bash
mkdir Nvec_Nvec2
cd Nvec_Nvec2

ln -s /nethome/kxw755/genomes/GCF_932526225.1_Nvec/nvec_protein.faa
ln -s /nethome/kxw755/genomes/Nvec_200/NV2g.20240221.protein.fa

nano of_nvec_nvec.job
```

```bash
#!/bin/bash

#BSUB -J OF_nvec_nvec
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o OF_nvec_nvec_%J.out
#BSUB -e OF_nvec_nvec_%J.err
#BSUB -B
#BSUB -N
###################################################################

# Parameters
GENOME_DIR="/scratch/projects/dark_genes/orthofinder_comps/Nvec_Nvec2"  # Update this with the path to your genome/proteome files
OUTPUT_DIR="/scratch/projects/dark_genes/orthofinder_comps/Nvec_Nvec2/results"  # Update this with the desired output directory
THREADS=16  # Adjust based on the number of available CPU cores

# Run OrthoFinder
echo "Running OrthoFinder..."
orthofinder -f "$GENOME_DIR" -o "$OUTPUT_DIR" -t $THREADS -a $THREADS

# Explanation of parameters:
# -f : Path to the folder containing genome/protein files
# -o : Path to the output folder
# -t : Number of threads for OrthoFinder to use
# -a : Number of threads for DIAMOND alignment

echo "OrthoFinder analysis complete. Results saved to $OUTPUT_DIR."
```

`bsub < of_nvec_nvec.job`

`scp -r kxw755@pegasus.ccs.miami.edu:/scratch/projects/dark_genes/orthofinder_comps/Nvec_Nvec2/results/Results_Dec27/Orthogroups/Orthogroups.tsv ./Nvec_Nvec2_Orthogroups.tsv`
