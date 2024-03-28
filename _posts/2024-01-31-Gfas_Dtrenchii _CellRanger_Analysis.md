---
layout: post
title: Gfas_Syms_CellRanger_Analysis
date: '2024-01-31'
categories: Analysis
tags: Cnidarin immunity lab
---

# Goal: To re-analyze the 10X single cell data with the combined genomic references of Galaxea fasicularis and potential symbionts 

**Reference Genomes:**
- Galaxea fasicularis v1 
    - From Reef Genomics
    - [Download site](http://gfas.reefgenomics.org/)
- Durisdinum trenchii (SCF082)
    - [Dougan et al. 2022](https://www.biorxiv.org/content/10.1101/2022.04.10.487810v1.full.pdf)
    - [Download site](https://espace.library.uq.edu.au/view/UQ:27da3e7)
- Breviolum (Symbiodinium minutum)
    - [Download site](https://marinegenomics.oist.jp/symb/viewer/download?project_id=21)
- Cladocopium 
    - [Download site](https://marinegenomics.oist.jp/symb/viewer/download?project_id=40)
    - alternate [Download site](https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=Clago1)
        - Reference: Liu H, Stephens TG, González-Pech RA, Beltran VH, Lapeyre B, Bongaerts P, Cooke I, Aranda M, Bourne DG, Forêt S, Miller DJ, van Oppen MJH, Voolstra CR, Ragan MA, Chan CX Symbiodinium genomes reveal adaptive evolution of functions related to coral-dinoflagellate symbiosis. Commun Biol. 2018;1():95. doi: 10.1038/s42003-018-0098-3
- Symbiodinium microadriaticum
    - [Download site](https://marinegenomics.oist.jp/symb/viewer/download?project_id=37)
    - alternate [Download Site](https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=Symmic1)
        - Reference: Aranda M, Li Y, Liew YJ, Baumgarten S, Simakov O, Wilson MC, Piel J, Ashoor H, Bougouffa S, Bajic VB, Ryu T, Ravasi T, Bayer T, Micklem G, Kim H, Bhak J, LaJeunesse TC, Voolstra CR Genomes of coral dinoflagellate symbionts highlight evolutionary adaptations conducive to a symbiotic lifestyle. Sci Rep. 2016 Dec 22;6():39734. doi: 10.1038/srep39734

 
# File upload and organization

```bash
mkdir sym_genomes
cd sym_genomes/

scp -r symA3_symb.gff.gz kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/sym_genomes/symA3_symb.gff.gz
scp -r symC_40.gff.gz kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/sym_genomes/symC_40.gff.gz
scp -r symbB.v1.2.augustus.gff3.gz kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/sym_genomes/symbB.v1.2.augustus.gff3.gz

scp -r symA3_37.fasta.gz kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/sym_genomes/symA3_37.fasta.gz
scp -r symbB.v1.0.genome.fa.gz kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/sym_genomes/symbB.v1.0.genome.fa.gz
scp -r symC_scaffold_40.fasta.gz kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/sym_genomes/symC_scaffold_40.fasta.gz
```

the sym C and A gff files did not have exon reads and cell ranger could not make the reference genome. I will download the following to see if these work

```bash

scp -r Clago1_AssemblyScaffolds_Repeatmasked.fasta.gz kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/sym_genomes/Clago1_AssemblyScaffolds_Repeatmasked.fasta.gz
scp -r Clago1_GeneCatalog_genes_20200812.gff3.gz kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/sym_genomes/Clago1_GeneCatalog_genes_20200812.gff3.gz

scp -r Symmic1_AssemblyScaffolds_Repeatmasked.fasta.gz kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/sym_genomes/Symmic1_AssemblyScaffolds_Repeatmasked.fasta.gz
scp -r Symmic1_all_genes_20180603.gff3.tgz kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/sym_genomes/Symmic1_all_genes_20180603.gff3.tgz

```

# Unzip and convert to gff3

```bash
gunzip symA3_symb.gff.gz
mv symA3_symb.gff symA3_symb.gff3

gunzip symbB.v1.2.augustus.gff3.gz

gunzip symC_40.gff.gz
mv symC_40.gff symC_40.gff3

gunzip Dtrenchii_SCF082_ANNOT_gff.gz
mv Dtrenchii_SCF082_ANNOT_gff Dtrenchii_SCF082_ANNOT.gff3

gunzip symA3_37.fasta.gz
gunzip symbB.v1.0.genome.fa.gz
gunzip symC_scaffold_40.fasta.gz
gunzip Dtrenchii_SCF082_ASSEMBLY_fasta.gz

mv Dtrenchii_SCF082_ASSEMBLY_fasta Dtrenchii_SCF082_ASSEMBLY.fasta
```

```bash
gunzip Clago1_AssemblyScaffolds_Repeatmasked.fasta.gz
gunzip Clago1_GeneCatalog_genes_20200812.gff3.gz

tar -xvzf Symmic1_all_genes_20180603.gff3.tgz
cd /nethome/kxw755/sym_genomes/global/projectb/sandbox/fungal/data/Symbiodinium_microadriaticum/Symmic1/download
cp Symmic1.ExternalModels.gff3 ../../../../../../../../

gunzip Symmic1_AssemblyScaffolds_Repeatmasked.fasta.gz
```

these gff3 files seem to have all the genomic information rather than just CDS

# Converting this GFF3 file to GTF:

`nano convert.job`

```bash
#!/bin/bash

#BSUB -J scSeq_convertGTF
#BSUB -q general
#BSUB -P dark_genes
#BSUB -n 6
#BSUB -W 12:00
#BSUB -u kxw755@earth.miami.edu
#BSUB -o dtre_convertGTF.out
#BSUB -e dtre_convertGTF.err

###################################################################

module load cufflinks/2.2.1 

gffread Dtrenchii_SCF082_ANNOT.gff3 -T -o Dtrenchii_SCF082_ANNOT.gtf
```

`bsub < convert.job`


`nano convert2.job`

```bash
#!/bin/bash

#BSUB -J scSeq_convertGTF
#BSUB -q general
#BSUB -P dark_genes
#BSUB -n 6
#BSUB -W 12:00
#BSUB -u kxw755@earth.miami.edu
#BSUB -o syms_convertGTF.out
#BSUB -e syms_convertGTF.err

###################################################################

module load cufflinks/2.2.1 

gffread symA3_symb.gff3 -T -o symA3_symb.gtf
gffread symbB.v1.2.augustus.gff3 -T -o symbB.v1.2.augustus.gtf
gffread symC_40.gff3 -T -o symC_40.gtf
```

`bsub < convert2.job`


`nano convert3.job`

```bash
#!/bin/bash

#BSUB -J scSeq_convertGTF
#BSUB -q general
#BSUB -P dark_genes
#BSUB -n 6
#BSUB -W 12:00
#BSUB -u kxw755@earth.miami.edu
#BSUB -o syms_convertGTF3.out
#BSUB -e syms_convertGTF3.err

###################################################################

module load cufflinks/2.2.1 

gffread Clago1_GeneCatalog_genes_20200812.gff3 -T -o Clago1_GeneCatalog_genes_20200812.gtf
gffread Symmic1.ExternalModels.gff3 -T -o Symmic1.ExternalModels.gtf
```

`bsub < convert3.job`


# Making combined genome reference

`nano mkgtf_syms.job`

```bash
#!/bin/bash

#BSUB -J scSeq_mkgtf_syms
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o gfas_syms_mkgtf_v2.out
#BSUB -e gfas_syms_mkgtf_v2.err

###################################################################

cellranger mkref \
--genome=gfas_1.0 --fasta=/nethome/kxw755/Gfas_v1/gfas_final_1.0.fasta --genes=/nethome/kxw755/Gfas_v1/gfas_1.0.filtered.gtf \
--genome=symA --fasta=/nethome/kxw755/sym_genomes/Symmic1_AssemblyScaffolds_Repeatmasked.fasta --genes=/nethome/kxw755/sym_genomes/Symmic1.ExternalModels.gtf \
--genome=symB --fasta=/nethome/kxw755/sym_genomes/symbB.v1.0.genome.fa --genes=/nethome/kxw755/sym_genomes/symbB.v1.2.augustus.gtf \
--genome=symC --fasta=/nethome/kxw755/sym_genomes/Clago1_AssemblyScaffolds_Repeatmasked.fasta --genes=/nethome/kxw755/sym_genomes/Clago1_GeneCatalog_genes_20200812.gtf \
--genome=symD --fasta=/nethome/kxw755/sym_genomes/Dtrenchii_SCF082_ASSEMBLY.fasta --genes=/nethome/kxw755/sym_genomes/Dtrenchii_SCF082_ANNOT.gtf \
--memgb 20
```

`bsub < mkgtf_syms.job`

error: 

```bash
Fatal LIMIT error: the number of junctions to be inserted on the fly =2833619 is larger than the limitSjdbInsertNsj=1000000
Fatal LIMIT error: the number of junctions to be inserted on the fly =2833619 is larger than the limitSjdbInsertNsj=1000000
SOLUTION: re-run with at least --limitSjdbInsertNsj 2833619

Mar 26 16:12:54 ...... FATAL ERROR, exiting

```

Maybe lets try this without symA and see if this frees up the memory/RAM usage

`nano mkgtf_syms2.job`

```bash
#!/bin/bash

#BSUB -J scSeq_mkgtf_syms2
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o gfas_syms_mkgtf_v2.out
#BSUB -e gfas_syms_mkgtf_v2.err

###################################################################

cellranger mkref \
--genome=gfas_1.0 --fasta=/nethome/kxw755/Gfas_v1/gfas_final_1.0.fasta --genes=/nethome/kxw755/Gfas_v1/gfas_1.0.filtered.gtf \
--genome=symB --fasta=/nethome/kxw755/sym_genomes/symbB.v1.0.genome.fa --genes=/nethome/kxw755/sym_genomes/symbB.v1.2.augustus.gtf \
--genome=symC --fasta=/nethome/kxw755/sym_genomes/Clago1_AssemblyScaffolds_Repeatmasked.fasta --genes=/nethome/kxw755/sym_genomes/Clago1_GeneCatalog_genes_20200812.gtf \
--genome=symD --fasta=/nethome/kxw755/sym_genomes/Dtrenchii_SCF082_ASSEMBLY.fasta --genes=/nethome/kxw755/sym_genomes/Dtrenchii_SCF082_ANNOT.gtf 
```

`bsub < mkgtf_syms2.job`

```bash
Fatal LIMIT error: the number of junctions to be inserted on the fly =1992803 is larger than the limitSjdbInsertNsj=1000000
Fatal LIMIT error: the number of junctions to be inserted on the fly =1992803 is larger than the limitSjdbInsertNsj=1000000
SOLUTION: re-run with at least --limitSjdbInsertNsj 1992803

Mar 27 11:04:18 ...... FATAL ERROR, exiting
```

okay - according to this [link](https://kb.10xgenomics.com/hc/en-us/articles/360003877352-How-can-I-modify-the-STAR-alignment-parameters-in-Cell-Ranger) I need to adjust a parameter in the STAR file, but I don't think I have access. I am going to just map the symbionts with out the Galaxea genome and see what hits. 


`nano mkgtf_syms_only.job`

```bash
#!/bin/bash

#BSUB -J scSeq_mkgtf_symsonly
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o syms_mkgtf_only.out
#BSUB -e syms_mkgtf_only.err

###################################################################

cellranger mkref \
--genome=symA --fasta=/nethome/kxw755/sym_genomes/Symmic1_AssemblyScaffolds_Repeatmasked.fasta --genes=/nethome/kxw755/sym_genomes/Symmic1.ExternalModels.gtf \
--genome=symB --fasta=/nethome/kxw755/sym_genomes/symbB.v1.0.genome.fa --genes=/nethome/kxw755/sym_genomes/symbB.v1.2.augustus.gtf \
--genome=symC --fasta=/nethome/kxw755/sym_genomes/Clago1_AssemblyScaffolds_Repeatmasked.fasta --genes=/nethome/kxw755/sym_genomes/Clago1_GeneCatalog_genes_20200812.gtf \
--genome=symD --fasta=/nethome/kxw755/sym_genomes/Dtrenchii_SCF082_ASSEMBLY.fasta --genes=/nethome/kxw755/sym_genomes/Dtrenchii_SCF082_ANNOT.gtf 
```

`bsub < mkgtf_syms_only.job`



# Rerunning Count with the combined genome

## Run 1

cd /nethome/kxw755/20231004_SingleCell_DG/

`nano count_W045_deep_syms.job`

```bash
#BSUB -J count_W045_deep_syms
#BSUB -q general
#BSUB -P dark_genes
#BSUB -n 6
#BSUB -W 120:00
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count.out
#BSUB -e count.err
#BSUB -B
#BSUB -N
###################################################################

cellranger count \
 --id=W-045_1_deep_syms \
 --transcriptome=/nethome/kxw755/sym_genomesgfas_1.0_and_symA3_and_symB_and_symC_and_symD \
 --fastqs=/nethome/kxw755/20231004_SingleCell_DG \
 --sample=AndradeRodriguez-15275-001_GEX3
```

`bsub < count_W045_deep_syms.job`

Export:

```bash
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20231004_SingleCell_DG/W-045_1_deep_combgenome/outs/filtered_feature_bc_matrix.h5 /Users/kevinwong/MyProjects/DarkGenes_Bleaching_Comparison/output/CellRanger/20240326_W045_combgeno_filt_feature_bc_matrix.h5

scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20231004_SingleCell_DG/W-045_1_deep_combgenome/outs/web_summary.html /Users/kevinwong/MyProjects/DarkGenes_Bleaching_Comparison/output/CellRanger/20240326_W045_combgeno_web_summary.html
```
