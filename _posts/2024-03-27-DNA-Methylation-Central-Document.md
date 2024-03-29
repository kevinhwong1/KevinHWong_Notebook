---
layout: post
title: DNA Methylation Central Document
date: '2024-03-27'
categories: Protocols
tags: Central Document
---

# DNA Methylation Sequencing Analysis Central Document

This will be my working document to connect all the datasets and different analysis pipelines together. 

## Projects

* [*Porties astreoides* Thermal Transplant Molecular](https://github.com/kevinhwong1/Thermal_Transplant_Molecular)
    * How does pervious bleaching events impact adult and larval coral DNA methyation patters?
    * How does parental thermal history influence coral larval DNA methylation patterns and phenotypes?
* [*Crassostrea virginica* nutition multiomic integration](https://github.com/hputnam/Cvir_Nut_Int)
    * How does nutirition enrichment alter the gut microbiome, gene expression, and DNA methylation patterns of the eastern oyster *Crassostrea virginica*?

## While Genome Bisulfite library prepration

Whole genome bisulfite sequencing (WGBS) is an approach to obtain single base-pair resolution of DNA methylation within your sample. I have developed a library prepration protocol with the Zymo Pico Methyl Seq kit for the coral *Porites astreoides* linked [here](https://kevinhwong1.github.io/KevinHWong_Notebook/Thermal-Transplant-WGBS-PicoMethyl-Protocol/).

## Bioinformatic pipelines

* DNA methylation QC and Quantification
    * This workflow describes how to use [NF-core Methylseq](https://nf-co.re/methylseq/2.6.0) to quanitfy methylation from raw reads, filter for 5X or 10X coverage, and prepare data for further analyses in R. 
    * The HPC workflow is linked [here](https://github.com/kevinhwong1/Thermal_Transplant_Molecular/blob/main/scripts/Past_WGBS_Workflow.md)
* DNA Methylation Characterization
    * This workflow describes how to characterize DNA methylation in various genomic features
    * The HPC workflow is linked [here](https://github.com/kevinhwong1/Thermal_Transplant_Molecular/blob/main/scripts/Past_Genomic_Feature_Analysis_20221014.md)
    * The R workflow is linked [here](https://github.com/kevinhwong1/Thermal_Transplant_Molecular/blob/main/scripts/Methylation_Characterization.md)
* SNP extraction from Whole Genome Bisulfite Sequencing (WGBS) data
    *  [BS-SNPer](https://github.com/hellbelly/BS-Snper) workflow is linked [here](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/_posts/2022-11-21-Testing-BS-SNPer-on-WGBS-data.md)
    * [Bis-SNP](https://people.csail.mit.edu/dnaase/bissnp2011/) workflow is linked [here](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/_posts/2023-08-04-BisSNP-Analysis.md)
    * [EpiDiverse](https://github.com/EpiDiverse/snp/blob/master/docs/usage.md) workflow is linked [here](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/_posts/2023-03-15-EpiDiverse-SNP-Analysis.md)
* Plotting, Statsitics, WGCNA, Functional enrichment
    * This R script is linked [here](https://github.com/kevinhwong1/Thermal_Transplant_Molecular/blob/main/scripts/WGBS_GM.Rmd)