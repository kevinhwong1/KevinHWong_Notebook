---
layout: post
title: Citrate Synthase Troubleshooting
date: '2019-03-14'
categories: Protocols, Processing
tags: citrate synthase
---

Troubleshooting Citrate Synthase Protocol

## Goal
* To create a protocol to quantify citrate synthase activities in *Astrangia poculata*

## Samples
* Sampled extended polyps with scissors and foreceps
* Sampled one symbiotic (#1) and aposymbiotic (#2) coral

## Overall workflow
1. Collect/thaw samples
2. Prepare and lyse samples
3. Quantify protein
4. Prepare samples
5. Add reagents to samples and +/- controls
6. Read on spectrophotometer (baseline)
7. Add Oxaloacetic Acid and read on spectrophotometer
8. Calculate

## Preparation

This protocol is modified from [Sigma-Aldrich Citrate Synthase Protocol](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/SA_Citrate_Synthase_Assay.pdf), [Hawkins et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4834350/), and [Spinazzi et al. 2012](https://www.ncbi.nlm.nih.gov/pubmed/22653162).

#### Kits
* [Sigma-Aldrich Citrate Synthase Kit](https://www.sigmaaldrich.com/catalog/product/sigma/cs0720?lang=en&region=US)
* [Pierce BCA Protein Assay Kit](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/Pierce_BCA_Protein_Assay.pdf)

#### Equipment
* Vortex
* 96 well plate (flat bottom)
* Plate spectrophotometer
* Tabletop and larger centrifuges for 1.5mL tubes

#### Reagents/Supplies not included in kits

* DI water
* 100% ethanol
* Lysis Buffer
* ZR BBashing lysis tubes with 0.1mm and 0.5mm beads

#### Reagent Preparation

* **Lysis buffer (25 mM Tris (pH 7.8), 1 mM EDTA, 10% gylcerol [v/v])**
  * To make 100 mL of lysis buffer, add:
    * 2.5 mL 1 M Tris (pH 8)
    * 0.2 mL 0.5 M EDTA
    * 10 mL 100% glycerol [v/v]
    * 87.3 mL DI water
  * Store at room temperature
* **Assay buffer (1X)**
  * Solution comes as 5X in the kit, we need to dilute to 1X
  * Today I needed ~3.5 mL of assay buffer (for reaction, Acetyl CoA solution, and citrate synthase dilution)
    * 0.7 mL of 5X Assay buffer + 2.8 mL DI water = 3.5 mL 1X Assay buffer
* **Acetyl CoA solution (30 mM)**
  * Dissolve entire vial with 1 mL DI water
  * Mix until homogenous
  * Aliquot into 5 vials of 200 &mu;L and store at -20&deg;C
* **Oxaloacetate (OAA) solution (10 mM)**
  * Dissolve 1.3 mg of OAA powder in 1 mL of 1X Assay Buffer
  * Store at -20&deg;C for maximum 1 week
* **DNTB solution (10 mM)**
  * Dissolve entire vial with 1 mL 100% Ethanol
  * Mix until homogenous
  * Aliquot into 5 vials of 200 &mu;L and store at -20&deg;C
* **Citrate Synthase**
  * Diluted 2 &mu;L of Citrate Synthase solution with 60 &mu;L of 1X Assay buffer

## Protocol

**1. Collect/Thaw samples**
  * Collected one polyp from a dark (symbiotic) and light (aposymbiotic) *A. poculata* colonies with scissors and foreceps
  * Placed immediately into Bead tube with 500 &mu;L of lysis Buffer

**2. Prepare and lyse samples**
  * Vortex samples for ~2 minutes
  * Remove liquid and transfer into a new, labelled 1.5 mL tube
  * Centrifuge for 1.5 minutes at 16000 rcf
  * Transfer supernatant to a new, labelled 1.5 mL tube, discard pellet
  * Aliquot 80 &mu;L for protein quantification and freeze the rest at -80 &deg;C

**3. Quantify protein**
  * I used the protocol from the [Pierce BCA Protein Assay Kit](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/Pierce_BCA_Protein_Assay.pdf)
  * My analysis for protein content (mg/mL) per sample can be found [here](https://github.com/kevinhwong1/methods_testing/blob/master/RAnalysis/Scripts/TProtein_20190314.R)
  * Protein content in each sample
    * Sample 1: 0.774 mg/mL
    * Sample 2: 0.275 mg/mL

**4. Prepare samples**
  * According to Hawkins et al. 2016, we need to dilute to a yield of 2-6 &mu;g of protein in 20 &mu;L of sample per well
    * Sample 1:
      * 0.774 mg/mL = 0.774 &mu;g/&mu;L
      * 0.774 &mu;g/&mu;L x 20 &mu;L = 15.4 &mu;g of Protein
      * Need to dilute by 3X
        * 6.6 &mu;L of Sample 1 + 13.5 DI water = 20 &mu;L with ~5.1 &mu;g of protein
    * Sample 2:
      * 0.275 mg/mL = 0.275 &mu;g/&mu;L
      * 0.275 &mu;g/&mu;L x 20 &mu;L = 5.2 &mu;g of Protein
      * No need to dilute
  * The above calculations are for 1 well, we need to multiply everything by 3 to have 3 replicate wells per sample
    * In 1.5 mL eppindorf tubes:
      * Sample 1: 19.8 &mu;L Sample 1 + 40.2 DI water
      * Sample 2: 60 &mu;L Sample 2
    * Store on ice
  * For the negative control, add 60 &mu;L of lysis buffer to a new, labeled 1.5 mL tube
  * For the positive control, use the diluted citrate synthase tube from the reagent preparation step

**5. Add reagents to samples and +/- controls**
  * Prepare a master mix of the reagents prior to adding it to the sample tubes
    *  In each well, we need 170 &mu;L of master mix with the following concentrations:
      * 3.3 &mu;L of 30 mM Acetyl CoA solution (making it 588 &mu;M)
      * 5 &mu;L of 10 mM DNTB solution (making it 294 &mu;M)
      * 161.7 &mu;L of 1X Assay buffer
    * We have 2 samples, 1 negative control (lysis buffer) and 1 positive control (citrate synthase), with three replicates for each
      * Therefore we need enough master mix for 12 wells
      * I made enough for 15 wells for a total of 2550 &mu;L of master mix
        * 49.5 &mu;L of 30 mM Acetyl CoA solution (making it 588 &mu;M)
        * 75 &mu;L of 10 mM DNTB solution (making it 294 &mu;M)
        * 2425.5 &mu;L of 1X Assay buffer
  * Mix well by vortexing
  * Add 510 &mu;L of master mix (170 &mu;L x 3 replicates) to the sample and control tubes
  * Mix well by vortexing
  * Add 190 &mu;L of sample into the corresponding wells on the 96 well plate

**6. Read on spectrophotometer (baseline activity)**
  * On kinetic mode, I read these samples at 412 nm for 3 mins with a 49 second interval, resulting in 4 data points

**7. Add Oxaloacetic Acid and read on spectrophotometer**
  * Add 10 &mu;L of OAA solution to each well
  * Use a multi-channel pipette to start all reactions simultaneously
  * Mix well by pipetting
  * Re-read the plate on the same settings as Step 6

**8. Calculate**
* I calculated the citrate synthase activity using RStudio and my script can be found [here](https://github.com/kevinhwong1/methods_testing/blob/master/RAnalysis/Scripts/Citrate_Synthase_20190315.R)

## Conclusions
* This protocol worked for *A. poculata*, with citrate synthase activities of:
  * Sample 1 (symbiotic) = **68.86 mU/mg**
  * Sample 2 (aposymbiotic) = **101.51 mU/mg**
* These values fall in the same range of citrate synthase activities in symbiotic *Exaiptasia pallida* sea anemones (Hawkins et al. 2016)
* However, our results differ from Hawkins et al. 2016, with our aposymbiotic colonies resulting in higher citrate synthase activities than the symbiotic colonies
  * After this preliminary test, we will test again with larger samples sizes to see if this trend holds true
