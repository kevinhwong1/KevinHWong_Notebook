---
layout: post
title: Zymo-DNA-RNA-Extract-P.astreoides-Genome
date: '2019-03-13'
categories: Processing, Protocols
tags: P.asteroides
---

# RNA/DNA extractions for P. astreoides genome samples

Kit: ZYMO Quick-DNA/RNA Miniprep Plus Kit

## Goal
* To extract DNA and RNA for an Illumina RNAseq library prep workshop

## Experimental design
Table 1: Sample descriptions for this extraction

| Sample # | Sample Type | Coral ID | Size of Sample |     Treatment   |
|:--------:|:-----------:|:--------:|:--------------:|:----------------:
|     1    |   Fragment  |   PG #1  | ~10mm diameter |  23&deg;C + Air |
|     2    |   Fragment  |   PG #2  | ~10mm diameter | -22&deg;C + Dark|

* I sampled duplicates of each sample, preserved them in 500&mu;L DNA/RNA shield in ZR BBashing lysis tubes with 0.1mm and 0.5mm beads and stored at -80&deg;C

#### Overall workflow
1. Prepare samples
2. Extract and prep DNA with column (elute RNA and Protein)
3. Extract and prep RNA with column (elute Protein)
4. Quantify DNA, RNA, and protein extracts

## Protocol preparation

* Using slightly modified [Zymo Duet DNA/RNA extraction protocol](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/Zymo_quick-dna-rna_miniprep_plus_kit.pdf) which will extract both DNA and RNA at the same time (Below are summary steps)

#### Reagents and supplies
* RNase-free Water
* 100% ethanol, ACS grade or better
* 10mM Tris HCl pH 8.0 made with RNase-free water
* ZR BBashing lysis tubes with 0.1mm and 0.5mm beads

#### Equipment
* Rocking oven that can be set to 70°C
* RNase away and a designated RNase free space
* Tabletop and larger centrifuges for 1.5mL and 50mL tubes capable of 12,000 x g

#### Clipper Sterilization
* Rub down clippers with:
  1. 10% Bleach solution
  2. DI water
  3. 90% ethanol
  4. RNAse free water

#### Notes before starting
* Wipe down benchtop with RNase away and have the spray bottle and kimwipes on-hand to use frequently

## Sample preparation

#### Fragment preparation
* Sterilize clippers (as outlined above)
* Clip off ~10mm in diameter of tissue and skeleton and place into ZR BBashing lysis tube
* Add 5000 &mu;L of DNA/RNA shield (may be too much)
* Vortex for ~1 minute
* Remove supernatant and transfer into 1.5mL tube
  * I was able to remove ~450&mu;L

## DNA Extraction
1. Add equal volumes of DNA lysis buffer as sample volume to the 5 mL tube with sample
  * 450 &mu;L for the fragment samples
  * Vortex to mix
2. Transfer 700 &mu;L into DNA filter column (yellow)
  * Centrifuge at 16,000 rcf for 30 seconds
  * Remove flow through liquid and transfer into a new 5 mL tube labeled for RNA
3. Repeat steps 2 until all liquid is gone
4. Add 400 &mu;L of DNA/RNA prep buffer to column
  * Centrifuge at 16,000 rcf for 30 seconds
  * Remove the flow through and transfer to waste
5. Add 700 &mu;L of DNA/RNA wash buffer
  * Centrifuge at 16,000 rcf for 30 seconds
  * Remove the flow through and transfer to waste
6. Add 400 &mu;L of DNA/RNA wash buffer
  * Centrifuge at 16,000 rcf for 2 minutes
  * Remove the flow through and transfer to waste
7. Removed column and place into a new sterile 1.5 mL tube
8. Add 50 &mu;L of 10 mM Tris HCL (warmed to 55&deg;C) directly to filter
  * Incubate at room temperature for 5 minutes
  * Centrifuge at 16,000 rcf for 30 seconds
  * Keep flow through in tube
9. Add another 50 &mu;L of 10 mM Tris HCL (warmed to 55&deg;C) directly to filter
  * Incubate at room temperature for 5 minutes
  * Centrifuge at 16,000 rcf for 30 seconds
  * Discard column and keep flow through in tube
10. After quantifying on the Qubit, I aliquoted into two 1.5 eppindorf tubes (10 &mu;L and the remaining 89 &mu;L) and stored at -20&deg;C
  * Tubes were labelled: *DNA #1 KW 03/13* and *DNA #2 KW 03/13*

## RNA Extraction
11. Add equal volumes of 100% Ethanol as sample volume to the 5 mL tube with sample
  * 900 &mu;L for the fragment samples
  * Vortex to mix
12. Transfer 700 &mu;L into RNA filter column (green)
  * Centrifuge at 16,000 rcf for 30 seconds
  * Remove flow through liquid and transfer into a new 5 mL tube labeled for Protein and store in -80&deg;C
13. Repeat steps 12 until all liquid is gone
14. Add 400 &mu;L of DNA/RNA wash buffer to column
  * Centrifuge at 16,000 rcf for 30 seconds
  * Remove the flow through and transfer to waste
15. Make DNase I master mix: [75μL DNA digestion buffer and 5μL DNase] X n (sample #)
  * Two samples, 150μL DNA digestion buffer and 10μL DNase
16. Add 80μL of DNase I master mix directly to the filter of each column. Incubate at room temperature for ~15 minutes
17. Add 400 &mu;L of DNA/RNA prep buffer to column
  * Centrifuge at 16,000 rcf for 30 seconds
  * Remove the flow through and transfer to waste
18. Add 700 &mu;L of DNA/RNA wash buffer
  * Centrifuge at 16,000 rcf for 30 seconds
  * Remove the flow through and transfer to waste
19. Add 400 &mu;L of DNA/RNA wash buffer
  * Centrifuge at 16,000 rcf for 2 minutes
  * Remove the flow through and transfer to waste
20. Removed column and place into a new sterile 1.5 mL tube
21. Add 50 &mu;L of RNAase-free water (warmed to 55&deg;C) directly to filter
  * Incubate at room temperature for 5 minutes
  * Centrifuge at 16,000 rcf for 30 seconds
  * Keep flow through in tube
22. Add another 50 &mu;L of RNAase-free water (warmed to 55&deg;C) directly to filter
  * Incubate at room temperature for 5 minutes
  * Centrifuge at 16,000 rcf for 30 seconds
  * Discard column and keep flow through in tube
23. Aliquot into appropriate tubes for storage at -80&deg;C

## Quantification

### Qubit

* For the Qubit readings, 10 &mu;L of each standard were used and 1 &mu;L of each sample were used.

Table 2: Qubit readings for the extracted DNA and RNA samples.

|     | Standard 1 | Standard 2 | #1 (ng/uL) | #2 (ng/uL) |
|:---:|:----------:|:----------:|:----------:|:----------:|
| DNA |   199.37   |  21533.38  |    10.8    |    29.2    |
| RNA |   401.46   |  11418.56  |    23.6    |    34.2    |

### Tapestation
* Only RNA samples were run on the tape station

#### Tapestation results
Sample #1:
![Sample#1]({{ site.baseurl}}/images/20190313_Sample1.png "Sample#1")

Sample #2:
![Sample#2]({{ site.baseurl}}/images/20190313_Sample2.png "Sample#2")

## Conclusions
* Samples are good for the Illumina RNA-Seq workshop
* Changes in Ethanol addition in RNA extraction step potentially produced higher concentrations of RNA 
