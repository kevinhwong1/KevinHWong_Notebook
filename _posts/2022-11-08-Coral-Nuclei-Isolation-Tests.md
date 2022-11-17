---
layout: post
title: Coral Nuclei Isolation Tests
date: '2022-11-08'
categories: Protocols
tags: Coral, nuclei
---

This post is to troubleshoot coral nuclei isolations for single cell ATAC-seq. I will be troubleshooting protocols from multiple papers and on a few species.

## Trial 1

This trial is based directly from [Roquis et al. 2022](https://wellcomeopenresearch.org/articles/6-195/v2).

I am also going to think about preserving the nuclei using [this protocol](http://collaslab.org/wp-content/uploads/2011/07/Isolation_somatic_cell_nuclei.pdf)

### Materials

Equipment

* Cold centrifuge with 15mL rotor
* Air Brush
* 15mL tubes
* Sterile whirlpak
* Dounce homogenizer

Chemicals
* KCl
* NaCl
* MgCl<sub>2</sub>
* Tris/Cl pH 7.4
* Sucrose
* Sodium buryrate
* Dithiotreitol (DTT)
* Phenylmethanesulfonyl fluoride (PMSF)
* cOmplete protease inhibitor (Roche cat #11836145001)
* NP-40
* CaCl<sub>2</sub>
* RNAse inhibitor 


### Methods

#### Make buffers

Table 1: Chemical concentrations of each buffer. I will break down the volumes in the following tables but this is from the paper.

|                                                          | **Buffer 1** | **Buffer 2** | **Buffer 3** | **MNase Buffer** |
|:--------------------------------------------------------:|:------------:|:------------:|:------------:|:----------------:|
|                          **KCl**                         |     30 mM    |     30 mM    |     30 mM    |                  |
|                         **NaCl**                         |    500 mM    |    500 mM    |    500 mM    |      500 mM      |
|                         **MgCl2**                        |    2.5 mM    |    2.5 mM    |    2.5 mM    |       40 mM      |
|                    **Tris/Cl pH 7.4**                    |     15 mM    |     15 mM    |     15 mM    |       10 mM      |
|                        **Sucrose**                       |      1 M     |      1 M     |     1.5 M    |        1 M       |
|                    **Sodium buryrate**                   |     5 mM     |     5 mM     |     5 mM     |       5 mM       |
|                  **Dithiotreitol (DTT)**                 |     25 mM    |     25 mM    |     25 mM    |                  |
|         **Phenylmethanesulfonyl fluoride (PMSF)**        |    0.5 mM    |    0.5 mM    |    0.5 mM    |      0.1 mM      |
| **cOmplete protease inhibitor** |      1 X     |      1 X     |      1 X     |                  |
|                         **NP-40**                        |              |     0.3 %    |              |                  |
|                         **CaCl2**                        |              |              |              |       10 mM      |



#### Tissue dissociation (Adapted from Grace Snyder)

* Mechanical Dissociation: 
    * Remove cells/tissue from skeleton with air brush and cold cell media
    * “Cell media” refers to 3.3X PBS (salinity ~ 35 ppt; Mg- and Ca-free), 2% FBS, 20 mM Hepes, pH 7.4
    * Mechanical dissociation is preferred due to shorter preparation time, but enzymatic dissociations are also optional; best for coral are dispase and liberase (~45 minutes for 1 cm2  fragment of a bouldering coral)
* Cell filtration: 
    * filter tissue slurry collection through 40um cell strainer (P. damicornis). 
    * For organisms with thicker mucus secretions, use a larger pore strainer first, such as 70 or 100um, followed by the 40um. Keep on ice. 
* Optional: 
    * dose with appropriate volume of RNAse inhibitor ([Millepore Sigma Protector RNAse Inhibitor](https://www.sigmaaldrich.com/US/en/product/roche/rnainhro))
    * 1 uL or 40 units per 1 mL of cell slurry

### Nuclei isolation

* Centrifuge tubes at 800 g at 4°C for 10 minutes. 
* Carefully remove the supernatant and gently resuspend the pellets in 1 mL of buffer 1 and 1 mL of buffer 2. 
* Transferred to a Dounce homogenizer and grind on ice for 4 minutes with pestle A, and let it rest on ice for 7 minutes. 
* Transfer liquid to corex tubes containing 8 mL of buffer 3, in a manner that the homogenate would form a layer on top of buffer 3. 
    * The two solutions have a different density, which allow them to stay one on top of the other without mixing. 
* Centrifuge for 20 mins at 4°C and 7,800 g, with the lowest break speed possible. 
    * This step allows separating coral nuclei, which will form a pellet at the bottom of the tube, from cell debris and Symbiodiniaceae cells, which stay at the interphase between the two buffers. 
* Completely remove the supernatant by pouring out of the tubes and then by aspiration with micropipette. 
* The pellets now contain coral nuclei and chromatin.

### Quality Control

Preliminary protocol taken from [Nadelmann et al. 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8191490/#S7title).

* FACS

    * Proceed to FACS (which should occur at least fifteen minutes after staining with NucBlue). 
    * Follow instructions for flow activated cell sorting to separate Hoechst-positive nuclei from small-particle debris 
    * Sort through the entire sample to maximize the yield of nuclei. 
    * Each sample should run through the FACS machine for approximately thirty minutes

* Microscopy

    * A standard light microscope can be used with a 10x or 20x lens to view the nuclei. 
    * Microscopes with automated platforms, such as a Keyence, can also be used, to observe the presence of the fluorescent DAPI label in the nuclei, as the sample had been previously stained with NucBlue. 
    * Check nuclei morphology to confirm presence of nuclei and removal of debris. An intact nucleus should be round and 8–12 μm in diameter (e.g. nuclei from heart tissue)

* RNA Extraction
    * Sion from the UM Oncogenomics core suggested to perform an RNA extraction of some of the isolated nuclei to ensure we have high quality RNA. 
    * We should also sequence a few too.