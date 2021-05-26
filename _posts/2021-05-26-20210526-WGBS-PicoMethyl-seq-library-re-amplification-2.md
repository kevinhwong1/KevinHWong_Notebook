---
layout: post
title: 20210526 WGBS PicoMethyl-seq library re-amplification 2
date: '2021-05-26'
categories: Processing
tags: Porites astreoides, WGBS, PicoMethyl-seq
---

**Project:** Thermal Transplant Molecular, Holobiont Integration (ES)

### Goal

Re-amplify PicoMethyl-seq libraries that had low output (<5ng/ul) from initial prep using [this protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/Thermal-Transplant-WGBS-PicoMethyl-Protocol/). I am also re-amplifying samples for E. Strand that had low output from the PicoMethyl-seq prep.

### Protocol

Repeating the last two steps of [this protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/Thermal-Transplant-WGBS-PicoMethyl-Protocol/):

1) Second amplification with **4 cycles** instead of 10 cycles

2) 1X bead clean up

For all samples, I used the following inputs for the second amplification:
- 12 uL of the PicoMethyl-seq library
- 13 uL of Library Amp Master Mix (2x)
- 1 uL of the combined i5 and i7 index primers (I used the same indices as the original library prep)

#### Samples, indices, original Qubit quantification, post re-amplification Qubit quantification


| KW/ES |   Sample #  |  Vial  | i5/i7 index # | Original Library  (ng/uL) | DNA 1  (ng/uL) | DNA 2  (ng/uL) |
|:-----:|:-----------:|:------:|:-------------:|:-------------------------:|:--------------:|:--------------:|
|   KW  |      49     | 18-358 |       49      |            5.96           |      20.4      |      20.4      |
|   KW  |      52     |  L-924 |       52      |            2.50           |      26.6      |      26.4      |
|   ES  |      7      |  1707  |       7       |            2.78           |      12.7      |      12.6      |
|   ES  |      8      |  2212  |       8       |            5.94           |      25.0      |      24.8      |
|   ES  |      10     |  2861  |       10      |            LOW            |       LOW      |       LOW      |
|   ES  |      15     |  1103  |       15      |            LOW            |       LOW      |       LOW      |
|       | Qubit std.1 |        |               |                           |     178.04     |                |
|       | Qubit std.2 |        |               |                           |    19128.66    |                |


#### Tapestation Results

[D5000 Tapestation protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/DNA-Tapestation/)

[Full Tapestation Report](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/Tapestation_Results/2021-05-26_tapestation.pdf)


![#49]({{ site.baseurl}}/images/20210526_49.png "Summary")

![#52]({{ site.baseurl}}/images/20210526_52.png "Summary")

![#7]({{ site.baseurl}}/images/20210526_7.png "Summary")

![#8]({{ site.baseurl}}/images/20210526_8.png "Summary")
