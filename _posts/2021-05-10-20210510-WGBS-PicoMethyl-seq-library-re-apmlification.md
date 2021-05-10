---
layout: post
title: 20210510 WGBS PicoMethyl-seq library re-apmlification
date: '2021-05-10'
categories: Processing, Protocols
tags: Porites astreoides, WGBS, PicoMethyl-seq
---

**Project:** Thermal Transplant Molecular

### Goal

Re-amplify PicoMethyl-seq libraries that had low output (<5ng/ul) from initial prep using [this protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/Thermal-Transplant-WGBS-PicoMethyl-Protocol/).

### Protocol

Repeating the last two steps of [this protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/Thermal-Transplant-WGBS-PicoMethyl-Protocol/):

1) Second amplification with **4 cycles** instead of 10 cycles

2) 1X bead clean up


For all samples, I used the following inputs for the second amplification:
- 12 uL of the PicoMethyl-seq library
- 13 uL of Library Amp Master Mix (2x)
- 1 uL of the combined i5 and i7 index primers (I used the same indices as the original library prep)


#### Samples, indices, original Qubit quantification, post re-amplification Qubit quantification

|   Sample #  | i5/i7 Index # | Original Library (ng/uL) | Re-amp DNA #1 (ng/uL) | Re-amp DNA #2 (ng/uL) |
|:-----------:|:-------------:|:------------------------:|:---------------------:|:---------------------:|
|      15     |       15      |            LOW           |          9.56         |          9.32         |
|      17     |       17      |           5.81           |          14.2         |          13.9         |
|      19     |       19      |            LOW           |          LOW          |          LOW          |
|      20     |       20      |           3.24           |          7.12         |          7.00         |
|      21     |       21      |           4.69           |          16.5         |          16.1         |
|      26     |       26      |           3.54           |          19.1         |          18.6         |
|      27     |       27      |           5.44           |          20.6         |          19.8         |
|      28     |       28      |           2.86           |          17.5         |          17.2         |
|      32     |       32      |           5.57           |          LOW          |          LOW          |
|      45     |       45      |            5.4           |          21.2         |          20.8         |
| Qubit std.1 |               |                          |         178.09        |                       |
| Qubit std.2 |               |                          |        19930.40       |                       |

#### Tapestation Results

[D5000 Tapestation protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/DNA-Tapestation/)

[Full Tapestation Report 1](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/Tapestation_Results/2021-05-10_tapestation1.pdf)

[Full Tapestation Report 2](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/Tapestation_Results/2021-05-10_tapestation2.pdf)


![#15]({{ site.baseurl}}/images/20210510_15.png "Summary")

![#17]({{ site.baseurl}}/images/20210510_17.png "Summary")

![#20]({{ site.baseurl}}/images/20210510_20.png "Summary")

![#21]({{ site.baseurl}}/images/20210510_21.png "Summary")

![#26]({{ site.baseurl}}/images/20210510_26.png "Summary")

![#27]({{ site.baseurl}}/images/20210510_27.png "Summary")

![#28]({{ site.baseurl}}/images/20210510_28.png "Summary")

![#45]({{ site.baseurl}}/images/20210510_45.png "Summary")
