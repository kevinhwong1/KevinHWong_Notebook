---
title: "20210426 WGBS PicoMethyl-seq library prep for Thermal Transplant Trial 4"
date: 2021-04-27
type: "Lab work"
project: "Porites Thermal Transplant"
categories: ["Lab work", "Porites astreoides", "DNA methylation"]
protocols: ["wgbs-library-prep-picomethyl"]
aliases:
  - "/20210426-WGBS-PicoMethyl-seq-library-prep-for-Thermal-Transplant-Trial-4/"
---
**Project:** Thermal Transplant Molecular

### Goal

Fourth trial of WGBS library prep using the Zymo Pico Methyl-seq library prep kit and [this protocol](/notebook/2021-04-16-thermal-transplant-wgbs-picomethyl-protocol.md).

Sample dilution and mix calculations can be found on [this google sheet](https://docs.google.com/spreadsheets/d/1kthTxfiwn0RAWAQLLW3-pWBg5MBleQaFaEdjgEvHr58/edit#gid=0)

### Protocol

I am using a modified version of [this protocol](/notebook/2021-04-16-thermal-transplant-wgbs-picomethyl-protocol.md).

For all samples, I used the following modifications:
- Use the **combined i5 and i7 index primers (an input of 1uL)**
- Reduce the Library Amp Master Mix (2x) to **13uL** (was originally 14uL)


#### Samples and Indices

| Sample | Coral ID | Life Stage | i5 Index # | i7 Index # |
|:------:|:--------:|:----------:|:----------:|:----------:|
| 18-227 |  R-11-A  |    Adult   |     15     |     15     |
| 18-370 |  P-15-A  |    Adult   |     16     |     16     |
| 18-130 |  R-5-A   |    Adult   |     17     |     17     |
| L-1263 |  R-7-A   |   Larvae   |     18     |     18     |
| L-924  |  P-9-A   |   Larvae   |     19     |     19     |
| L-562  |  R-5-A   |   Larvae   |     20     |     20     |
| L-1059 |  R-15-A  |   Larvae   |     21     |     21     |
| 18-67  |  P-14-B  |    Adult   |     22     |     22     |

#### Qubit Results
To test DNA quantity: [Qubit](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2019-05-31-Qubit-Protocol.md)  

I used the DNA broad range assay.

I ran the qubit and tapestation on the following day (20210427). After the bead clean up step I stored the libraries in the fridge at 4C for ~16 hours.

|    Sample    | DNA 1 (ng/ul) | DNA 2 (ng/ul) |
|:------------:|:-------------:|:-------------:|
|  Standard 1  |     232.77    |               |
|  Standard 2  |    19299.53   |               |
|  15          |      LOW      |      LOW      |
|  16          |      21.0     |      20.2     |
|  17          |      5.94     |      5.68     |
|  18          |      14.3     |      14.2     |
|  19          |      LOW      |      LOW      |
|  20          |      3.36     |      3.12     |
|  21          |      4.74     |      4.64     |
|  22          |      19.3     |      18.6     |


#### Tapestation Results

[D5000 Tapestation protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/DNA-Tapestation/)

[Full Tapestation Report](/images/Tapestation_Results/2021-04-27_tapestation.pdf)

![#16](/images/20210427_16.png "Summary")

![#17](/images/20210427_17.png "Summary")

![#18](/images/20210427_18.png "Summary")

![#21](/images/20210427_21.png "Summary")

![#22](/images/20210427_22.png "Summary")
