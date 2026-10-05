---
title: "20210512 WGBS PicoMethyl-seq library prep for Thermal Transplant Trial 8"
date: 2021-05-12
type: "Lab work"
project: "Porites Thermal Transplant"
categories: ["Lab work", "Porites astreoides", "DNA methylation"]
protocols: ["wgbs-library-prep-picomethyl"]
aliases:
  - "/20210512-WGBS-PicoMethyl-seq-library-prep-for-Thermal-Transplant-Trial-8/"
---
**Project:** Thermal Transplant Molecular

### Goal

Eighth trial of WGBS library prep using the Zymo Pico Methyl-seq library prep kit and [this protocol](/notebook/2021-04-16-thermal-transplant-wgbs-picomethyl-protocol.md). I am re-doing 3 samples that had a PCR bubbles due to over amplification from [this post](/notebook/2021-04-01-20210401-wgbs-picomethyl-seq-library-prep-for-thermal-transplant-trial-1.md) and 2 samples that did not reamplify from [this post](/notebook/2021-05-10-20210510-wgbs-picomethyl-seq-library-re-amplification.md).

Sample dilution and mix calculations can be found on [this google sheet](https://docs.google.com/spreadsheets/d/1kthTxfiwn0RAWAQLLW3-pWBg5MBleQaFaEdjgEvHr58/edit#gid=0)

### Protocol

I am using a modified version of [this protocol](/notebook/2021-04-16-thermal-transplant-wgbs-picomethyl-protocol.md).

For all samples, I used the following modifications:
- Use the **combined i5 and i7 index primers (an input of 1uL)**
- Reduce the Library Amp Master Mix (2x) to **13uL** (was originally 14uL)


#### Samples and Indices

| Sample | Coral ID | Life Stage | i5 Index # | i7 Index # |
|:------:|:--------:|:----------:|:----------:|:----------:|
| 18-358 |  P-10-A  |    Adult   |     49     |     49     |
| 18-20  |  P-6-B   |    Adult   |     50     |     50     |
| L-933  |  P-10-A  |    Larvae  |     51     |     51     |
| L-924  |  P-9-A   |    Larvae  |     52     |     52     |
| L-1257 |  P-19-A  |    Larvae  |     53     |     53     |


#### Qubit Results
To test DNA quantity: [Qubit](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2019-05-31-Qubit-Protocol.md)  

I used the DNA broad range assay.

|    Sample    | DNA 1 (ng/ul) | DNA 2 (ng/ul) |
|:------------:|:-------------:|:-------------:|
|  Standard 1  |     180.74    |               |
|  Standard 2  |    21047.36   |               |
|  49          |      5.96     |      5.80     |
|  50          |      11.9     |      11.7     |
|  51          |      18.1     |      18.2     |
|  52          |      2.50     |      2.46     |
|  53          |      20.2     |      19.7     |


#### Tapestation Results

[D5000 Tapestation protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/DNA-Tapestation/)

[Full Tapestation Report](/images/Tapestation_Results/2021-05-12_tapestation.pdf)

![#49](/images/20210512_49.png "Summary")

![#50](/images/20210512_50.png "Summary")

![#51](/images/20210512_51.png "Summary")

![#52](/images/20210512_52.png "Summary")

![#53](/images/20210512_53.png "Summary")
