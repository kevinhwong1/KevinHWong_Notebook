---
layout: post
title: 20210712 AH Tag-seq Plate Prep
date: '2021-07-12'
categories: Processing
tags: Tagseq, Mcap2020
---

**Project:** Mcap2020

### Goal

Dilute extracted RNA to 13 ug/uL in 50 uL and consolidate into one plate for Tag-seq. Molecular QC information can be found [here](https://docs.google.com/spreadsheets/d/1Ew0AOs88n1i6wEyuqbv3tsRXwMjhqyYEeZX5sACFdDY/edit#gid=1004871473).

### Samples

|  Sample | RNA ng/ul average | Sample Input for 13ng/ul (uL) | Dilutent Vol (uL) | Total Vol (uL) | Well  | Sample Name |
|:-------:|:-----------------:|:-----------------------------:|:-----------------:|:--------------:|:-----:|:-----------:|
|    D1   |        180        |              3.61             |       46.39       |       50       |   A1  |     AH1     |
|   D10   |        126        |              5.16             |       44.84       |       50       |   B1  |     AH2     |
|   D11   |        148        |              4.39             |       45.61       |       50       |   C1  |     AH3     |
|   D12   |        84.2       |              7.72             |       42.28       |       50       |   D1  |     AH4     |
|   D13   |        94.7       |              6.86             |       43.14       |       50       |   E1  |     AH5     |
|   D14   |        88.8       |              7.32             |       42.68       |       50       |   F1  |     AH6     |
|   D15   |        66.4       |              9.79             |       40.21       |       50       |   G1  |     AH7     |
|   D16   |         96        |              6.77             |       43.23       |       50       |   H1  |     AH8     |
|   D17   |        49.8       |             13.05             |       36.95       |       50       |   A2  |     AH9     |
|   D18   |        55.1       |             11.80             |       38.20       |       50       |   B2  |     AH10    |
|   D19   |        49.6       |             13.10             |       36.90       |       50       |   C2  |     AH11    |
|    D2   |        250        |              2.60             |       47.40       |       50       |   D2  |     AH12    |
|   D20   |        68.1       |              9.54             |       40.46       |       50       |   E2  |     AH13    |
|   D21   |        53.6       |             12.13             |       37.87       |       50       |   F2  |     AH14    |
|   D22   |        94.6       |              6.87             |       43.13       |       50       |   G2  |     AH15    |
|   D23   |        29.1       |             22.34             |       27.66       |       50       |   H2  |     AH16    |
|   D24   |        32.1       |             20.25             |       29.75       |       50       |   A3  |     AH17    |
|   D29   |        28.3       |             22.97             |       27.03       |       50       |   B3  |     AH18    |
|   D3    |        561        |              1.16             |       48.84       |       50       |   C3  |     AH19    |
|   D30   |        22.8       |             28.51             |       21.49       |       50       |   D3  |     AH20    |
|   D31   |        19.9       |             32.66             |       17.34       |       50       |   E3  |     AH21    |
|   D32   |       14.35       |             45.30             |        4.70       |       50       |   F3  |     AH22    |
|   D36   |         13        |             50.00             |        0.00       |       50       |   G3  |     AH23    |
|   D37   |        36.1       |             18.01             |       31.99       |       50       |   H3  |     AH24    |
|   D38   |        18.1       |             35.91             |       14.09       |       50       |   A4  |     AH25    |
|    D4   |        277        |              2.35             |       47.65       |       50       |   B4  |     AH26    |
|   D40   |        17.6       |             36.93             |       13.07       |       50       |   C4  |     AH27    |
|    D5   |       158.5       |              4.10             |       45.90       |       50       |   D4  |     AH28    |
|    D6   |        140        |              4.64             |       45.36       |       50       |   E4  |     AH29    |
|    D7   |        177        |              3.67             |       46.33       |       50       |   F4  |     AH30    |
|    D8   |        118        |              5.51             |       44.49       |       50       |   G4  |     AH31    |
|    D9   |        105        |              6.19             |       43.81       |       50       |   H4  |     AH32    |
| Plug_10 |        53.3       |             12.20             |       37.80       |       50       |   A5  |     AH33    |
| Plug_11 |        31.9       |             20.38             |       29.62       |       50       |   B5  |     AH34    |
|  Plug_2 |        40.9       |             15.89             |       34.11       |       50       |   C5  |     AH35    |
|  Plug_4 |        80.1       |              8.11             |       41.89       |       50       |   D5  |     AH36    |
|  Plug_5 |        77.5       |              8.39             |       41.61       |       50       |   E5  |     AH37    |
|  Plug_6 |        27.2       |             23.90             |       26.10       |       50       |   F5  |     AH38    |
|  Plug_9 |         75        |              8.67             |       41.33       |       50       |   G5  |     AH39    |

### Plate Map

![Plate Map]({{ site.baseurl}}/images/20210712_PlateMap.png "Plate Map")

### Protocol

1) Remove samples from the -80&deg;C freezer and place on a tube rack in an ice bucket to thaw.

2) Fill each well of the 96-well PCR plate with the corresponding diluent volume (I used Ultra Pure Water).

![Dilutent]({{ site.baseurl}}/images/ultra-pure-h2o-RNA.jpg "Dilutent")

3) Place plate onto a block frozen at -80&deg;C.

4) Once samples are thawed, spun down, and pipette mixed the sample before adding the appropriate volume into the corresponding well.

5) Seal plate with aluminum foil tape and seal each well. Will re-label the plate once we have the JA number and shipping information.

![Plate label]({{ site.baseurl}}/images/20210712_PlateLabel.jpg "Plate label")

6) Place remaining samples and plate back into the -80&deg;C freezer.

![Plate loc]({{ site.baseurl}}/images/20210712_PlateLoc.jpg "Plate loc")
