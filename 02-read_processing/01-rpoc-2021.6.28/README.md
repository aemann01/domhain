# Raw data processing

## Setup

Download raw data from ////

```bash
cd ~/domhain/02-read_processing/01-rpoc-2021.6.28 # on pickles
mkdir raw
cd raw
wget ////
cd ..
```

### 1. Install R packages (v4.1.0)

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("magrittr")
install.packages("stringr")
BiocManager::install("dada2")
install.packages("data.table")
install.packages("broom")
install.packages("qualpalr")
install.packages("seqinr")
```

NOTE: Problems with the stringi package? Try first installing it with conda (with the domhain environment activated)

```bash
conda install -c r r-stringi
```

Then update the conda environment

```bash
conda --update all
```

Load up R and try to load the stringi package

```R
library(stringi)
```

### 2. Load required libraries

```R
library(dada2)
library(stringr)
library(data.table)
library(broom)
library(qualpalr)
library(ShortRead)
library(Biostrings)
library(seqinr)
```

### 3. File path setup

```R
rawpath <- "raw"
wdpath <- "/home/allie/domhain/02-read_processing/01-rpoc-2021.6.28/" # change to where git repository was cloned
fnFs <- sort(list.files(rawpath, pattern="_R1_001.fastq.gz", full.names=T))
fnRs <- sort(list.files(rawpath, pattern="_R2_001.fastq.gz", full.names=T))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
```

```text
  [1] "DM000013V1PQ83"   "DM00001V1PQ1"     "DM00002V1PQ"
  [4] "DM00004V1PQ55"    "DM00004V1PQ65"    "DM00005V1PQ36"
  [7] "DM00005V1PQ46"    "DM00006V1PQ1"     "DM00007V1PQ1"
 [10] "DM00008V1PQ16"    "DM00008V1PQ26"    "DM00009V1PQ55"
 [13] "DM00009V1PQ65"    "DM00010V1PQ54"    "DM00010V1PQ64-65"
 [16] "DM00011V1PQ55"    "DM00012V1PQ1"     "DM00013V1PQ65"
 [19] "DM00013V1PQ73"    "DM00014V1PQ3"     "DM00014V1PQ64"
 [22] "DM00014V1PQ84"    "DM00015V1PQ"      "DM00016V1PQ16"
 [25] "DM00016V1PQ26"    "DM00017V1PQ16"    "DM00017V1PQ31"
 [28] "DM00018V1PQ46"    "DM00018V1PQ74"    "DM00018V1PQ75"
 [31] "DM00019V1PQ16"    "DM00019V1PQ36"    "DM00020V1PQ16"
 [34] "DM00020V1PQ26"    "DM00021V1PQ11-12" "DM00021V1PQ31-41"
 [37] "DM00021V1PQ55"    "DM00021V1PQ65"    "DM00021V1PQ74"
 [40] "DM00021V1PQ75"    "DM00021V1PQ84"    "DM00021V1PQ85"
 [43] "DM00022V1PQ16"    "DM00022V1PQ26"    "DM00023V1PQ16"
 [46] "DM00023V1PQ26"    "DM00024V1PQ16"    "DM00024V1PQ41"
 [49] "DM00024V1PQ46"    "DM00024V1PQ54"    "DM00024V1PQ55"
 [52] "DM00024V1PQ64"    "DM00024V1PQ65"    "DM00024V1PQ75"
 [55] "DM00025V1PQ31-41" "DM00026V1PQ65"    "DM00026V1PQ85"
 [58] "DM00027V1PQ55"    "DM00027V1PQ65"    "DM00028V1PQ16"
 [61] "DM00028V1PQ26"    "DM00029V1PQ26"    "DM00032V1PQ55"
 [64] "DM00032V1PQ84"    "DM00032V1PQ85"    "DM00033V1PQ26"
 [67] "DM00033V1PQ65"    "DM00034V1PQ55"    "DM00034V1PQ65"
 [70] "DM00035V1PQ16"    "DM00035V1PQ26"    "DM00036V1PQ51"
 [73] "DM00036V1PQ52"    "DM00036V1PQ55"    "DM00036V1PQ61"
 [76] "DM00036V1PQ62"    "DM00036V1PQ64"    "DM00036V1PQ65"
 [79] "DM00036V1PQ74"    "DM00036V1PQ81"    "DM00036V1PQ84"
 [82] "DM00036V1PQ85"    "DM00037V1PQ11"    "DM00037V1PQ42"
 [85] "DM00038V1PQ36"    "DM00038V1PQ85"    "DM00039V1PQ11"
 [88] "DM00039V1PQ41"    "DM00040V1PQ51"    "DM00040V1PQ81"
 [91] "DM00041V1PQ26"    "DM00041V1PQ36"    "DM00042V1PQ55-65"
 [94] "DM00042V1PQ75-85" "DM00043V1PQ63"    "DM00043V1PQ85"
 [97] "DM00044V1PQ16"    "DM00044V1PQ26"    "DM00044V1PQ55"
[100] "DM00045V1PQ65"    "DM00045V1PQ75"    "DM00046V1PQ26"
[103] "DM00046V1PQ45"    "DM00047V1PQ55"    "DM00047V1PQ65"
[106] "DM00048V1PQ51"    "DM00048V1PQ54"    "DM00048V1PQ62"
[109] "DM00048V1PQ74"    "DM00048V1PQ84"    "DM00049V1PQ11"
[112] "DM00049V1PQ41"    "DM00049V1PQ55"    "DM00049V1PQ65"
[115] "DM00050V1PQ11"    "DM00050V1PQ54"    "DM00050V1PQ55"
[118] "DM00050V1PQ65"    "DM00050V1PQ75"    "DM00050V1PQ85"
[121] "DM00051V1PQ11"    "DM00051V1PQ41"    "DM00052V1PQ55"
[124] "DM00053V1PQ75"    "DM00054V1PQ31"    "DM00054V1PQ63"
[127] "DM00055V1PQ55"    "DM00055V1PQ65"    "DM00055V1PQ75"
[130] "DM00055V1PQ85"    "DM00056V1PQ36"    "DM00057V1PQ16"
[133] "DM00057V1PQ75"    "DM00058V1PQ36"    "DM00058V1PQ74"
[136] "DM00059V1PQ55"    "DM00059V1PQ65"    "DM00059V1PQ74"
[139] "DM00060V1PQ26"    "DM00061V1PQ16"    "DM00061V1PQ65"
[142] "DM00062V1PQ55"    "DM00062V1PQ75"    "DM00062V1PQ84"
[145] "DM00063V1PQ55"    "DM00063V1PQ74"    "DM00063V1PQ75"
[148] "DM00064V1PQ16"    "DM00064V1PQ26"    "DM00065V1PQ16"
[151] "DM00065V1PQ26"    "DM00065V1PQ36"    "DM00066V1PQ16"
[154] "DM00066V1PQ26"    "DM00067V1PQ11"    "DM00067V1PQ41"
[157] "DM00067V1PQ75"    "DM00067V1PQ83"    "DM00067V1PQ84"
[160] "DM00068V1PQ53"    "DM00068V1PQ75"    "DM00069V1PQ53"
[163] "DM00069V1PQ83"    "DM00070V1PQ16"    "DM00070V1PQ42"
[166] "DM00071V1PQ16"    "DM00071V1PQ26"    "DM00072V1PQ11"
[169] "DM00072V1PQ31"    "DM00072V1PQ46"    "PCRBLANK1"
[172] "PCRBLANK2"        "PCRBLANK4"        "PCRBLANK5"
[175] "PCRBLANK6"        "PCRBLANK8"
```

### 4. Plot quality scores

```R
system("mkdir img")
pdf(paste(wdpath, "img/", "forward_quality_plot.pdf", sep=""))
plotQualityProfile(fnFs[10:20])
dev.off()
pdf(paste(wdpath, "img/", "reverse_quality_plot.pdf", sep=""))
plotQualityProfile(fnRs[10:20])
dev.off()
```

### 5. Preliminary filter (removes sequences with N's)

```R
fnFs.filtN <- file.path(rawpath, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(rawpath, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress = TRUE)
```

### 6. Primer removal

```R
cutadapt <- as.character(system("which cutadapt", intern=T))
system("cutadapt --version")
path.cut <- file.path(rawpath, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
FWD.RC <- dada2:::rc("MAYGARAARMGNATGYTNCARGA")
REV.RC <- dada2:::rc("GMCATYTGRTCNCCRTCRAA")
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", "MAYGARAARMGNATGYTNCARGA", "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", "GMCATYTGRTCNCCRTCRAA", "-A", FWD.RC) 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,"-o", fnFs.cut[i], "-p", fnRs.cut[i], fnFs.filtN[i], fnRs.filtN[i]))
}
cutFs <- sort(list.files(path.cut, pattern = "R1", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2", full.names = TRUE))
```

### 7. Filter and trim reads

```R
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, minLen = c(200,200), maxN=c(0,0), maxEE=c(3,4), rm.phix=TRUE, matchIDs=TRUE, compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
retained
```

```text
                                             percentage_retained
DM000013V1PQ83_2_S19_L001_R1_001.fastq.gz              92.631082
DM00001V1PQ1_S1_L001_R1_001.fastq.gz                   87.664474
DM00002V1PQ_S2_L001_R1_001.fastq.gz                    89.437535
DM00004V1PQ55_S3_L001_R1_001.fastq.gz                  85.945946
DM00004V1PQ65_S4_L001_R1_001.fastq.gz                  50.000000
DM00005V1PQ36_1_S5_L001_R1_001.fastq.gz                87.014726
DM00005V1PQ46_2_S6_L001_R1_001.fastq.gz                88.780237
DM00006V1PQ1_S7_L001_R1_001.fastq.gz                   71.237026
DM00007V1PQ1_S8_L001_R1_001.fastq.gz                   85.736329
DM00008V1PQ16_2_S9_L001_R1_001.fastq.gz                92.938502
DM00008V1PQ26_1_S10_L001_R1_001.fastq.gz               92.484485
DM00009V1PQ55_S11_L001_R1_001.fastq.gz                 86.700216
DM00009V1PQ65_S12_L001_R1_001.fastq.gz                 87.627529
DM00010V1PQ54_S13_L001_R1_001.fastq.gz                 88.356857
DM00010V1PQ64-65_S14_L001_R1_001.fastq.gz              87.350338
DM00011V1PQ55_1_S15_L001_R1_001.fastq.gz               86.418446
DM00012V1PQ1_S16_L001_R1_001.fastq.gz                  92.601655
DM00013V1PQ65_3_S18_L001_R1_001.fastq.gz               86.494772
DM00013V1PQ73_1_S17_L001_R1_001.fastq.gz               93.162001
DM00014V1PQ3_S20_L001_R1_001.fastq.gz                  89.714086
DM00014V1PQ64_2_S21_L001_R1_001.fastq.gz               90.870840
DM00014V1PQ84_1_S22_L001_R1_001.fastq.gz               90.687717
DM00015V1PQ_S23_L001_R1_001.fastq.gz                   88.601917
DM00016V1PQ16_S25_L001_R1_001.fastq.gz                 90.934368
DM00016V1PQ26_S26_L001_R1_001.fastq.gz                 90.173539
DM00017V1PQ16_S27_L001_R1_001.fastq.gz                 88.090264
DM00017V1PQ31_S28_L001_R1_001.fastq.gz                 88.565780
DM00018V1PQ46_1_S29_L001_R1_001.fastq.gz               90.793173
DM00018V1PQ74_1_S30_L001_R1_001.fastq.gz               88.489011
DM00018V1PQ75_1_S31_L001_R1_001.fastq.gz               87.570167
DM00019V1PQ16_S32_L001_R1_001.fastq.gz                 87.538604
DM00019V1PQ36_1_S33_L001_R1_001.fastq.gz               88.671409
DM00020V1PQ16_1_S34_L001_R1_001.fastq.gz               89.959027
DM00020V1PQ26_1_S35_L001_R1_001.fastq.gz               90.534769
DM00021V1PQ11-12_1_S37_L001_R1_001.fastq.gz            92.484016
DM00021V1PQ31-41_1_S38_L001_R1_001.fastq.gz            90.949921
DM00021V1PQ55_1_S39_L001_R1_001.fastq.gz               86.885398
DM00021V1PQ65_1_S40_L001_R1_001.fastq.gz               90.115282
DM00021V1PQ74_1_S41_L001_R1_001.fastq.gz               93.292303
DM00021V1PQ75_1_S36_L001_R1_001.fastq.gz               89.647318
DM00021V1PQ84_1_S42_L001_R1_001.fastq.gz               92.717977
DM00021V1PQ85_1_S43_L001_R1_001.fastq.gz               92.385014
DM00022V1PQ16_1_S44_L001_R1_001.fastq.gz               88.147626
DM00022V1PQ26_1_S45_L001_R1_001.fastq.gz               92.133274
DM00023V1PQ16_1_S46_L001_R1_001.fastq.gz               91.261154
DM00023V1PQ26_1_S47_L001_R1_001.fastq.gz               86.877287
DM00024V1PQ16_1_S49_L001_R1_001.fastq.gz               89.015781
DM00024V1PQ41_1_S50_L001_R1_001.fastq.gz               93.034754
DM00024V1PQ46_1_S51_L001_R1_001.fastq.gz               84.291498
DM00024V1PQ54_1_S52_L001_R1_001.fastq.gz               87.378424
DM00024V1PQ55_1_S53_L001_R1_001.fastq.gz               90.926361
DM00024V1PQ64_1_S54_L001_R1_001.fastq.gz               91.984428
DM00024V1PQ65_1_S55_L001_R1_001.fastq.gz               87.695883
DM00024V1PQ75_1_S56_L001_R1_001.fastq.gz               89.801903
DM00025V1PQ31-41_1_S57_L001_R1_001.fastq.gz            90.008846
DM00026V1PQ65_1_S58_L001_R1_001.fastq.gz               90.764517
DM00026V1PQ85_1_S59_L001_R1_001.fastq.gz               82.332155
DM00027V1PQ55_1_S60_L001_R1_001.fastq.gz               83.793055
DM00027V1PQ65_1_S61_L001_R1_001.fastq.gz               88.734743
DM00028V1PQ16_1_S62_L001_R1_001.fastq.gz               88.581829
DM00028V1PQ26_1_S63_L001_R1_001.fastq.gz               82.725986
DM00029V1PQ26_1_S64_L001_R1_001.fastq.gz               89.578044
DM00032V1PQ55_1_S73_L001_R1_001.fastq.gz               91.912121
DM00032V1PQ84_1_S74_L001_R1_001.fastq.gz               89.694183
DM00032V1PQ85_1_S75_L001_R1_001.fastq.gz               88.991579
DM00033V1PQ26_S76_L001_R1_001.fastq.gz                 84.500337
DM00033V1PQ65_S77_L001_R1_001.fastq.gz                 89.662648
DM00034V1PQ55_S78_L001_R1_001.fastq.gz                 87.956332
DM00034V1PQ65_S79_L001_R1_001.fastq.gz                 82.510135
DM00035V1PQ16_S81_L001_R1_001.fastq.gz                 92.565730
DM00035V1PQ26_S80_L001_R1_001.fastq.gz                 88.372769
DM00036V1PQ51_S82_L001_R1_001.fastq.gz                 90.057320
DM00036V1PQ52_S158_L001_R1_001.fastq.gz                88.145874
DM00036V1PQ55_S149_L001_R1_001.fastq.gz                89.323326
DM00036V1PQ61_S139_L001_R1_001.fastq.gz                88.407752
DM00036V1PQ62_S130_L001_R1_001.fastq.gz                87.657505
DM00036V1PQ64_S121_L001_R1_001.fastq.gz                89.475772
DM00036V1PQ65_S111_L001_R1_001.fastq.gz                88.719229
DM00036V1PQ74_S102_L001_R1_001.fastq.gz                89.664668
DM00036V1PQ81_S83_L001_R1_001.fastq.gz                 85.792464
DM00036V1PQ84_S84_L001_R1_001.fastq.gz                 81.296024
DM00036V1PQ85_S159_L001_R1_001.fastq.gz                89.478231
DM00037V1PQ11_1_S150_L001_R1_001.fastq.gz              84.296454
DM00037V1PQ42_1_S140_L001_R1_001.fastq.gz              80.808474
DM00038V1PQ36_1_S131_L001_R1_001.fastq.gz              90.562065
DM00038V1PQ85_1_S122_L001_R1_001.fastq.gz              89.671774
DM00039V1PQ11_1_S112_L001_R1_001.fastq.gz              91.694946
DM00039V1PQ41_1_S103_L001_R1_001.fastq.gz              79.926082
DM00040V1PQ51_1_S93_L001_R1_001.fastq.gz               89.833128
DM00040V1PQ81_1_S92_L001_R1_001.fastq.gz               85.979711
DM00041V1PQ26_1_S151_L001_R1_001.fastq.gz              89.917185
DM00041V1PQ36_1_S160_L001_R1_001.fastq.gz              90.886483
DM00042V1PQ55-65_1_S141_L001_R1_001.fastq.gz           90.244805
DM00042V1PQ75-85_1_S132_L001_R1_001.fastq.gz           89.039130
DM00043V1PQ63_1_S123_L001_R1_001.fastq.gz              77.668012
DM00043V1PQ85_1_S113_L001_R1_001.fastq.gz              86.913691
DM00044V1PQ16_S104_L001_R1_001.fastq.gz                84.498092
DM00044V1PQ26_S94_L001_R1_001.fastq.gz                 81.943309
DM00044V1PQ55_S85_L001_R1_001.fastq.gz                 87.465723
DM00045V1PQ65_S152_L001_R1_001.fastq.gz                89.938132
DM00045V1PQ75_S142_L001_R1_001.fastq.gz                85.370207
DM00046V1PQ26_S133_L001_R1_001.fastq.gz                85.490742
DM00046V1PQ45_S124_L001_R1_001.fastq.gz                83.747873
DM00047V1PQ55_S114_L001_R1_001.fastq.gz                90.034532
DM00047V1PQ65_S105_L001_R1_001.fastq.gz                90.996986
DM00048V1PQ51_S95_L001_R1_001.fastq.gz                 82.336177
DM00048V1PQ54_S86_L001_R1_001.fastq.gz                 88.323486
DM00048V1PQ62_S153_L001_R1_001.fastq.gz                91.521526
DM00048V1PQ74_S143_L001_R1_001.fastq.gz                87.875657
DM00048V1PQ84_S134_L001_R1_001.fastq.gz                91.488695
DM00049V1PQ11_S125_L001_R1_001.fastq.gz                89.128561
DM00049V1PQ41_S115_L001_R1_001.fastq.gz                83.540133
DM00049V1PQ55_S106_L001_R1_001.fastq.gz                88.880597
DM00049V1PQ65_S97_L001_R1_001.fastq.gz                 89.462772
DM00050V1PQ11_S87_L001_R1_001.fastq.gz                 81.577769
DM00050V1PQ54_S154_L001_R1_001.fastq.gz                91.645287
DM00050V1PQ55_S145_L001_R1_001.fastq.gz                89.212883
DM00050V1PQ65_S135_L001_R1_001.fastq.gz                85.875878
DM00050V1PQ75_S126_L001_R1_001.fastq.gz                81.836244
DM00050V1PQ85_S116_L001_R1_001.fastq.gz                84.245551
DM00051V1PQ11_S107_L001_R1_001.fastq.gz                89.247312
DM00051V1PQ41_S98_L001_R1_001.fastq.gz                 86.265092
DM00052V1PQ55_S88_L001_R1_001.fastq.gz                 86.334780
DM00053V1PQ75_S155_L001_R1_001.fastq.gz                86.133069
DM00054V1PQ31_S146_L001_R1_001.fastq.gz                88.635432
DM00054V1PQ63_S136_L001_R1_001.fastq.gz                87.531958
DM00055V1PQ55_S117_L001_R1_001.fastq.gz                85.263055
DM00055V1PQ65_S127_L001_R1_001.fastq.gz                85.638756
DM00055V1PQ75_S108_L001_R1_001.fastq.gz                82.799016
DM00055V1PQ85_S99_L001_R1_001.fastq.gz                 80.034090
DM00056V1PQ36_S89_L001_R1_001.fastq.gz                 86.174532
DM00057V1PQ16_S156_L001_R1_001.fastq.gz                88.611485
DM00057V1PQ75_S147_L001_R1_001.fastq.gz                86.984111
DM00058V1PQ36_S137_L001_R1_001.fastq.gz                91.724670
DM00058V1PQ74_S128_L001_R1_001.fastq.gz                92.239754
DM00059V1PQ55_S118_L001_R1_001.fastq.gz                89.681076
DM00059V1PQ65_S109_L001_R1_001.fastq.gz                85.464171
DM00059V1PQ74_S100_L001_R1_001.fastq.gz                76.034115
DM00060V1PQ26_S90_L001_R1_001.fastq.gz                 90.277915
DM00061V1PQ16_S157_L001_R1_001.fastq.gz                87.022655
DM00061V1PQ65_S148_L001_R1_001.fastq.gz                86.514548
DM00062V1PQ55_S138_L001_R1_001.fastq.gz                84.896339
DM00062V1PQ75_S129_L001_R1_001.fastq.gz                90.649608
DM00062V1PQ84_S119_L001_R1_001.fastq.gz                85.217025
DM00063V1PQ55_S110_L001_R1_001.fastq.gz                83.078270
DM00063V1PQ74_S101_L001_R1_001.fastq.gz                78.302262
DM00063V1PQ75_S91_L001_R1_001.fastq.gz                 89.739535
DM00064V1PQ16_S169_L001_R1_001.fastq.gz                90.205388
DM00064V1PQ26_S170_L001_R1_001.fastq.gz                89.966661
DM00065V1PQ16_S171_L001_R1_001.fastq.gz                88.360165
DM00065V1PQ26_S172_L001_R1_001.fastq.gz                87.429591
DM00065V1PQ36_S173_L001_R1_001.fastq.gz                81.868031
DM00066V1PQ16_S174_L001_R1_001.fastq.gz                86.518013
DM00066V1PQ26_S175_L001_R1_001.fastq.gz                87.834692
DM00067V1PQ11_S176_L001_R1_001.fastq.gz                81.344226
DM00067V1PQ41_S177_L001_R1_001.fastq.gz                83.578659
DM00067V1PQ75_S178_L001_R1_001.fastq.gz                82.168150
DM00067V1PQ83_S179_L001_R1_001.fastq.gz                81.490433
DM00067V1PQ84_S180_L001_R1_001.fastq.gz                87.896290
DM00068V1PQ53_S181_L001_R1_001.fastq.gz                81.144112
DM00068V1PQ75_S182_L001_R1_001.fastq.gz                85.918433
DM00069V1PQ53_S183_L001_R1_001.fastq.gz                90.274900
DM00069V1PQ83_S184_L001_R1_001.fastq.gz                88.340238
DM00070V1PQ16_S185_L001_R1_001.fastq.gz                89.009532
DM00070V1PQ42_S186_L001_R1_001.fastq.gz                89.509517
DM00071V1PQ16_S187_L001_R1_001.fastq.gz                87.372922
DM00071V1PQ26_S188_L001_R1_001.fastq.gz                91.618056
DM00072V1PQ11_S189_L001_R1_001.fastq.gz                84.979424
DM00072V1PQ31_S190_L001_R1_001.fastq.gz                89.751939
DM00072V1PQ46_S191_L001_R1_001.fastq.gz                86.819681
PCRBLANK1_S24_L001_R1_001.fastq.gz                     61.538462
PCRBLANK2_S48_L001_R1_001.fastq.gz                     83.333333
PCRBLANK4_S96_L001_R1_001.fastq.gz                      6.811146
PCRBLANK5_S120_L001_R1_001.fastq.gz                    80.273973
PCRBLANK6_S144_L001_R1_001.fastq.gz                     2.898551
PCRBLANK8_S192_L001_R1_001.fastq.gz                     8.571429
```

### 8. Learn and plot error rates

```R
set.seed(12349)
errF <- learnErrors(filtFs, multithread=T, random=T)
errR <- learnErrors(filtRs, multithread=T, random=T)
png(paste(wdpath, "img/", "error_plot.png", sep=""))
plotErrors(errF, nominalQ=TRUE) 
dev.off()
```

### 9. Dereplication

```R
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

### 10. Sample inference

```R
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

### 11. Filter out samples with fewer than 5000 reads

```R
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(derepFs, getN), sapply(derepRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN))
samples_to_keep <- track[,4] > 5000
samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)]
paste(samples_to_remove)
```

```text
 [1] "DM00002V1PQ"      "DM00004V1PQ55"    "DM00004V1PQ65"    "DM00005V1PQ36"
 [5] "DM00005V1PQ46"    "DM00009V1PQ55"    "DM00010V1PQ64-65" "DM00013V1PQ73"
 [9] "DM00024V1PQ46"    "DM00026V1PQ85"    "DM00059V1PQ74"    "PCRBLANK1"
[13] "PCRBLANK2"        "PCRBLANK4"        "PCRBLANK5"        "PCRBLANK6"
[17] "PCRBLANK8"
```

### 12. Merge paired end reads

```R
mergers <- mergePairs(dadaFs[samples_to_keep], derepFs[samples_to_keep], dadaRs[samples_to_keep], derepRs[samples_to_keep], verbose=T)
```

### 13. Construct sequence table

```R
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

```text
[1]   159 14893
```

### 14. Sequence length distribution plot

```R
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab))))
png(paste(wdpath, "img/", "length_hist.png", sep=""))
plot(x=length.histogram[,1], y=length.histogram[,2])
dev.off()
```

### 15. Remove chimeras

```R
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=T, verbose=T)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
```

```text
[1]   159 11657
[1] 0.9555989
```

### 16. Processing summary

```R
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep,], sapply(dadaFs[samples_to_keep], getN), sapply(dadaRs[samples_to_keep], getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "non_chimeric")
rownames(track) <- sample.names[samples_to_keep]
track
```

### 17. Save output

```R
write.table(data.frame("row_names"=rownames(track),track),"read_retention.txt", row.names=FALSE, quote=F, sep="\t")
uniquesToFasta(seqtab.nochim, "rep_set.fa")
system("awk '/^>/{print \">ASV\" ++i; next}{print}' < rep_set.fa > rep_set_fix.fa")
system("mv rep_set_fix.fa rep_set.fa")
```

### 18. Clean up ASV names

```R
my_otu_table <- t(as.data.frame(seqtab.nochim)) 
ASV.seq <- as.character(unclass(row.names(my_otu_table))) 
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') 
colnames(seqtab.nochim) <- ASV.num 
write.table(data.frame("row_names"=rownames(seqtab.nochim),seqtab.nochim),"sequence_table.merged.txt", row.names=FALSE, quote=F, sep="\t")
```

### 19. Assign taxonomy

Get database

```bash
mkdir refdb
cd refdb
update_blast.pl 
```

```bash
cat rep_set.fa | parallel --block 100k --recstart '>' --pipe blastn -db ~/refdb/nt/nt -evalue 1e-10 -outfmt 6 -query - > rep_set.blast.out
```

Taxonomy assigned from blast output with MEGAN Community Edition (v.6). Select -> All Nodes, File -> Export -> Export to CSV -> readName_to_taxonPath, assigned, tab

Clean up taxonomy file

```bash
sed -i 's/"NCBI;cellular organisms;//' taxonomy.txt
sed -i 's/"//g' taxonomy.txt
```

Remove unwanted ASVs from taxonomy, sequence table

```bash
grep -v "Bacteria" taxonomy.txt | awk '{print $1}' > unwanted.ids
grep "Bacteria" taxonomy.txt | awk '{print $1}' > wanted.ids
seqtk subseq rep_set.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt > taxonomy_bac.txt
sed -i 's/;$//' taxonomy_bac.txt
python fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
```

Build a tree from the filtered representative sequences

```bash
mafft rep_set.filt.fa > rep_set.align.fa
fasttree -nt rep_set.align.fa > rep_set.align.tre
```

