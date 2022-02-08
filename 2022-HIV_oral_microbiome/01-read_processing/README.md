# Raw data processing

## Setup

Download raw data from ////

```bash
# cd /home/allie/domhain/02-read_processing # on hillary
# mkdir 02-raw_fastq && cd 02-raw_fastq
# wget ////
# cd ..
```

Activate conda environment

```bash
conda activate domhain
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
rawpath <- "02-raw_fastq"
wdpath <- "/home/allie/domhain/02-read_processing/" # change to where git repository was cloned
fnFs <- sort(list.files(rawpath, pattern="_R1_001.fastq.gz", full.names=T))
fnRs <- sort(list.files(rawpath, pattern="_R2_001.fastq.gz", full.names=T))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
# check for duplicate sample names
duplicated(sample.names)
```

### 4. Plot quality scores

```R
system("mkdir img")
pdf(paste(wdpath, "img/", "forward_quality_plot.pdf", sep=""))
plotQualityProfile(fnFs[10:25])
dev.off()
pdf(paste(wdpath, "img/", "reverse_quality_plot.pdf", sep=""))
plotQualityProfile(fnRs[10:25])
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
  system2(cutadapt, args = c("--cores=0", R1.flags, R2.flags, "-n", 2,"-o", fnFs.cut[i], "-p", fnRs.cut[i], fnFs.filtN[i], fnRs.filtN[i]))
}
cutFs <- sort(list.files(path.cut, pattern = "R1", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2", full.names = TRUE))
```

### 7. Filter and trim reads

```R
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, trimRight=25, maxN=c(0,0), maxEE=c(4,6), rm.phix=TRUE, matchIDs=TRUE, compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
retained
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
# reassign sample names
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

### 10. Sample inference

```R
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

<!-- Do this at the analysis step - keep all samples for now
### 11. Filter out samples with fewer than 1000 reads

```R
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(derepFs, getN), sapply(derepRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN))
samples_to_keep <- track[,4] > 1000
samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)]
paste(samples_to_remove)
```

```text
 [1] "DM00004V1PQ55"    "DM00004V1PQ65"    "DM00026V1PQ85-1"  "DM00148V1PQ75"
 [5] "DM00153V1PQ75"    "DM00165V1PQ75"    "DM00189V1PQ53"    "DM00217V1PQ65"
 [9] "DM00241V1PQ46"    "DM00241V1PQ74"    "DM00395V1PQ55"    "DM00445V1PQ54"
[13] "DM00446V1PQ26-65" "DM00447V1PQ16-55" "DM00448V1PQ55"    "DM00449V1PQ16"
[17] "DM00450V1PQ51"    "DM00451V1PQ16"    "DM00452V1PQ26"    "DM00453V1PQ16"
[21] "DM00454V1PQ32"    "DM00455V1PQ55"    "DM00456V1PQ16"    "DM00457V1PQ16"
[25] "DM00458V1PQ17"    "DM00459V1PQ55"    "DM00459V1PQ74"    "DM00459V1PQ84"
[29] "DM00460V1PQ82"    "DM00461V1PQ55"    "DM004623V1PQ75"   "DM00462V1PQ36"
[33] "DM00463V1PQ26"    "DM00463V1PQ74"    "DM00464V1PQ65"    "DM00464V1PQ75"
[37] "DM00465V1PQ55"    "DM00467V1PQ16"    "DM00468V1PQ26"    "DM00469V1PQ16"
[41] "DM00470V1PQ16"    "DM00543V1PQ75"    "DM00544V1PQ16"    "PCRBlank1"
[45] "PCRBlank10"       "PCRBlank11"       "PCRBlank14"       "PCRBlank16"
[49] "PCRBlank17"       "PCRBlank19"       "PCRBlank2"        "PCRBlank20"
[53] "PCRBlank21"       "PCRBlank23"       "PCRBlank24"       "PCRBlank25"
[57] "PCRBlank26"       "PCRBlank27"       "PCRBlank29"       "PCRBlank3"
[61] "PCRBlank31"       "PCRBlank4"        "PCRBlank5"        "PCRBlank6"
[65] "PCRBlank7"        "PCRBlank8"        "PCRBlank9"
``` -->

### 12. Merge paired end reads

```R
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=T)
```

### 13. Construct sequence table

```R
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

```text
[1]   969 67535
```

### 14. Length filter

```R
table(nchar(colnames(seqtab)))
```

```text
   56    62    63    64    65    66    67    68    69    70    71    72    73
    1    39    68    59    17    40    67    45    40    85    20    66    51
   74    75    76    77    78    79    80    81    82    83    84    85    86
   41    52    69    35    92    54    46    29    68    32    69   106    52
   87    88    89    90    91    92    93    94    95    96    97    98    99
   59    80    69    60   110    34    76   100    41    97    57   131    62
  100   101   102   103   104   105   106   107   108   109   110   111   112
   49    42    67   105    39    76    58    52    86    53    43    59    85
  113   114   115   116   117   118   119   120   121   122   123   124   125
   42    95    82    50    66    87    38    30    39    32   114   114    62
  126   127   128   129   130   131   132   133   134   135   136   137   138
   42   157    95    79    72    29    50    44    41    30    77   121    79
  139   140   141   142   143   144   145   146   147   148   149   150   151
  131    69    51    75    47    69    39    35   119   112    31    52   104
  152   153   154   155   156   157   158   159   160   161   162   163   164
   50    47    53    30    91   129    34    65    64    51    32    79    53
  165   166   167   168   169   170   171   172   173   174   175   176   177
   77    54    53    63    58    29    92   128    26    28   140    19    34
  178   179   180   181   182   183   184   185   186   187   188   189   190
   85    37    57    67    39    65    98    38    39    50    50    48   108
  191   192   193   194   195   196   197   198   199   200   201   202   203
   43   105    63    33    90    90    31    24    73    27    70    88    25
  204   205   206   207   208   209   210   211   212   213   214   215   216
   32    43    24    21    34    30    41    38    23    48    52    52    46
  217   218   219   220   221   222   223   224   225   226   227   228   229
  158    23    41   256    20    59   152    11    52    47    20    56    39
  230   231   232   233   234   235   236   237   238   239   240   241   242
   30    23    23    30   144    45    14    14    55    22    40    40    11
  243   244   245   246   247   248   249   250   251   252   253   254   255
   22    94    28    22    44    16    16    31    67    87    63    23    21
  256   257   258   259   260   261   262   263   264   265   266   267   268
  187    15    13    15     2    35    36    17    65    23    18    15    56
  269   270   271   272   273   274   275   276   277   278   279   280   281
   21     8    10    15    40    14    18    16    18    12    31    17     7
  282   283   284   285   286   287   288   289   290   291   292   293   294
    9    67    14    30    27    17    19    74    16    27    30     8    21
  295   296   297   298   299   300   301   302   303   304   305   306   307
   27     4    24    21    19    28    17    23    19    18    17    13    20
  308   309   310   311   312   313   314   315   316   317   318   319   320
   23    48    25    21    11    57    18    17    26     6     8    43    55
  321   322   323   324   325   326   327   328   329   330   331   332   333
    8    18    17    12    18    41    71    26    11    31    20     7    18
  334   335   336   337   338   339   340   341   342   343   344   345   346
   12     8     7    15     8    37    27    15    14    27    16    29    29
  347   348   349   350   351   352   353   354   355   356   357   358   359
    8    17    10    19    23    20    14    23    18    33    37     8    15
  360   361   362   363   364   365   366   367   368   369   370   371   372
   39    32    57    14     4    11    19    31     7    17     9     8     2
  373   374   375   376   377   378   379   380   381   382   383   384   385
   31    31     5    24    13    20     7     8    20    12     6    37    35
  386   387   388   389   390   391   392   393   394   395   396   397   398
   15    56    15    58    10    16     5     7    23     3    16    13    13
  399   400   401   402   403   404   405   406   407   408   409   410   411
   16    14    11    12    31    36    20    40     8    50     3    55    16
  412   413   414   415   416   417   418   420   421   422   423   424   425
   11     4    15     9    51    25    15   127    13    16     7    12    11
  426   427   428   429   430   431   432   433   434   435   436   437   438
  228    16    28    25     7     6    15    40    13   166     6    16    12
  439   440   441   442   443   444   445   446   447   448   449   450   451
   23    13    37     6    34    12    26    25    13     7     7    32     9
  452   453   454   455   456   457   458   459   460   461   462   463   464
    7     8     1     3    18    53     5    32    13    25    33     3     6
  465   466   467   468   469   470   471   472   473   474   475   476   477
    2    15     4     4     5    13    40    19     4     5     9     8    40
  478   479   480   481   482   483   484   485   486   487   488   489   490
49392     9     9    63    35    14   453     2     7   843     5     7    14
  491   492   493   494   495   496   497
    3     1   209     8     4     5     1
```

Filter out any ASVs under 450bp in length

```R
seqlens <- nchar(getSequences(seqtab))
seqtab.filt <- seqtab[,seqlens >= 450]
dim(seqtab.filt)
```

```text
[1]   969 51500
```

### 15. Sequence length distribution plot post filter

```R
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab.filt))))
png(paste(wdpath, "img/", "length_hist.png", sep=""))
plot(x=length.histogram[,1], y=length.histogram[,2])
dev.off()
```

### 15. Remove chimeras

```R
seqtab.nochim <- removeBimeraDenovo(seqtab.filt, method="pooled", multithread=T, verbose=T)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab.filt)
```

```text
Identified 28743 bimeras out of 51500 input sequences.
[1]   969 22757
[1] 0.921333
```

### 16. Processing summary

```R
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochimeras")
rownames(track) <- sample.names
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

Run kraken2

```bash
kraken2 --db ~/refdb/kraken_rpoc/ --threads 8 --use-names --output rep_set.kraken.out rep_set.fa --unclassified-out rep_set.unclassified.kraken.out --confidence 0.01
```

```text
22757 sequences (10.88 Mbp) processed in 0.643s (2124.7 Kseq/m, 1016.06 Mbp/m).
  22196 sequences classified (97.53%)
  561 sequences unclassified (2.47%)
```

Get taxonomic lineage information from tax ids

```bash
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
# cd refdb && mkdir ncbi_taxonomy
# cd ncbi_taxonomy
# wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
# tar xzf new_taxdump.tar.gz
# cd ../..
awk -F"\t" '{print $3}' rep_set.kraken.out | sed 's/^.*(//' | sed 's/taxid //' | sed 's/)//' > query
```

<!-- Add unidentified to taxonomy list (only do this once)

```bash
sed -i '1 i \0       |       unidentified    |               |' ~/refdb/ncbi_taxonomy/fullnamelineage.dmp
```

We also have to change taxid 81850 to 33958 (was merged on March 11, 2021). Open query file in nano and do CTL+W and replace string (only found once in file). NOTE! Depending on the version of the ncbi taxonomy file you pull and kraken version, you may have some issues with taxon names not being in the file. -->

Use queries to pull full lineage information (takes a little bit of time)

```bash
cat query | while read line; do grep -w -m 1 ^$line ~/refdb/ncbi_taxonomy/fullnamelineage.dmp | awk -F"|" '{print $3, $2}' | sed 's/\t//g' | sed 's/  / /g' | sed 's/cellular organisms; //' | sed 's/; /;/g' | sed 's/ /_/g'; done > lineages
```

Merge with ASV name to get taxonomy table

```bash
awk '{print $2}' rep_set.kraken.out > asvids
paste asvids lineages > taxonomy.txt
```

Remove unwanted ASVs from taxonomy, sequence table

```bash
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" | awk '{print $1}' > wanted.ids
seqtk subseq rep_set.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" > taxonomy_bac.txt
python ../00-scripts/fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
sed -i 's/;/\t/g' taxonomy_bac.txt
```

Build a tree from the filtered representative sequences

```bash
mafft --thread 55 rep_set.filt.fa > rep_set.align.fa
fasttree -nt rep_set.align.fa > rep_set.align.tre
```
