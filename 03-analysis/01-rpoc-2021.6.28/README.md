# Data analysis

Running on pickles

### Install required libraries

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("philr")
BiocManager::install("ggtree")
BiocManager::install("ALDEx2")
install.packages("vegan")
install.packages("RColorBrewer")
install.packages("UpSetR")
install.packages("ggdendro")
install.packages("ape")
install.packages("ggfortify")
install.packages("randomForest")
install.packages("rfUtilities")
install.packages("phytools")
install.packages("gridExtra")
```

### Load required libraries

```R
library(philr)
library(ggtree)
library(ALDEx2)
library(vegan)
library(RColorBrewer)
library(UpSetR)
library(ggdendro)
library(ape)
library(ggfortify)
library(randomForest)
library(rfUtilities)
library(phytools)
library(phyloseq)
library(gridExtra)
library(microbiome)
```

### Load data into R

```R
map <- read.table("~/domhain/01-metadata/sweets_fermented_foods.txt", sep="\t", header=T)
seqtab <- read.table("~/domhain/02-read_processing/01-rpoc-2021.6.28/sequence_table.merged.txt", header=T, row.names=1)
tax <- read.table("~/domhain/02-read_processing/01-rpoc-2021.6.28/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
tree <- read.tree("~/domhain/02-read_processing/01-rpoc-2021.6.28/rep_set.align.tre")
tree.root <- midpoint.root(tree)
```

### Check if metadata and sequence table have same samples

```R
notinmeta <- setdiff(row.names(seqtab), map$short_id)
notinraw <- setdiff(map$short_id, row.names(seqtab))
notinmeta
```

```text
character(0)
```

```R
notinraw
```

```text
 [1] "DM00002V1PQ"   "DM00004V1PQ55" "DM00004V1PQ65" "DM00005V1PQ36"
 [5] "DM00005V1PQ46" "DM00009V1PQ55" "DM00010V1PQ64" "DM00013V1PQ73"
 [9] "DM00024V1PQ46" "DM00026V1PQ85" "DM00029V1PQ31" "DM00030V1PQ31"
[13] "DM00030V1PQ55" "DM00031V1PQ26" "DM00031V1PQ46" "DM00032V1PQ26"
[17] "DM00032V1PQ46" "DM00045V1PQ55" "DM00048V1PQ61" "DM00050V1PQ36"
[21] "DM00053V1PQ55" "DM00056V1PQ85" "DM00059V1PQ74" "DM00060V1PQ75"
[25] "DM00063V1PQ84"
```

### Create phyloseq object

```R
rownames(map) <- map$short_id
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=F), sample_data(map), tax_table(as.matrix(tax)), tree.root)
ps.dat
```

```text
otu_table()   OTU Table:         [ 7539 taxa and 159 samples ]
sample_data() Sample Data:       [ 159 samples by 20 sample variables ]
tax_table()   Taxonomy Table:    [ 7539 taxa by 1 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 7539 tips and 7538 internal nodes ]
```

### PHILR transformation

```R
philr.dat <- transform_sample_counts(ps.dat, function(x) x+1) # add pseudocount of one to ASVs to avoid log-ratios calculated from zero
is.rooted(phy_tree(philr.dat)) # check that tree is rooted
# [1] TRUE
is.binary.tree(phy_tree(philr.dat)) #check that multichotomies are resolved in tree
# [1] TRUE
phy_tree(philr.dat) <- makeNodeLabel(phy_tree(philr.dat), method="number", prefix="n")
asv.table <- otu_table(philr.dat)
tree <- phy_tree(philr.dat)
metadata <- sample_data(philr.dat)
tax <- tax_table(philr.dat)
philr.t <- philr(asv.table, tree, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")
```

### Heirarchical cluster dendrogram from transformed data

```R
system("mkdir img")
cols <- brewer.pal(7, "Set2") # color palette
hc <- hclust(dist(philr.t), method="complete")
df2 <- data.frame(cluster=cutree(hc,5), states=factor(hc$labels, levels=hc$labels[hc$order])) # get cluster assocaited with each sample
write.table(df2, "philr_cluster.txt", quote=F, sep="\t", col.names=NA)
hcd <- as.dendrogram(hc)
dend_data <- dendro_data(hcd, type="rectangle")
p1 <- ggplot(dend_data$segments) + geom_segment(aes(x=x,y=y, xend=xend, yend=yend)) + theme_classic() + geom_text(data = dend_data$labels, aes(x, y, label = label, hjust = 1, angle = 90)) + ylim(-2,30) + xlab("") + ylab("") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
merge <- merge(df2, map, by.x=c("states"), by.y=c("short_id"))
p2 <- ggplot(merge, aes(states, y=1, fill=factor(merge$aliquot_type))) + geom_tile() + scale_fill_manual(values=cols) + scale_y_continuous(expand=c(0,0)) + theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), legend.position="none")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)
maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
pdf("img/philr_dendrogram_aliquot_type.pdf")
grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5,1/5))
dev.off()
# by hiv status
p2 <- ggplot(merge, aes(states, y=1, fill=factor(merge$study_group))) + geom_tile() + scale_fill_manual(values=cols) + scale_y_continuous(expand=c(0,0)) + theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), legend.position="none")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)
maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
pdf("img/philr_dendrogram_study_group.pdf")
grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5,1/5))
dev.off()
```

### PCA of PHILR distances

```R
philr.dist <- dist(philr.t, method="euclidean")
pca <- prcomp(as.matrix(philr.dist))
pdf("img/philr_screeplot.pdf")
screeplot(pca)
dev.off()
pdf("img/pca_aliquot_type.pdf")
autoplot(pca, data=sample_data(ps.dat), colour="aliquot_type") + theme_minimal() 
dev.off()
pdf("img/pca_study_group.pdf")
autoplot(pca, data=sample_data(ps.dat), colour="study_group") + theme_minimal() 
dev.off()
pdf("img/pca_sex.pdf")
autoplot(pca, data=sample_data(ps.dat), colour="sex") + theme_minimal() 
dev.off()
```

### Betadiversity 

PCA of PHILR distances

```R
philr.dist <- dist(philr.t, method="euclidean")
pca <- prcomp(as.matrix(philr.dist))
pdf("img/philr_screeplot.pdf")
screeplot(pca)
dev.off()
pdf("img/pca_aliquot_type.pdf")
autoplot(pca, data=sample_data(ps.dat), colour="aliquot_type") + theme_minimal() 
dev.off()
pdf("img/pca_study_group.pdf")
autoplot(pca, data=sample_data(ps.dat), colour="study_group") + theme_minimal() 
dev.off()
pdf("img/pca_sex.pdf")
autoplot(pca, data=sample_data(ps.dat), colour="sex") + theme_minimal() 
dev.off()
```

Bray curtis

```R
ord <- ordinate(ps.dat, "PCoA", "bray")
pdf("img/bdiv_bray.pdf")
ord.plt <- plot_ordination(ps.dat, ord, type="samples", color="aliquot_type", shape="study_group")
ord.plt + geom_point() + theme_minimal()
dev.off()
```

```R
ord <- ordinate(ps.dat, "PCoA", "unifrac")
pdf("img/bdiv_unifrac.pdf")
ord.plt <- plot_ordination(ps.dat, ord, type="samples", color="aliquot_type", shape="study_group")
ord.plt + geom_point() + theme_minimal()
dev.off()
```

### Upset plot

```R
merged <- merge(seqtab, map, by="row.names")
n <- ncol(seqtab) + 1

# by aliquot type
agg <- aggregate(merged[,2:n], by=list(merged$aliquot_type), FUN=sum) 
#remove columns with only zeros
agg <- agg[,colSums(agg !=0) > 0]
rownames(agg) <- agg$Group.1
#convert to presence absence table 
agg[agg>1] <- 1
agg <- data.frame(t(agg[,-1]))
pdf("img/upset_aliquot_type.pdf")
upset(agg, order.by="freq", mainbar.y.label="Number of ASVs", sets.x.label="Shared ASVs", mb.ratio = c(0.55, 0.45))
dev.off()

# by study group
agg <- aggregate(merged[,2:n], by=list(merged$study_group), FUN=sum) 
#remove columns with only zeros
agg <- agg[,colSums(agg !=0) > 0]
rownames(agg) <- agg$Group.1
#convert to presence absence table 
agg[agg>1] <- 1
agg <- data.frame(t(agg[,-1]))
pdf("img/upset_study_group.pdf")
upset(agg, order.by="freq", mainbar.y.label="Number of ASVs", sets.x.label="Shared ASVs", mb.ratio = c(0.55, 0.45))
dev.off()

# by sex
agg <- aggregate(merged[,2:n], by=list(merged$sex), FUN=sum) 
#remove columns with only zeros
agg <- agg[,colSums(agg !=0) > 0]
rownames(agg) <- agg$Group.1
#convert to presence absence table 
agg[agg>1] <- 1
agg <- data.frame(t(agg[,-1]))
pdf("img/upset_sex.pdf")
upset(agg, order.by="freq", mainbar.y.label="Number of ASVs", sets.x.label="Shared ASVs", mb.ratio = c(0.55, 0.45))
dev.off()
```

### Alpha diversity

```R
pdf("img/adiv_aliqout_type.pdf")
plot_richness(ps.dat, measures=c("Observed", "Shannon"), color="aliquot_type") + theme_minimal()
dev.off()
```

### PERMANOVA

```R
metadata <- as(sample_data(ps.dat), "data.frame")
adonis(philr.dist ~ aliquot_type, data=metadata)
```

```text
Call:
adonis(formula = philr.dist ~ aliquot_type, data = metadata)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

              Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
aliquot_type   6     201.5  33.581  1.5161 0.05647  0.001 ***
Residuals    152    3366.6  22.149         0.94353
Total        158    3568.1                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
adonis(philr.dist ~ sex, data=metadata)
```

```text
Call:
adonis(formula = philr.dist ~ sex, data = metadata)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
sex         2      55.7  27.847  1.2368 0.01561  0.069 .
Residuals 156    3512.4  22.515         0.98439
Total     158    3568.1                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
adonis(philr.dist ~ study_group, data=metadata)
```

```text
Call:
adonis(formula = philr.dist ~ study_group, data = metadata)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
study_group   2     108.2  54.120  2.4402 0.03034  0.001 ***
Residuals   156    3459.8  22.178         0.96966
Total       158    3568.1                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
adonis(philr.dist ~ age_y, data=metadata)
```

```text
Call:
adonis(formula = philr.dist ~ age_y, data = metadata)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
age_y       8     296.9  37.113  1.7018 0.08321  0.001 ***
Residuals 150    3271.2  21.808         0.91679
Total     158    3568.1                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
adonis(philr.dist ~ fizzy_drinks, data=metadata)
```

```text
Call:
adonis(formula = philr.dist ~ fizzy_drinks, data = metadata)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

              Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
fizzy_drinks   2      80.4  40.215  1.7988 0.02254  0.001 ***
Residuals    156    3487.7  22.357         0.97746
Total        158    3568.1                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Dispersion test

```R
dispr <- vegan::betadisper(philr.dist, phyloseq::sample_data(ps.dat)$aliquot_type)
dispr
```
```text
	Homogeneity of multivariate dispersions

Call: vegan::betadisper(d = philr.dist, group =
phyloseq::sample_data(ps.dat)$aliquot_type)

No. of Positive Eigenvalues: 158
No. of Negative Eigenvalues: 0

Average distance to median:
 CA-PD  CA-PE  CA-PF CAE-PE CAE-PF  CE-PF  CF-PF
 4.553  4.414  4.290  2.826  4.316  3.776  4.669

Eigenvalues for PCoA axes:
(Showing 8 of 158 eigenvalues)
 PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8
230.77 194.82 146.91 128.25  98.39  90.43  82.90  75.08
```

```R
pdf("img/betadispr_aliquot_type.pdf")
boxplot(dispr)
dev.off()
```











### Differential abundance 

First need to format for qiime2 (transpose sequence table beforehand, columns = samples, rows = ASV IDs)

```R
system("biom convert -i ../01-raw_data_processing/sequence_table.16s.filtered.tr.txt -o sequence_table.16s.filtered.biom --table-type="OTU table" --to-hdf5")
system("biom summarize-table -i sequence_table.16s.filtered.biom")
```
```text
Num samples: 61
Num observations: 2,989
Total count: 7,928,033
Table density (fraction of non-zero values): 0.047

Counts/sample summary:
 Min: 2,394.000
 Max: 328,813.000
 Median: 119,915.000
 Mean: 129,967.754
 Std. dev.: 66,955.954
 Sample Metadata Categories: None provided
 Observation Metadata Categories: None provided

Counts/sample detail:
NegCtrl: 2,394.000
S48A: 16,000.000
W07A: 44,641.000
W12E: 47,275.000
W11A: 50,547.000
S13E: 52,958.000
S07E: 59,115.000
S24A: 61,119.000
S33E: 64,722.000
W04E: 66,235.000
S09E: 67,009.000
W14E: 71,210.000
S31E: 77,538.000
S25E: 80,651.000
W24E: 80,996.000
W03A: 83,347.000
W13A: 85,035.000
S08A: 87,589.000
S10A: 92,886.000
S15E: 95,415.000
W28E: 97,112.000
S32A: 99,606.000
S49E: 100,502.000
Negctrl: 101,417.000
W30E: 107,583.000
S06A: 109,901.000
S16A: 110,355.000
S12A: 110,985.000
W09A: 116,494.000
S11E: 117,369.000
S26A: 119,915.000
W27A: 122,649.000
S04A: 123,121.000
W17A: 124,667.000
W20E: 126,772.000
W16E: 129,583.000
W31A: 131,180.000
S02E: 132,168.000
W26E: 135,470.000
W06E: 135,878.000
W29A: 138,578.000
W25A: 147,346.000
W10E: 147,675.000
S27E: 152,616.000
S19E: 162,642.000
S05E: 163,124.000
S30A: 166,900.000
W23A: 173,031.000
S20A: 182,490.000
S14A: 189,895.000
W32E: 190,175.000
S21E: 201,080.000
S17E: 210,433.000
S01A: 211,940.000
S18A: 219,883.000
W15A: 229,232.000
S37E: 249,545.000
Blank: 249,745.000
S35E: 278,400.000
S36A: 295,081.000
S34A: 328,813.000
```

```R
system("qiime tools import --input-path sequence_table.16s.filtered.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path sequence_table.16s.filtered.qza")
```

Now can run ALDEx plugin through qiime

```R
system("qiime feature-table filter-samples --i-table sequence_table.16s.filtered.qza --m-metadata-file map.txt --p-where "[Sample-type]='swab'" --o-filtered-table swab-feature-table.qza")
system("qiime aldex2 aldex2 --i-table swab-feature-table.qza --m-metadata-file map.txt --m-metadata-column Season --output-dir season_aldex")
system("qiime aldex2 effect-plot --i-table season_aldex/differentials.qza --o-visualization season_aldex/season")
system("qiime tools view season_aldex/season.qzv")
```

![aldex season](https://github.com/aemann01/necrobiome/blob/master/02-analysis/season_aldex/effect_plot.png)


```R
system("qiime aldex2 extract-differences --i-table season_aldex/differentials.qza --o-differentials season_aldex/season --p-sig-threshold 0.1 --p-effect-threshold 0 --p-difference-threshold 0")
system("qiime tools export --input-path season_aldex/season.qza --output-path season_aldex/")
system("awk '{print $1}' season_aldex/differentials.tsv | grep "ASV" | while read line; do grep -w $line tax_for_phyloseq.txt ; done > season_aldex/differentials.taxonomy.txt")
system("head season_aldex/differentials.taxonomy.txt")
```

```text
ASV4	Bacteria	Actinobacteria	Actinobacteria_c	Corynebacteriales	Corynebacteriaceae	Corynebacterium	Corynebacterium_unknown
ASV5	Bacteria	Firmicutes	Bacilli	Bacillales	Planococcaceae	Planococcaceae_unknown	Planococcaceae_unknown
ASV6	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV7	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV10	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV16	Bacteria	Proteobacteria	Betaproteobacteria	Burkholderiales	Burkholderiaceae	Paraburkholderia	Paraburkholderia_unknown
ASV17	Bacteria	Firmicutes	Clostridia	Clostridiales	Peptostreptococcaceae	Peptostreptococcaceae_unknown	Peptostreptococcaceae_unknown
ASV20	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV22	Bacteria	Actinobacteria	Actinobacteria_c	Corynebacteriales	Corynebacteriaceae	Corynebacterium	Corynebacterium_unknown
ASV23	Bacteria	Firmicutes	Bacilli	Bacillales	Planococcaceae	Planococcaceae_unknown	Planococcaceae_unknown
```

Low vs high temperature -- first need to filter out the na sample from mapping file (delete in excel) and then filter from qza

```R
system("qiime feature-table filter-samples --i-table swab-feature-table.qza --m-metadata-file map.filt.txt --o-filtered-table swab-feature-table.filt.qza")
system("qiime aldex2 aldex2 --i-table swab-feature-table.filt.qza --m-metadata-file map.filt.txt --m-metadata-column Temp_group_binary --output-dir temp-low-hi_aldex")
system("qiime aldex2 effect-plot --i-table temp-low-hi_aldex/differentials.qza --o-visualization temp-low-hi_aldex/temp-low-hi")
system("qiime tools view temp-low-hi_aldex/temp-low-hi.qzv")
```

![aldex season](https://github.com/aemann01/necrobiome/blob/master/02-analysis/temp-low-hi_aldex/effect_plot.png)

```R
system("qiime aldex2 extract-differences --i-table temp-low-hi_aldex/differentials.qza --o-differentials temp-low-hi_aldex/temp-low-hi --p-sig-threshold 0.1 --p-effect-threshold 0 --p-difference-threshold 0")
system("qiime tools export --input-path temp-low-hi_aldex/temp-low-hi.qza --output-path temp-low-hi_aldex/")
system("awk '{print $1}' temp-low-hi_aldex/differentials.tsv | grep "ASV" | while read line; do grep -w $line tax_for_phyloseq.txt ; done > temp-low-hi_aldex/differentials.taxonomy.txt")
system("head temp-low-hi_aldex/differentials.taxonomy.txt")
```

```text
ASV4	Bacteria	Actinobacteria	Actinobacteria_c	Corynebacteriales	Corynebacteriaceae	Corynebacterium	Corynebacterium_unknown
ASV5	Bacteria	Firmicutes	Bacilli	Bacillales	Planococcaceae	Planococcaceae_unknown	Planococcaceae_unknown
ASV6	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV10	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV16	Bacteria	Proteobacteria	Betaproteobacteria	Burkholderiales	Burkholderiaceae	Paraburkholderia	Paraburkholderia_unknown
ASV20	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV22	Bacteria	Actinobacteria	Actinobacteria_c	Corynebacteriales	Corynebacteriaceae	Corynebacterium	Corynebacterium_unknown
ASV23	Bacteria	Firmicutes	Bacilli	Bacillales	Planococcaceae	Planococcaceae_unknown	Planococcaceae_unknown
ASV26	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV32	Bacteria	Proteobacteria	Gammaproteobacteria	Enterobacterales	Morganellaceae	Providencia	Providencia_unknown
```

Clostridium seems to be very high in high temp, very low in low, plot these along temperature gradient. First merge with taxonomy. Open in excel and get values from clostridium.

```R
seqtab.nochim <- read.table("../01-raw_data_processing/sequence_table.16s.filtered.tr.txt", header=T, row.names=1)
taxa <- read.table("../01-raw_data_processing/taxonomy_L7.txt", header=F, row.names=1)
merged <- merge(seqtab.nochim, taxa, by=0)
write.table(data.frame("row_names"=rownames(merged),merged),"sequence_taxonomy_table.16s.merged.txt", row.names=FALSE, quote=F, sep="\t")
```

Test plot differentially abundant ASVs by temperature

```R
test <- otu_table(philr.dat)[,"ASV4"]
test.m <- merge(rawmetadata, test, by=0)
test.m$Temperature_C <- as.numeric(as.character(test.m$Temperature_C))
test.m <- test.m[!is.na(test.m$Temperature_C),]
test.m$Temp_C_bin <- cut(test.m$Temperature_C, breaks=10)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}
df2 <- data_summary(test.m, varname="ASV4", groupnames=c("Temperature_C"))
png(paste("imgs/", "test_ASV4.png", sep=""))
ggplot(df2, aes(x=as.factor(Temperature_C), y=ASV4)) + geom_bar(stat="identity", color="black", fill="white") + theme_minimal() + xlab("Temperature C") + ylab("IRL Transformed Read Counts") + geom_errorbar(aes(ymin=ASV4-sd, ymax=ASV4+sd), width=.2) + geom_point(test.m, mapping=aes(x=as.factor(Temperature_C), y=ASV4)) + geom_jitter()
dev.off()
```

![asv4 IRL counts temp](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/test_ASV4.png)

Now do a real plot based off of differentially abundant clades in phylofactor by temp group

```R
# get list of asvs to filter data
factor1 <- c("ASV5", "ASV19", "ASV2", "ASV80", "ASV23", "ASV85", "ASV212", "ASV122")
factor2 <- c("ASV7", "ASV125", "ASV20", "ASV335", "ASV471", "ASV10", "ASV194", "ASV6", "ASV93", "ASV26", "ASV35", "ASV38", "ASV34")
factor3 <- c("ASV14", "ASV9", "ASV1", "ASV56", "ASV25", "ASV3", "ASV67")
# convert philr.dat transformed otu table to dataframe
df <- otu_table(philr.dat)
df <- as.data.frame(df)
# filter data
df.fac1 <- df[,which((names(df) %in% factor1)==TRUE)]
df.fac2 <- df[,which((names(df) %in% factor2)==TRUE)]
df.fac3 <- df[,which((names(df) %in% factor3)==TRUE)]
# sum across all rows
df.fac1 <- as.data.frame(rowSums(df.fac1))
df.fac2 <- as.data.frame(rowSums(df.fac2))
df.fac3 <- as.data.frame(rowSums(df.fac3))
# merge data with metadata
df.fac1.m <- merge(rawmetadata, df.fac1, by=0)
df.fac2.m <- merge(rawmetadata, df.fac2, by=0)
df.fac3.m <- merge(rawmetadata, df.fac3, by=0)
# prep temperature category
df.fac1.m$Temperature_C <- as.numeric(as.character(df.fac1.m$Temperature_C))
df.fac2.m$Temperature_C <- as.numeric(as.character(df.fac2.m$Temperature_C))
df.fac3.m$Temperature_C <- as.numeric(as.character(df.fac3.m$Temperature_C))
df.fac1.m <- df.fac1.m[!is.na(df.fac1.m$Temperature_C),]
df.fac2.m <- df.fac2.m[!is.na(df.fac2.m$Temperature_C),]
df.fac3.m <- df.fac3.m[!is.na(df.fac3.m$Temperature_C),]
df.fac1.m$Temp_C_bin <- cut(df.fac1.m$Temperature_C, breaks=10)
df.fac2.m$Temp_C_bin <- cut(df.fac2.m$Temperature_C, breaks=10)
df.fac3.m$Temp_C_bin <- cut(df.fac3.m$Temperature_C, breaks=10)
# summarize data into bins
df.fac1.sum <- data_summary(df.fac1.m, varname="rowSums(df.fac1)", groupnames=c("Temperature_C"))
df.fac2.sum <- data_summary(df.fac2.m, varname="rowSums(df.fac2)", groupnames=c("Temperature_C"))
df.fac3.sum <- data_summary(df.fac3.m, varname="rowSums(df.fac3)", groupnames=c("Temperature_C"))
# clean up
colnames(df.fac1.sum) <- c("Temperature_C", "factor", "sd")
colnames(df.fac2.sum) <- c("Temperature_C", "factor", "sd")
colnames(df.fac3.sum) <- c("Temperature_C", "factor", "sd")
# plot
png(paste("imgs/", "fact1_barplot_tempG.png", sep=""))
ggplot(df.fac1.sum, aes(x=as.factor(Temperature_C), y=factor)) + geom_bar(stat="identity", color="black", fill="white") + theme_minimal() + xlab("Temperature C") + ylab("ILR Transformed Read Counts") + geom_errorbar(aes(ymin=factor-sd, ymax=factor+sd), width=.2) + geom_point(df.fac1.sum, mapping=aes(x=as.factor(Temperature_C), y=factor)) + geom_jitter()
dev.off()
png(paste("imgs/", "fact2_barplot_tempG.png", sep=""))
ggplot(df.fac2.sum, aes(x=as.factor(Temperature_C), y=factor)) + geom_bar(stat="identity", color="black", fill="white") + theme_minimal() + xlab("Temperature C") + ylab("ILR Transformed Read Counts") + geom_errorbar(aes(ymin=factor-sd, ymax=factor+sd), width=.2) + geom_point(df.fac2.sum, mapping=aes(x=as.factor(Temperature_C), y=factor)) + geom_jitter()
dev.off()
png(paste("imgs/", "fact3_barplot_tempG.png", sep=""))
ggplot(df.fac3.sum, aes(x=as.factor(Temperature_C), y=factor)) + geom_bar(stat="identity", color="black", fill="white") + theme_minimal() + xlab("Temperature C") + ylab("ILR Transformed Read Counts") + geom_errorbar(aes(ymin=factor-sd, ymax=factor+sd), width=.2) + geom_point(df.fac3.sum, mapping=aes(x=as.factor(Temperature_C), y=factor)) + geom_jitter()
dev.off()
```

![fact1 ILR counts temp](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/fact1_barplot_tempG.png)
![fact2 ILR counts temp](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/fact2_barplot_tempG.png)
![fact3 ILR counts temp](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/fact3_barplot_tempG.png)

### Random forest

```R
otu_table <- read.table("../01-raw_data_processing/sequence_table.16s.filtered.txt", sep="\t", header=T, row.names=1, stringsAsFactors=F, comment.char="")
otu_table <- t(otu_table)
metadata <- read.table("map.txt", sep="\t", header=T, row.names=1, stringsAsFactors=T, comment.char="")
metadata <- metadata[metadata$Season %in% c("winter", "summer"),]
metadata$Season <- factor(metadata$Season)
otu_nonzero_counts <- apply(otu_table, 1, function(y) sum(length(which(y > 0))))

remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

otu_table_rare_removed <- remove_rare(table=otu_table, cutoff_pro=0.1)
otu_table_rare_removed_norm <- sweep(otu_table_rare_removed, 2, colSums(otu_table_rare_removed), '/')*100
otu_table_scaled <- scale(otu_table_rare_removed_norm, center=T, scale=T)
otu_table_scaled_var <- data.frame(t(otu_table_scaled))
otu_table_scaled_var$var <- metadata[rownames(otu_table_scaled_var), "Season"]
set.seed(151)
otu_table_scaled_var <- otu_table_scaled_var %>% filter(!is.na(var))
rf_season <- randomForest(x=otu_table_scaled_var[,1:(ncol(otu_table_scaled_var)-1)], y=otu_table_scaled_var[, ncol(otu_table_scaled_var)], ntree=10000, importance=T, proximity=T)
rf_season
```

```text
Call:
 randomForest(x = otu_table_scaled_var[, 1:(ncol(otu_table_scaled_var) -      1)], y = otu_table_scaled_var[, ncol(otu_table_scaled_var)],      ntree = 10000, importance = T, proximity = T)
               Type of random forest: classification
                     Number of trees: 10000
No. of variables tried at each split: 14

        OOB estimate of  error rate: 6.9%
Confusion matrix:
       summer winter class.error
summer     32      2  0.05882353
winter      2     22  0.08333333
```

Lysing matrix

```R
metadata <- metadata[metadata$Matrix %in% c("E", "A"),]
metadata$Matrix <- factor(metadata$Matrix)
otu_table_scaled_var <- data.frame(t(otu_table_scaled))
otu_table_scaled_var$var <- metadata[rownames(otu_table_scaled_var), "Matrix"]
set.seed(151)
otu_table_scaled_var <- otu_table_scaled_var %>% filter(!is.na(var))
rf_matrix <- randomForest(x=otu_table_scaled_var[,1:(ncol(otu_table_scaled_var)-1)], y=otu_table_scaled_var[, ncol(otu_table_scaled_var)], ntree=10000, importance=T, proximity=T)
rf_matrix
```

```text
Call:
 randomForest(x = otu_table_scaled_var[, 1:(ncol(otu_table_scaled_var) -      1)], y = otu_table_scaled_var[, ncol(otu_table_scaled_var)],      ntree = 10000, importance = T, proximity = T)
               Type of random forest: classification
                     Number of trees: 10000
No. of variables tried at each split: 14

        OOB estimate of  error rate: 65.52%
Confusion matrix:
   A  E class.error
A 10 19   0.6551724
E 19 10   0.6551724
```

Insects

```R
metadata <- metadata[metadata$Insects %in% c("y", "n"),]
metadata$Insects <- factor(metadata$Insects)
otu_table_scaled_var <- data.frame(t(otu_table_scaled))
otu_table_scaled_var$var <- metadata[rownames(otu_table_scaled_var), "Insects"]
set.seed(151)
otu_table_scaled_var <- otu_table_scaled_var %>% filter(!is.na(var))
rf_insects <- randomForest(x=otu_table_scaled_var[,1:(ncol(otu_table_scaled_var)-1)], y=otu_table_scaled_var[, ncol(otu_table_scaled_var)], ntree=10000, importance=T, proximity=T)
rf_insects
```

```text
Call:
 randomForest(x = otu_table_scaled_var[, 1:(ncol(otu_table_scaled_var) -      1)], y = otu_table_scaled_var[, ncol(otu_table_scaled_var)],      ntree = 10000, importance = T, proximity = T)
               Type of random forest: classification
                     Number of trees: 10000
No. of variables tried at each split: 14

        OOB estimate of  error rate: 6.25%
Confusion matrix:
   n  y class.error
n 21  3       0.125
y  0 24       0.000
```

Temp group

```R
metadata <- metadata[metadata$Temp_group %in% c("20C", "30C", "40C", "50C"),]
metadata$Temp_group <- factor(metadata$Temp_group)
otu_table_scaled_var <- data.frame(t(otu_table_scaled))
otu_table_scaled_var$var <- metadata[rownames(otu_table_scaled_var), "Temp_group"]
set.seed(151)
otu_table_scaled_var <- otu_table_scaled_var %>% filter(!is.na(var))
rf_tgroup <- randomForest(x=otu_table_scaled_var[,1:(ncol(otu_table_scaled_var)-1)], y=otu_table_scaled_var[, ncol(otu_table_scaled_var)], ntree=10000, importance=T, proximity=T)
rf_tgroup
```

```text
Call:
 randomForest(x = otu_table_scaled_var[, 1:(ncol(otu_table_scaled_var) -      1)], y = otu_table_scaled_var[, ncol(otu_table_scaled_var)],      ntree = 10000, importance = T, proximity = T)
               Type of random forest: classification
                     Number of trees: 10000
No. of variables tried at each split: 14

        OOB estimate of  error rate: 21.74%
Confusion matrix:
    20C 30C 40C class.error
20C   1   6   1  0.87500000
30C   0  26   2  0.07142857
40C   0   1   9  0.10000000
```

What top 5 taxa are most important in the random forest models?

```R
png(paste("imgs/", "rf_season_importance.png", sep=""))
varImpPlot(rf_season)
dev.off()
pdf(paste("imgs/", "rf_season_importance.pdf", sep=""))
varImpPlot(rf_season)
dev.off()
```

![rf_season](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/rf_season_importance.png)

Season:
1. ASV6	Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium;Clostridium_unknown
2. ASV10	Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium;Clostridium_unknown
3. ASV2	Bacteria;Firmicutes;Bacilli;Bacillales;Planococcaceae;Sporosarcina;Sporosarcina_unknown
4. ASV4	Bacteria;Actinobacteria;Actinobacteria_c;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium_unknown
5. ASV25	Bacteria;Proteobacteria;Gammaproteobacteria;Ignatzschineria_o;Ignatzschineria_f;Ignatzschineria_f_unknown;Ignatzschineria_f_unknown

```R
png(paste("imgs/", "rf_insects_importance.png", sep=""))
varImpPlot(rf_insects)
dev.off()
pdf(paste("imgs/", "rf_insects_importance.pdf", sep=""))
varImpPlot(rf_insects)
dev.off()
```

![rf insects](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/rf_insects_importance.png)

Insects:
1. ASV3	Bacteria;Proteobacteria;Gammaproteobacteria;Ignatzschineria_o;Ignatzschineria_f;Ignatzschineria_f_unknown;Ignatzschineria_f_unknown
2. ASV1	Bacteria;Proteobacteria;Gammaproteobacteria;Ignatzschineria_o;Ignatzschineria_f;Ignatzschineria_f_unknown;Ignatzschineria_f_unknown
3. ASV73	Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus_unknown
4. ASV89	Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Morganellaceae;Providencia;Providencia_unknown
5. ASV53	Bacteria;Firmicutes;Bacilli;Lactobacillales;Carnobacteriaceae;Trichococcus;Trichococcus_unknown

Plot for figure

```R
impToPlot.season <- importance(rf_season, scale=F)[,3]
impToPlot.season <- sort(impToPlot.season)
short.imp <- tail(impToPlot.season, 10)
pdf(paste("imgs/", "rf_season_importance.filt.pdf", sep=""))
dotchart(short.imp, xlim=c(0.00, 0.06))
dev.off()
impToPlot.season <- importance(rf_season, scale=F)[,4]
impToPlot.season <- sort(impToPlot.season)
short.imp <- tail(impToPlot.season, 10)
pdf(paste("imgs/", "rf_season_importance.filt.gini.pdf", sep=""))
dotchart(short.imp, xlim=c(0.0, 2.5))
dev.off()
impToPlot.insect <- importance(rf_insects, scale=F)[,3]
impToPlot.insect <- sort(impToPlot.insect)
short.imp <- tail(impToPlot.insect, 10)
pdf(paste("imgs/", "rf_insects_importance.filt.pdf", sep=""))
dotchart(short.imp, xlim=c(0.00, 0.06))
dev.off()
impToPlot.insect <- importance(rf_insects, scale=F)[,4]
impToPlot.insect <- sort(impToPlot.insect)
short.imp <- tail(impToPlot.insect, 10)
pdf(paste("imgs/", "rf_insects_importance.filt.gini.pdf", sep=""))
dotchart(short.imp, xlim=c(0.0, 2.5))
dev.off()
```

### Correlation between temp/days in field and ratio between major phyla 

added a pseudo count of 1 to get around divide by zero issues

```R
png(paste("imgs/", "corr_actino-firmi_temp.png", sep=""))
qplot(as.numeric(as.character(seqtab.phylum$Temperature_C)), seqtab.phylum$actino.firmi, geom=c("point", "smooth")) + theme_minimal()
dev.off()
png(paste("imgs/", "corr_actino-proteo_temp.png", sep=""))
qplot(as.numeric(as.character(seqtab.phylum$Temperature_C)), seqtab.phylum$actino.proteo, geom=c("point", "smooth")) + theme_minimal()
dev.off()
png(paste("imgs/", "corr_firmi-proteo_temp.png", sep=""))
qplot(as.numeric(as.character(seqtab.phylum$Temperature_C)), seqtab.phylum$firmi.proteo, geom=c("point", "smooth")) + theme_minimal()
dev.off()
```

![corr_temp1](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/corr_actino-firmi_temp.png)
![corr_temp2](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/corr_actino-proteo_temp.png)
![corr_temp3](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/corr_firmi-proteo_temp.png)


```R
png(paste("imgs/", "corr_actino-firmi_days.png", sep=""))
qplot(seqtab.phylum$Days_in_Field, seqtab.phylum$actino.firmi, geom=c("point", "smooth")) + theme_minimal()
dev.off()
png(paste("imgs/", "corr_actino-proteo_days.png", sep=""))
qplot(seqtab.phylum$Days_in_Field, seqtab.phylum$actino.proteo, geom=c("point", "smooth")) + theme_minimal()
dev.off()
png(paste("imgs/", "corr_firmi-proteo_days.png", sep=""))
qplot(seqtab.phylum$Days_in_Field, seqtab.phylum$firmi.proteo, geom=c("point", "smooth")) + theme_minimal()
dev.off()
```

![corr_days1](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/corr_actino-firmi_days.png)
![corr_days2](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/corr_actino-proteo_days.png)
![corr_days3](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/corr_firmi-proteo_days.png)
