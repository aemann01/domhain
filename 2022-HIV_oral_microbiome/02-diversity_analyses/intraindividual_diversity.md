### Intraindividaul diversity

We want to see if there are patterns of microbial diversity in the same individual (same tooth type) that are indicative of HIV status. Three individuals with multiple sampling at the same time point with the same tooth type (CA-PD), one HEU (DM00021), one HUU (DM00036), and one HI (DM00024)

```R
library(phyloseq)
library(ggplot2)
# First need to filter out our samples from the larger phyloseq object
sub.ps.dat <- subset_samples(ps.dat.noUS, study_id == "DM00021" | study_id == "DM00036" | study_id == "DM00024")
sub.ps.dat <- subset_samples(sub.ps.dat, aliquot_type == "CA-PD")
sample_data(sub.ps.dat)['manifest_id'] <- row.names(sample_data(sub.ps.dat))

# ordinate and plot
ordcap <- ordinate(sub.ps.dat, "CAP", "bray", ~study_group)
pdf("img/intraindividual_div.capscale_plt.pdf")
plot_ordination(sub.ps.dat, ordcap, "samples", color="study_group", label="manifest_id") + theme_minimal()
dev.off()

# taxonomy plot
# first need to get relative abundance normalized data
rel.abund <- transform_sample_counts(sub.ps.dat, function(x) x/sum(x))
glom <- tax_glom(rel.abund, taxrank=rank_names(rel.abund)[8]) # collapse 
data <- psmelt(glom) # create dataframe from phyloseq object
data$Sample <- factor(data$Sample, levels=unique(data$Sample))
pdf("img/intraindividual_div.taxbar.pdf")
ggplot(data, aes(x=Sample, y=Abundance, fill=V4)) + geom_bar(aes(), stat="identity", position="stack") + theme_minimal() + theme(axis.text.x = element_text(angle = 90))
dev.off()

# streptococcus only 
strep.dat <- subset_taxa(rel.abund, V8 == "Streptococcus")
glom <- tax_glom(strep.dat, taxrank=rank_names(strep.dat)[9]) # collapse 
data <- psmelt(glom) # create dataframe from phyloseq object
data$Sample <- factor(data$Sample, levels=unique(data$Sample))
pdf("img/intraindividual_div.streptaxbar.pdf")
ggplot(data, aes(x=Sample, y=Abundance, fill=V9)) + geom_bar(aes(), stat="identity", position="stack") + theme_minimal() + theme(axis.text.x = element_text(angle = 90))
dev.off()
```






