library(phyloseq)
library(ggplot2)
load("../02-diversity_analyses/.RData")

# function to summarize data
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

# red complex bacteria
ps.dat.core <- microbiome::transform(ps.dat, "compositional")
redcomp <- subset_taxa(ps.dat.core, V10 == "Tannerella_forsythia" | V10 == "Porphyromonas_gingivalis" | V8 == "Treponema_denticola")
pdf("redcomplex_abund.pdf")
plot_bar(redcomp, "V12", facet_grid=study_group~aliquot_type, fill="study_group") + geom_bar(aes(color=study_group, fill=study_group), stat="identity", position="stack")
dev.off()

df <- psmelt(redcomp)
df <- data_summary(df, varname="Abundance", groupnames=c("study_group", "V2"))

pdf("redcomp_boxplot.pdf")
ggplot(df, aes(x=study_group, y=log10(Abundance), group=study_group, color=study_group)) + geom_boxplot(fill="light grey") + geom_violin(alpha=0.25) + geom_jitter(position=position_jitter(width=0.1)) + scale_color_manual(values=c("#fa78faff", "#8214a0ff", "#00a0faff", "#fae6beff"))
dev.off()

wilcox.test(df[df$study_group == "HI",]$Abundance, df[df$study_group == "HEU",]$Abundance)

# 	Wilcoxon rank sum test with continuity correction

# data:  df[df$study_group == "HI", ]$Abundance and df[df$study_group == "HEU", ]$Abundance
# W = 5566395, p-value = 0.9508
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(df[df$study_group == "HI",]$Abundance, df[df$study_group == "HUU",]$Abundance)

# 	Wilcoxon rank sum test with continuity correction

# data:  df[df$study_group == "HI", ]$Abundance and df[df$study_group == "HUU", ]$Abundance
# W = 5758941, p-value = 0.06721
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(df[df$study_group == "HEU",]$Abundance, df[df$study_group == "HUU",]$Abundance)

# 	Wilcoxon rank sum test with continuity correction

# data:  df[df$study_group == "HEU", ]$Abundance and df[df$study_group == "HUU", ]$Abundance
# W = 4386652, p-value = 0.09359
# alternative hypothesis: true location shift is not equal to 0

# orange complex bacteria
orgcomp <- subset_taxa(ps.dat.core, V9 == "Campylobacter_gracilis" | V9 == "Campylobacter_rectus" | V9 == "Campylobacter_showae" | V8 == "Fusobacterium_nucleatum" | V8 == "Fusobacterium_periodonticum" | V10 == "Prevotella_intermedia" | V10 == "Prevotella_nigrescens" | V10 == "Streptococcus_constellatus")
pdf("orgcomplex_abund.pdf")
plot_bar(orgcomp, "V12", facet_grid=study_group~aliquot_type, fill="study_group") + geom_bar(aes(color=study_group, fill=study_group), stat="identity", position="stack")
dev.off()

df <- psmelt(orgcomp)
df <- data_summary(df, varname="Abundance", groupnames=c("study_group", "V2"))

pdf("redcomp_boxplot.pdf")
ggplot(df, aes(x=study_group, y=log10(Abundance), group=study_group, color=study_group)) + geom_boxplot(fill="light grey") + geom_violin(alpha=0.25) + geom_jitter(position=position_jitter(width=0.1)) + scale_color_manual(values=c("#fa78faff", "#8214a0ff", "#00a0faff", "#fae6beff"))
dev.off()

wilcox.test(df[df$study_group == "HI",]$Abundance, df[df$study_group == "HEU",]$Abundance)












