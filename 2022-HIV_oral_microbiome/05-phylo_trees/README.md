# Treponema tree
```bash
grep "Treponema" ../01-read_processing/taxonomy_bac.filt.txt | awk '{print $1}' | sed 's/"//g' > treponema.ids
seqtk subseq ../01-read_processing/rep_set.filt.fa treponema.ids > treponema.fa
mafft --auto treponema.fa > treponema.align.fa
raxmlHPC-PTHREADS-SSE3 -T 25 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n treponema.tre -s treponema.align.fa
```

# get data formatted for R

```bash
python3 summarize_seqtab.py -c study_group
python3 summarize_seqtab.py -c aliquot_type
````

# R tree build
# clean up and root tree with figtree before continuing

```R
library(phytools)
library(ape)
tree <- read.tree('treponema.root.tre')
dat <- read.table("study_group_abund.txt", sep="\t", row.names=1, header=T)
col_order <- c("HI", "HEU", "HUU")
dat <- dat[,col_order]
# filt <- na.omit(dat[tree$tip.label,])

pdf("trep_barplot_study_group.pdf")
plotTree.barplot(tree, dat, args.barplot=list(legend.text=TRUE,space=c(0,1.2),args.legend=list(x=1,y=17),col=c("#8214a0ff","#fa78faff","#00a0faff")))
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
par(fg="black")
for(i in 1:Ntip(tree)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(tree)])),
    rep(obj$yy[i],2),lty="dotted")
dev.off()

dat <- read.table("aliquot_type_abund.txt", sep="\t", row.names=1, header=T)

pdf("trep_barplot_aliquot_type.pdf")
plotTree.barplot(tree, dat, args.barplot=list(legend.text=TRUE,space=c(0,1.2),args.legend=list(x=1,y=17),col=c("#aa0a3cff","#fa7850ff","#f0f032ff", "#005ac8ff", "#14d2dcff", "#0ab45aff")))
dev.off()
```

























# Streptococcus tree

grep "Streptococcus" ../01-read_processing/taxonomy_bac.filt.txt | awk '{print $1}' | sed 's/"//g' > strep.ids
seqtk subseq ../01-read_processing/rep_set.filt.fa strep.ids > strep.fa
mafft --auto strep.fa > strep.align.fa
raxmlHPC-PTHREADS-SSE3 -T 25 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n strep.tre -s strep.align.fa




# Streptococcus tree

grep "Streptococcus" ../01-read_processing/taxonomy_bac.filt.txt | awk '{print $1}' | sed 's/"//g' > strep.ids
seqtk subseq ../01-read_processing/rep_set.filt.fa strep.ids > strep.fa
mafft --auto strep.fa > strep.align.fa
raxmlHPC-PTHREADS-SSE3 -T 25 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n strep.tre -s strep.align.fa
