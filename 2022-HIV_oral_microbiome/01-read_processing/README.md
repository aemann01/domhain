# Raw data processing

### 1. Setup

Download raw data from ////

```bash
cd $HOME/domhain/2022-HIV_oral_microbiome/02-read_processing 
mkdir raw && cd raw
wget ////
cd ..
```

### 2. Run quality filtering, read merging, ASV generation, chimera removal 

To run jupyter notebook for DADA2 processing:

```bash
jupyter-notebook
```

### 3. Taxonomic assignment with Kraken2

3a. Run kraken2

```bash
kraken2 --db ~/refdb/kraken_rpoc/ --threads 8 --use-names --output rep_set.kraken.out rep_set.fa --unclassified-out rep_set.unclassified.kraken.out --confidence 0.01
```

```text
22766 sequences (10.89 Mbp) processed in 0.485s (2816.9 Kseq/m, 1347.04 Mbp/m).
  22205 sequences classified (97.54%)
  561 sequences unclassified (2.46%)
```

3b. Get taxonomic lineage information from tax ids (uncomment lines below when running for the first time)

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

We also have to change taxid 81850 to 33958 (was merged on March 11, 2021). Open query file in nano and do CTL+W and replace string (only found once in file). 

*NOTE! Depending on the version of the ncbi taxonomy file you pull and kraken version, you may have some issues with taxon names not being in the file.* -->

3c. Use queries to pull full lineage information (takes a little bit of time)

```bash
cat query | while read line; do grep -w -m 1 ^$line ~/refdb/ncbi_taxonomy/fullnamelineage.dmp | awk -F"|" '{print $3, $2}' | sed 's/\t//g' | sed 's/  / /g' | sed 's/cellular organisms; //' | sed 's/; /;/g' | sed 's/ /_/g'; done > lineages
```

3d. Merge with ASV name to get taxonomy table

```bash
awk '{print $2}' rep_set.kraken.out > asvids
paste asvids lineages > taxonomy.txt
```

3e. Remove unwanted ASVs from taxonomy, sequence table

```bash
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" | awk '{print $1}' > wanted.ids
seqtk subseq rep_set.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" > taxonomy_bac.txt
python fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
sed -i 's/;/\t/g' taxonomy_bac.txt
```

3f. Build a tree from the filtered representative sequences

```bash
mafft --thread 55 rep_set.filt.fa > rep_set.align.fa
fasttree -nt rep_set.align.fa > rep_set.align.tre
```
