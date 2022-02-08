# *rpoC* Oral Bacteria Database Build

### 1. Download bacterial assembly data from NCBI

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
```

### 2. Download all assemblies with coding sequences annotated

```bash
cat assembly_summary.txt | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/$/\/*cds_from_genomic.fna.gz/' > query
wget -i query &
# make sure you got them all
grep "pattern" wget-log -B 14 | grep "2021" | grep "ftp:" | awk '{print $3}' > missing
wget -i missing &
# remove duplicates
rm *fna.gz.1
# concatenate together
cat *fna.gz > all_genomes.fna.gz
gzip -d all_genomes.fna.gz
```

### 3. Pull rpoC annotated sequences

*NOTE: not pulling gamma subunits -- these are found in chloroplasts and are not the same as rpoC in bacteria. Will have 15690 rpoC loci*

```bash
grep "gene=rpoC" all_genomes.fna | grep -v "gamma" | awk '{print $1}' | sed 's/>//' > rpoc.ids
seqtk subseq all_genomes.fna rpoc.ids > rpoc.fna
rm all_genomes.fna
```

### 4. Generate custom kraken2 database from rpoC sequences

4a. Get each element that will be needed to generate the formatted fasta file

```bash
grep ">" rpoc.fna > headers
sed -i 's/|/_/' headers          
grep -v ">" rpoc.fna > seqs
awk -F"_" '{print $2}' headers > accessions
```

4b. Use accessions to get taxids for each entry (should have 15690)

```bash
cat accessions | while read line; do esummary -db nuccore -id $line | xtract -pattern DocumentSummary -element TaxId,Title; done > taxids 2>esummary.err
```

*NOTE: If there are connection issues, the esummary.err file will tell you which accessions you didn't successfully pull a taxid for*

4c. Generate master fasta file in kraken's format requirements

```bash
paste accessions taxids | sed 's/^/>/' | sed 's/\t/|kraken:taxid|/' | sed 's/\t/ /' > fixed_headers
paste fixed_headers seqs | sed 's/\t/\n/' > rpoc_ref.fa
```

4d. Create reference database

```bash
cd ~/refdb # change me to where you want to keep your database folder
mkdir kraken_rpoc
kraken2-build --download-taxonomy --db kraken_rpoc/
kraken2-build --add-to-library ~/domhain/00-database_build/rpoc_ref.fa --db kraken_rpoc
kraken2-build --build --max-db-size 8000000000 --db kraken_rpoc
```
