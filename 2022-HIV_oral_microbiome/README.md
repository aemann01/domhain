# HIV and the oral microbiome

This repository describes the analysis performed in the paper [*Mann et al. 2023. The impact of HIV on the oral microbiome of children living in sub-Saharan Africa using a rpoC gene-fragment metataxonomic approach. Microbiology Spectrum*](https://journals.asm.org/doi/10.1128/spectrum.00871-23). Associated data for this project is deposited in the European Nucleotide Archive under accession number [PRJEB60354](https://www.ebi.ac.uk/ena/browser/view/PRJEB60354)

[![DOI](https://zenodo.org/badge/457036815.svg)](https://zenodo.org/badge/latestdoi/457036815)

## Setup

This repository assumes you are running in a Unix environment (e.g., Mac OSX or Linux) and you have conda installed.

To get this repository:

- Install and set up anaconda or miniconda as described at the [bioconda
  documentation](https://bioconda.github.io/user/install.html), including
  setting up channels.
- Clone this repository to your machine and change into the directory with

```bash
git clone https://github.com/aemann01/domhain.git && cd domhain/2022-HIV_oral_microbiome
```

- Run the following command to install the environment

```bash
conda env create -f environment.yml
```

- To load the environment

```bash
conda activate 2022-HIV_oral_microbiome
```

- Add the R kernel to Jupyter by installing a kernel spec

```bash
R -e 'IRkernel::installspec()'
```

- To turn off the environment run

```bash
conda deactivate
```

## Repository structure

```
.
├── 00-database_build
│   ├── README.md
├── 01-read_processing
│   ├── DADA2_processing.Rmd
│   ├── DADA2_processing.ipynb
│   ├── README.md
│   └── fix_taxonomy.py
├── 02-diversity_analyses
│   ├── README.md
│   ├── diversity_analyses.Rmd
│   └── diversity_analyses.ipynb
├── 03-network_analyses
│   ├── Network_comparison.ipynb
│   ├── README.md
│   ├── network_analysis.Rmd
│   └── network_analysis.ipynb
├── 04-differential_abundance
│   ├── README.md
|   └── differential_abundance.ipynb
├── 05-phylo_trees
│   ├── README.md
|   └── summarize_seqtab.py
├── 06-pathogens
|   └── scratch.r
├── 07-age
|   └── age_analyses.ipynb
├── 08-low_smutans
│   ├── intraindividual_diversity.ipynb
|   └── low_smutans_capplot.ipynb
├── 09-ASV_analyses
|   └── ASV_analyses.ipynb
├── 10-siblings
│   ├── sibling_analyses.ipynb
|   └── siblings.txt
├── 12-revisions
|   └── revision_figures.ipynb
├── README.md
├── domhain_metadata.sql
└── environment.yml
```
