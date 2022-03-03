#!/usr/bin/env python3

'''
///////
'''

## TO DO: ONLY INCLUDE DOMHAIN SAMPLES


import pandas as pd 
import argparse
from Bio import Phylo

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--category', help="category in mapping file to summarize data over")
args = parser.parse_args()

# load data
metadat = pd.read_csv("map.txt", sep="\t")
biom = pd.read_csv("../01-read_processing/sequence_table.merged.txt", sep="\t")
# format data
merged_biom = pd.merge(biom, metadat, left_on="row_names", right_on="manifest_id")
merged_biom["row_names"] = merged_biom[args.category] # replace names with category
clean_biom = merged_biom.iloc[:,0:biom.shape[1]] # now only select columns you want

# group by category and get mean read count per group
grouped_biom = clean_biom.groupby("row_names").mean() # group by category, mean across groups
grouped_biom = grouped_biom.transpose() # transpose
normdf = grouped_biom.div(grouped_biom.sum(axis=1), axis=0) # get relative abundance

with open("%s_abund.txt" % args.category, "w") as outfile:
	normdf.to_csv(outfile, sep="\t")

