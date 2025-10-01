## Imports
import os
import glob
import argparse
import pandas as pd
from collections import Counter

## Parser
parser = argparse.ArgumentParser(prog='CountUniqueTaxids.py', description='Python script to count the unique taxids per eDNA sample.')
parser.add_argument('-i', '--infolder', dest='infolder', action='store', help='Filtered blast results folder.', required=True)
parser.add_argument('-c', '--categories', dest='categories', action='store', help='Categories (.dmp) file from NCBI (taxcat.tar.gz).', required = True)
parser.add_argument('-o', '--outfolder', dest='outfolder', action='store', help='Output folder to save results.', required=True)

## Functions

# List files
def list_files(infolder):
	# List files
	lst = sorted(glob.glob(f"{infolder}/*"))
	# Return
	return lst

# Load file
def load_file(infile):
	# Load file
	df = pd.read_table(infile, header=None)
	# Return dataset
	return df

# Extract taxids
def load_extract_reads(infile, categories):
	# Load file
	df = load_file(infile)
	# Reduce hits to reads
	df = df[[0,17]].drop_duplicates().rename(columns = {0: 'id',17: 'taxid'})
	# Add categories
	df = df.merge(categories, on="taxid", how="left").drop(['sp_taxid','taxid'], axis=1).fillna("M")
	# Filter duplicates categories per read
	df = df[['id','domain']].drop_duplicates()
	# Generate matrix
	df['value'] = 1
	df = df.pivot_table(columns='domain', index='id', values='value', fill_value=0)
	# Add missing columns
	cat_cols = sorted(list(categories.domain.unique())) + ["M"]
	for col in cat_cols:
		if col not in list(df.columns):
			df[col] = 0
	# Sort columns
	df = df.loc[:,cat_cols]
	# Reset index
	df = df.reset_index()
	# Return
	return df

# Save results
def save_results(df, outfile):
	df.to_csv(outfile, sep="\t", index=None)

## Exec
if __name__ == "__main__":
	# Define arguments
	args = parser.parse_args()
	infolder = args.infolder
	categories = args.categories
	outfolder = args.outfolder

	# List samples
	samples = list_files(infolder)

	# Read categories
	categories = pd.read_table(categories, header = None, names = ['domain','sp_taxid','taxid'])

	# Loop samples
	for sample in samples:
		# Define sample name
		name = os.path.basename(sample)
		# List blast results
		tsvs = list_files(sample)
		# Store results
		res = []
		# Loop files
		for tsv in tsvs:
			# Extract taxids
			cats = load_extract_reads(tsv, categories)
			# Append 
			res.append(cats)
		# Merge res
		matrix = pd.concat(res)
		# Save results
		if not os.path.exists(outfolder): os.makedirs(outfolder)
		save_results(matrix, f"{outfolder}/{name}.reads_cats.tsv")
