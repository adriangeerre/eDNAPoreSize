## Imports
import os
import glob
import argparse
import pandas as pd

## Parser
parser = argparse.ArgumentParser(prog='CountUniqueTaxids.py', description='Python script to count the unique taxids per eDNA sample.')
parser.add_argument('-i', '--infolder', dest='infolder', action='store', help='Filtered blast results folder.', required=True)
parser.add_argument('-c', '--categories', dest='categories', action='store', help='Categories (.dmp) file from NCBI (taxcat.tar.gz)', required = True)
parser.add_argument('-o', '--outfolder', dest='outfolder', action='store', help='Output folder to save results.', required=True)

## Functions

# Load file
def load_file(infile):
	# Load file
	df = pd.read_table(infile)
	# Return dataset
	return df

## Exec
if __name__ == "__main__":
	# Define arguments
	args = parser.parse_args()
	infolder = args.infolder
	categories = args.categories
	outfolder = args.outfolder

	# Create output folder
	if not os.path.isdir(outfolder): os.makedirs(outfolder)

	# Read categories
	categories = pd.read_table(categories, header = None, names = ['domain','sp_taxid','taxid'])

	# Define files
	files = glob.glob(f'{infolder}/*.taxid_counts.tsv')

	# Loop files
	for f in files:
		# Load file
		df = load_file(f)

		# Merge files
		df = df.merge(categories, left_on="Taxid", right_on="taxid", how="left").drop(['sp_taxid','taxid'], axis=1)

		# Fill NAs
		df['domain'] = df.domain.fillna('M')

		# Save file
		name = os.path.basename(f)
		outfile = name.replace(".tsv",".cat.tsv")
		df.to_csv(f"{outfolder}/{outfile}", index=None, sep="\t")

		# Summarize
		df = df.loc[:,['Count','domain']].groupby('domain').sum('Count').reset_index()

		# Save file
		outfile = f'{name.split(".")[0]}.domains.tsv'
		df.to_csv(f"{outfolder}/{outfile}", index=None, sep="\t")