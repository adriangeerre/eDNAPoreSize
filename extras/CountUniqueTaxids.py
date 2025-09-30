## Imports
import os
import glob
import argparse
import pandas as pd
from collections import Counter

## Parser
parser = argparse.ArgumentParser(prog='CountUniqueTaxids.py', description='Python script to count the unique taxids per eDNA sample.')
parser.add_argument('-i', '--infolder', dest='infolder', action='store', help='Filtered blast results folder.', required=True)
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
def load_extract_taxids(infile):
	# Load file
	df = load_file(infile)
	# Extract column
	taxids = df[17].to_list()
	# Return
	return taxids

# Counter to dataframe
def counter_to_pandas(d):
	# To pandas
	df = pd.DataFrame(d.items(), columns=['Taxid', 'Count'])
	# Order by count
	df = df.sort_values(by='Count', ascending=False).reset_index(drop=True)
	# Return
	return df

# Count taxids
def count_taxids(lsts):
	# Flatten list of lists
	flat = [j for i in lsts for j in i]
	# Count taxids
	counts = Counter(flat)
	# Counter to dataframe
	counts = counter_to_pandas(counts)
	# Return
	return counts

# Save results
def save_results(df, outfile):
	df.to_csv(outfile, sep="\t", index=None)

## Exec
if __name__ == "__main__":
	# Define arguments
	args = parser.parse_args()
	infolder = args.infolder
	outfolder = args.outfolder

	# List samples
	samples = list_files(infolder)

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
			taxids = load_extract_taxids(tsv)
			# Append 
			res.append(taxids)
			# Count taxids
			counts = count_taxids(res)
			# Save results
			if not os.path.exists(outfolder): os.makedirs(outfolder)
			save_results(counts, f"{outfolder}/{name}.taxid_counts.tsv")
