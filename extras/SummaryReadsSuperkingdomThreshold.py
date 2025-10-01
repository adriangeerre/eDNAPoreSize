## Imports
import os
import glob
import argparse
import pandas as pd

## Parser
parser = argparse.ArgumentParser(prog='SummaryReadsSuperkingdomThreshold.py', description='Python script to count the unique taxids per eDNA sample.')
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
	df = pd.read_table(infile).set_index('id')
	# Return dataset
	return df

# Summarize data
def count_reads(tsv, name, val):
	# Load data
	df = load_file(tsv)
	# Divide data
	spks = df.loc[:,['A','B','E','V']]
	# Clean data
	single = spks[spks.sum(axis=1) == 1]
	multiple = spks[spks.sum(axis=1) > 1]
	unknown = spks[spks.sum(axis=1) == 0]
	# Summarize
	counts = single.sum(axis=0).reset_index().rename(columns = {'index': 'superkingdom', 0: 'read_counts'})
	# Add multiple and unknown
	new_rows = pd.DataFrame({'superkingdom': ['M','U'], 'read_counts': [multiple.shape[0], unknown.shape[0]]})
	counts = pd.concat([counts, new_rows], axis=0)
	# Add sample and value
	counts['sample'] = name
	counts['threshold'] = val
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

	# List thresholds
	folders = list_files(infolder)

	# Loop thresholds
	for folder in folders:
		# Define value name
		val=os.path.basename(folder)
		# List read files
		tsvs = list_files(folder)
		# Store results
		res = []
		# Loop files
		for tsv in tsvs:
			# Define sample
			name = os.path.basename(tsv).split(".")[0]
			# Count reads per category
			counts = count_reads(tsv, name, val)
			# Append
			res.append(counts)
		# Merge res
		df = pd.concat(res, axis=0)
		# Save results
		if not os.path.exists(outfolder): os.makedirs(outfolder)
		save_results(df, f"{outfolder}/summary.{val}.tsv")
