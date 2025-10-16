## Imports
import os
import glob
import argparse
import pandas as pd

## Parser
parser = argparse.ArgumentParser(prog='Summarize18Sreads.py', description='Python script to summarize the 18S reads per group and threshold.')
parser.add_argument('-i', '--infolder', dest='infolder', action='store', help='18S ribosomal reads results folder.', required=True)
parser.add_argument('-b', '--blast', dest='blast', action='store', help='18S ribosomal accession Blast results folder (threshold: 70).', required=True)
parser.add_argument('-p', '--pident', dest='pident', action='store', help='Percentage of identity to filter blast results.', required=True, type=int)
parser.add_argument('-o', '--outfolder', dest='outfolder', action='store', help='Output folder to save results.', required=True)

## Functions

# List files
def list_files(infolder):
	# List files
	lst = sorted(glob.glob(f"{infolder}/*"))
	# Return
	return lst

# Load file
def load_file(infile: str, load_header: bool):
	# Load file
	if load_header:
		df = pd.read_table(infile)
	else:
		df = pd.read_table(infile, header=None)
	# Return dataset
	return df

# Extract reads
def extract_reads(infolder, pident):
	# Define files
	tsvs = sorted(glob.glob(f"{infolder}/*.other_euka.tsv"))
	# Store results
	res = []
	# Loop files
	for tsv in tsvs:
		try:
			# Extract reads
			reads = list(load_file(tsv, False)[0].unique())
			res.append(reads)
		except:
			print(f"Warning: {tsv}")
	# Flat list
	res = [j for i in res for j in i]
	# Return
	return res

# Count reads
def count_reads_r18S(infolder):
	# Define files
	tsvs = sorted(glob.glob(f"{infolder}/*.r18S.tsv"))
	# Store results
	res = []
	# Loop files
	for tsv in tsvs:
		try:
			count = len(load_file(tsv, False)[0].unique())
			res.append(count)
		except:
			print(f"Warning: {tsv}")
	# Sum individual counts
	res = sum(res)
	# Return
	return res

# Save results
def save_results(df, outfile):
	df.to_csv(outfile, sep="\t", index=None)

## Exec
if __name__ == "__main__":
	# Define arguments
	args = parser.parse_args()
	infolder = args.infolder
	blast = args.blast
	pident = args.pident
	outfolder = args.outfolder

	# Define samples
	samples = list_files(infolder)

	# Loop samples
	for sample in samples:
		# Define sample name
		name = os.path.basename(sample)
		# Define threshold
		threshold = os.path.basename(outfolder)
		# Create output folder
		if not os.path.isdir(f"{outfolder}"): os.makedirs(f"{outfolder}")
		# Load blast results
		hits = load_file(f'{blast}/{name}_r18S.tsv', False)
		# Filter blast results
		hits = hits[hits[2] > pident]
		# Extract all reads
		reads = extract_reads(sample, pident)
		# Filter blast to match reads
		hits = hits[hits[0].isin(reads)]
		# Count reads
		acce_reads = count_reads_r18S(sample)
		blast_reads = len(hits[0].unique())
		# Create dataframe
		df = pd.DataFrame({'type': ['Accession','Blast'], 'count': [acce_reads, blast_reads], 'sample': name, 'threshold': threshold})
		# Save results
		save_results(df, f"{outfolder}/{name}.tsv")