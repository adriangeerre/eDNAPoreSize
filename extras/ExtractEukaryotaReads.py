## Imports
import os
import glob
import argparse
import pandas as pd

## Parser
parser = argparse.ArgumentParser(prog='ExtractEukaryotaReads.py', description='Python script to count the unique taxids per eDNA sample.')
parser.add_argument('-i', '--infolder', dest='infolder', action='store', help='Filtered blast results folder.', required=True)
parser.add_argument('-t', '--taxids', dest='taxids', action='store', help='Summary file containing taxids and counts.', required = True)
parser.add_argument('-r', '--r18S', dest='r18S', action='store', help='NCBI 18S accessions.', required = True)
#parser.add_argument('-f', '--fastq', dest='fastq', action='store', help='Fastq file containing shotgun reads.', required = True)
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

# Extract taxids
def load_extract_reads_taxids(infile: str, load_header: bool):
	# Load file
	df = load_file(infile, load_header)
	# Extract column
	df = df[[0,1,17]]
	# Return
	return df

## Exec
if __name__ == "__main__":
	# Define arguments
	args = parser.parse_args()
	infolder = args.infolder
	taxids = args.taxids
	r18S = args.r18S
	outfolder = args.outfolder

	# Read taxid file
	taxids = pd.read_table(taxids)

	# Filter taxid
	taxids = taxids[taxids['domain'] == "E"].Taxid.to_list()

	# 18S accessions
	r18S = pd.read_table(r18S, header=None)[0].to_list()

	# Define samples
	samples = list_files(infolder)

	# Loop samples
	for sample in samples:
		# Define sample name
		name = os.path.basename(sample)
		# Create output folder
		if not os.path.isdir(f"{outfolder}/{name}"): os.makedirs(f"{outfolder}/{name}")
		# List filtered blast results
		tsvs = list_files(sample)
		# Loop files
		for tsv in tsvs:
			# Extract taxids
			df = load_extract_reads_taxids(tsv, False)
			# Filter Eukaryota taxids
			df = df[df[17].isin(taxids)]
			# Divide data: 18S / U
			acce = df[df[1].isin(r18S)]
			unkn = df[~df[1].isin(r18S)]
			# Save 18S reads
			outfile = os.path.basename(tsv).replace(".tsv",".r18S.tsv")
			acce.to_csv(f"{outfolder}/{name}/{outfile}", index=None, header=None, sep="\t")
			# Save Unknown reads
			outfile = os.path.basename(tsv).replace(".tsv",".other_euka.tsv")
			unkn.to_csv(f"{outfolder}/{name}/{outfile}", index=None, header=None, sep="\t")
