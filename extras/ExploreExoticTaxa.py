## Imports
import os
import re
import glob
import argparse
import pandas as pd

## Parser
parser = argparse.ArgumentParser(prog='ExploreExoticTaxa.py', description='Python script to explore the exotic taxa for shotgun sequencing.')
parser.add_argument('-i', '--infolder', dest='infolder', action='store', help='Filtered blast results folder.', required=True)
parser.add_argument('-t', '--tax_file', dest='tax_file', help='Tab-delimited table for the relation between taxid and taxonomy.')
parser.add_argument('--superkingdom', action='store', dest='s', help='Superkingdom to filter for taxids (multiple as comma-separated).', required=False, default=None, type=lambda x: re.split(r'[,\s]+', x.strip()))
parser.add_argument('--kingdom', action='store', dest='k', help='Kingdom to filter for taxids (multiple as comma-separated).', required=False, default=None, type=lambda x: re.split(r'[,\s]+', x.strip()))
parser.add_argument('--phylum', action='store', dest='p', help='Phylum to filter for taxids (multiple as comma-separated).', required=False, default=None, type=lambda x: re.split(r'[,\s]+', x.strip()))
parser.add_argument('--class', action='store', dest='c', help='Class to filter for taxids (multiple as comma-separated).', required=False, default=None, type=lambda x: re.split(r'[,\s]+', x.strip()))
parser.add_argument('--order', action='store', dest='o', help='Order to filter for taxids (multiple as comma-separated).', required=False, default=None, type=lambda x: re.split(r'[,\s]+', x.strip()))
parser.add_argument('--family', action='store', dest='f', help='Family to filter for taxids (multiple as comma-separated).', required=False, default=None, type=lambda x: re.split(r'[,\s]+', x.strip()))
parser.add_argument('--genus', action='store', dest='g', help='Genus to filter for taxids (multiple as comma-separated).', required=False, default=None, type=lambda x: re.split(r'[,\s]+', x.strip()))
parser.add_argument('-o', '--outfolder', dest='outfolder', action='store', help='Output folder to save results.', required=True)

## Functions

# List files
def list_files(infolder):
	# List files
	lst = sorted(glob.glob(f"{infolder}/*"))
	# Return
	return lst

# Load file
def load_file(infile: str):
	# Load file
	df = pd.read_table(infile, header=None)
	# Add header
	cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'slen', 'qlen', 'qcovs', 'qcovhsp', 'qcovus', 'taxid']
	df.columns = cols
	# Return dataset
	return df

# Extract taxids
def extract_taxids(tax_file, s, k, p, c, o, f, g):
	# Filter table
	if s:
		tax_file = tax_file[tax_file['superkingdom'].isin(s)]
	if k:
		tax_file = tax_file[tax_file['kingdom'].isin(k)]
	if p:
		tax_file = tax_file[tax_file['phylum'].isin(p)]
	if c:
		tax_file = tax_file[tax_file['class'].isin(c)]
	if o:
		tax_file = tax_file[tax_file['order'].isin(o)]
	if f:
		tax_file = tax_file[tax_file['family'].isin(f)]
	if g:
		tax_file = tax_file[tax_file['genus'].isin(g)]
	# Extract taxids
	taxids = list(set(tax_file['taxid'].to_list()))
	# Return
	return taxids

# Filter reads
def filter_results(tsv, taxids):
	# Load file
	df = load_file(tsv)
	# Select taxids
	df = df[df['taxid'].isin(taxids)]
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
	tax_file = args.tax_file
	s = args.s
	k = args.k
	p = args.p
	c = args.c
	o = args.o
	f = args.f
	g = args.g
	outfolder = args.outfolder

	# Read taxonomy file
	tax_file = pd.read_table(tax_file)

	# Extract taxids
	taxids = extract_taxids(tax_file, s, k, p, c, o, f, g)

	# Define samples
	samples = list_files(infolder)

	# Loop samples
	for sample in samples:
		# Report
		print(sample)
		# Define sample name
		name = os.path.basename(sample)
		# List filtered blast results
		tsvs = list_files(sample)
		# Loop and store results
		res = [filter_results(tsv, taxids) for tsv in tsvs]
		# To dataframe
		df = pd.concat(res)
		# Save results
		if not os.path.isdir(f"{outfolder}/{name}"): os.makedirs(f"{outfolder}/{name}")
		taxa = '_'.join([i[0] for i in [s,k,p,c,o,f,g] if not isinstance(i, type(None)) and len(i) == 1])
		save_results(df, f"{outfolder}/{name}/{taxa}.tsv")