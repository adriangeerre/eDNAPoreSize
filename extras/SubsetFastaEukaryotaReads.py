## Imports
import os
import glob
import argparse
import pandas as pd
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed


## Parser
parser = argparse.ArgumentParser(prog='SubsetFastaEukaryotaReads.py', description='Python script to count the unique taxids per eDNA sample.')
parser.add_argument('-i', '--infolder', dest='infolder', action='store', help='Filtered blast results folder.', required=True)
parser.add_argument('-f', '--fasta', dest='fasta', action='store', help='Merged fasta file containing shotgun reads.', required = True)
parser.add_argument('-o', '--outfolder', dest='outfolder', action='store', help='Output folder to save results.', required=True)

## Functions

# Return entries in reads
def filter_entries(chunk, reads_set):
	# Return
	return [entry for entry in chunk if entry.id in reads_set]

# Parse fasta chunk
def parse_in_chunks(fasta, chunk_size):
	# Store results
	chunk = []
	# Parse fasta
	for num, entry in enumerate(SeqIO.parse(fasta, "fasta")):
		chunk.append(entry)
		# Report chunk
		if num % chunk_size == 0:
			print(f"Processed: {num} lines...")
			yield chunk
			chunk = []
	if chunk:
		yield chunk

# Parralel parse fasta
def parallel_filter(fasta, reads, workers, chunk_size):
	# Store results
	reads_set = set(reads)
	results = []
	# Process data
	with ProcessPoolExecutor(max_workers=workers) as executor:
		futures = [executor.submit(filter_entries, chunk, reads_set)
			for chunk in parse_in_chunks(fasta, chunk_size)]
		# Save chunks into single
		for i, future in enumerate(as_completed(futures), 1):
			res = future.result()
			results.extend(res)
	# Return
	return results

## Exec
if __name__ == '__main__':
	# Define arguments
	args = parser.parse_args()
	infolder = args.infolder
	fasta = args.fasta
	outfolder = args.outfolder

	# Define name
	name = os.path.basename(infolder)

	# Create output folder
	if not os.path.isdir(f"{outfolder}"): os.makedirs(f"{outfolder}")

	# Define files
	files = glob.glob(f"{infolder}/*.other_euka.tsv")

	# Load files
	dfs = [pd.read_table(i, header=None) for i in files]

	# Merge all files
	df = pd.concat(dfs)

	# Read names
	reads = df[0].unique()
	print(f"Number of reads to extract: {len(reads)}")

	# Delete unwanted variables
	del files; del dfs; del df

	# Parse fasta
	res = parallel_filter(fasta, reads, workers=4, chunk_size=1000000)


	# Save results
	outfile = f"{outfolder}/{name}.other_euka.fasta"
	SeqIO.write(res, outfile, "fasta")