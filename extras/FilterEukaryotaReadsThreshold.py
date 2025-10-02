## Imports
import os
import glob
import argparse
import pandas as pd

## Parser
parser = argparse.ArgumentParser(prog='FilterEukaryotaReadsThreshold.py', description='Python script to define 18S reads at different thresholds.')
parser.add_argument('-i', '--infolder', dest='infolder', action='store', help='Filtered blast results folder.', required=True)
parser.add_argument('-r', '--r18Sfolder', dest='r18Sfolder', action='store', help='18S ribosomal reads lists folder.', required=True)
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
	# Return dataset
	return df

# Save results
def save_results(df, outfile):
	df.to_csv(outfile, sep="\t", index=None)

## Exec
if __name__ == "__main__":
	# Define arguments
	args = parser.parse_args()
	infolder = args.infolder
	r18Sfolder = args.r18Sfolder
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
		if not os.path.isdir(f"{outfolder}/{name}"): os.makedirs(f"{outfolder}/{name}")
		# List filtered blast results
		tsvs = list_files(sample)
		# Loop files
		for tsv in tsvs:
			# Define outfile name
			outname = "_".join(os.path.basename(tsv).split("_")[:4])
			# Load filtered reads
			filtered_reads = list(load_file(tsv)[0].unique())
			# Load 18S reads file
			try:
				# Load dataset
				r18S = load_file(f"{r18Sfolder}/{name}/{outname}_70.r18S.tsv")
				# Filter datasets
				r18S = r18S[r18S[0].isin(filtered_reads)]
				# Save results
				save_results(r18S, f"{outfolder}/{name}/{outname}_{threshold}.r18S.tsv")
			except:
				print(f"Warning: {outname}_70.r18S.tsv") # Normally, empty files
			# Load other reads file
			try:
				# Load dataset
				other_euka = load_file(f"{r18Sfolder}/{name}/{outname}_70.other_euka.tsv")
				# Filter dataset				
				other_euka = other_euka[other_euka[0].isin(filtered_reads)]
				# Save results
				save_results(other_euka, f"{outfolder}/{name}/{outname}_{threshold}.other_euka.tsv")
			except:
				print(f"Warning: {outname}_70.other_euka.tsv")


