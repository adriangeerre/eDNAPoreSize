## Imports
import os
import glob
import argparse
import pandas as pd

## Parser
parser = argparse.ArgumentParser(prog='CountPhylum18S.py', description='Python script to count the number of 18S reads per phylum.')
parser.add_argument('-i', '--infile', dest='infile', action='store', help='File containing the LCA results.', required=True)
parser.add_argument('-a', '--acces', dest='acces', action='store', help='Folder containing the 18S ribosomal reads obtained from Blast accession (initial results).', required=True)
parser.add_argument('-b', '--blast', dest='blast', action='store', help='18S ribosomal reads obtained from Blast 18S mapping results.', required=True)
# parser.add_argument('-p', '--pident', dest='pident', action='store', help='Percentage of identity to filter 18S blast results.', required=True, type=int)
parser.add_argument('-o', '--outfolder', dest='outfolder', action='store', help='Output folder to save results.', required=True)

### -i: 06-mergeLCA/A0_2.lca.tsv
### -a: 10-Estimate18S/90/A0_2
### -b: 11-Blast_r18S/blast_results/A0_2_r18S.tsv

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

# Group accession results
def merge_acces_results(lst):
	# Filter values
	lst = [i for i in lst if i.endswith(".r18S.tsv")]
	# Load files
	dfs = [load_file(i, True) for i in lst]
	# Merge datasets
	df = pd.concat(dfs, axis=0)
	# Return
	return df

# Filter and count reads
def count_r18S_reads_per_phylum(lca, acces, blast, sample):
	# Filter data
	one = lca[lca['V1'].isin(list(acces['0'].unique()))]
	two = lca[lca['V1'].isin(list(blast[0].unique()))]
	# Merge data and remove redundancy
	df = pd.concat([one, two], axis=0).reset_index().drop('index', axis=1).drop_duplicates()
	# Count reads per phylum
	df = df.groupby(['superkingdom','phylum'], dropna=False).count().reset_index().rename(columns = {'V1': 'r18S'}).fillna('Unknown')
	# Count total reads
	total = lca.groupby(['superkingdom', 'phylum'], dropna=False).count().reset_index().rename(columns = {'V1': 'total'}).fillna('Unknown')
	# Merge total with r18S phylum count
	df = total.merge(df, on=['superkingdom', 'phylum'], how='left').fillna(value=0).astype({'r18S': 'int', 'total': 'int'})
	df['non-r18S'] = df['total'] - df['r18S']
	df = df.loc[:,['superkingdom','phylum','r18S','non-r18S','total']]
	# Add sample name
	df['sample'] = sample
	# Return
	return df

# Save results
def save_results(df, outfile):
	df.to_csv(outfile, sep="\t", index=None)

## Exec
if __name__ == "__main__":
	# Define arguments
	args = parser.parse_args()
	infile = args.infile
	acces = args.acces
	blast = args.blast
	#pident = args.pident
	outfolder = args.outfolder

	# Define sample name
	sample = os.path.basename(infile).split(".")[0]
	
	# Load LCA results
	lca = load_file(infile, True).loc[:,['V1','superkingdom', 'phylum']]
	
	# Load 18S accession reads
	acces_files = list_files(acces)
	acces = merge_acces_results(acces_files)
	
	# Load 18S blast reads
	blast = load_file(blast, False)

	# Filter
	df = count_r18S_reads_per_phylum(lca, acces, blast, sample)

	# Save results
	if not os.path.exists(outfolder): os.makedirs(outfolder)
	save_results(df, f"{outfolder}/{sample}.phylum_r18S.tsv")
