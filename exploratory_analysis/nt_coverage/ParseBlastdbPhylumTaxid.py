## Imports
import os
import argparse
from ete4 import NCBITaxa
from collections import Counter

## Parser
parser = argparse.ArgumentParser(prog='ParseBlastbPhylumTaxid.py', description='Python script to obtain accession branching down of a taxid.')
parser.add_argument('-b','--blast_file', dest="blast_file", action='store', help='File extracted from blast database using: blastdbcmd.')
parser.add_argument('-t','--taxid', dest="taxid", action='store', type=int, help='Taxid to branch down.')
parser.add_argument('-o','--outdir', dest="outdir", action='store', help='Path to save output files')

## Functions

# Define descendant taxids
def define_descendants(ncbi, taxid):
	# Define descendants
	descendants = ncbi.get_descendant_taxa(taxid, intermediate_nodes=True)
	print(f"Found {len(descendants)} total descendant under taxid {taxid}")
	# Return
	return descendants

# Report
def report_rank(ncbi, descendants):
	# Extract taxonomic levels
	ranks = ncbi.get_rank(descendants)
	# Return
	print(Counter(ranks.values()))

# Save descendants
def save_results(lst, outfile):
	# Open file
	with open(outfile, "w") as w:
		# Loop taxids
		for i in lst:
			# Write
			w.write(f'{i}\n')

# Parse blast database
def parse_blastdb(blast_file):
	# Store results
	res = {}
	# Open
	with open(blast_file, "r") as f:
		for line in f:
			acce, tax = line.strip().split(" ")
			if int(tax) not in res.keys():
				res[int(tax)] = [acce]
			else:
				res[int(tax)].append(acce)
	# Report
	print(f"Database contains a total of {len(res.keys())} taxids")
	# Return
	return res

# Extract accessions per taxid in descendants
def extract_acces_descendant(descendants, bdb):
	# Save unmatched taxids
	unmatched = []
	# Loop descendants
	for ds in descendants:
		# Extract accessions
		if ds in bdb.keys():
			acces = bdb[ds]
			if not os.path.isdir(f'{outdir}/acces'): os.makedirs(f'{outdir}/acces')
			save_results(acces, f'{outdir}/acces/{ds}.tsv')
		else:
			unmatched.append(ds)
	# Save unmatched taxids
	save_results(unmatched, f'{outdir}/unmatched_taxids.tsv')

## Exec
if __name__ == "__main__":
	# Define arguments
	args = parser.parse_args()
	blast_file = args.blast_file
	taxid = args.taxid
	outdir = args.outdir

	# Load class
	ncbi = NCBITaxa()

	# Define descendants
	descendants = define_descendants(ncbi, taxid)

	# Report ranks
	report_rank(ncbi, descendants)

	# Save descendants
	if not os.path.isdir(outdir): os.makedirs(outdir)
	save_results(descendants, f'{outdir}/descendats.{taxid}.tsv')

	# Parse blast database
	bdb = parse_blastdb(blast_file)

	# Parse descendant against database
	extract_acces_descendant(descendants, bdb)