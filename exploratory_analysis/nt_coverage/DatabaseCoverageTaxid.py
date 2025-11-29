## Imports
import os
import argparse
from ete4 import NCBITaxa

## Parser
parser = argparse.ArgumentParser(prog='DatabaseCoverageTaxid.py', description='Python script to estimate coverage, as branching down of a taxid.')
parser.add_argument('-i','--infolder', dest="infolder", action='store', help='File extracted from blast database using: blastdbcmd.')
parser.add_argument('-t','--taxid_list', dest="taxid_list", action='store', help='List of taxids to estimate coverage.')
parser.add_argument('-c','--coverage_output', dest="coverage_output", action='store', help='Output file to save results.')
parser.add_argument('-o','--outdir', dest="outdir", action='store', help='Path to save output files.')

## Functions

# Load taxid list
def load_taxid_list(taxid_list):
	with open(taxid_list,"r") as f:
		lines = f.readlines()
		lst = [i.strip() for i in lines]
	return lst

# Define descendant taxids
def define_descendants(ncbi, taxid):
	# Define descendants
	descendants = ncbi.get_descendant_taxa(taxid, intermediate_nodes=True)
	print(f"Found {len(descendants)} total descendant under taxid {taxid}")
	# Return
	return descendants

# Gather existing accessions
def gather_existing_accessions(existing, taxid, outdir):
	# Create output directory
	if not os.path.isdir(outdir): os.makedirs(outdir)
	# Save into file
	with open(f'{outdir}/{taxid}.txt','w') as w:
		for p in existing:
			# Extract accessions
			with open(p, 'r') as r:
				lines = r.readlines()
			# Write results
			for a in lines:
				w.write(f'{p}\t{a.strip()}\n')

# Pipeline
def pipeline(infolder, ncbi, taxid_list, coverage_output, outdir):
	# Store results
	res = []
	# Loop taxids
	for taxid in taxid_list:
		# Define descendants
		descendants = define_descendants(ncbi, taxid)
		# Filter existing descendants
		existing = [f'{infolder}/{i}.tsv' for i in descendants if os.path.exists(f'{infolder}/{i}.tsv')]
		# Gather existing
		gather_existing_accessions(existing, taxid, outdir)
		# Count accessions
		nacces = sum([sum(1 for _ in open(i)) for i in existing])
		# Append results
		res.append([taxid, len(descendants),len(existing),nacces])
	# Save results
	save_coverage_results(res, coverage_output)
	# Return
	return res

# Save results
def save_coverage_results(res, coverage_output):
	# Create output directory
	outdir = os.path.basename(coverage_output)
	if not os.path.isdir(outdir): os.makedirs(outdir)
	# Save file
	with open(coverage_output,"w") as w:
		w.write(f'\ttaxid\ttotal_taxids\tdb_taxids\tno_of_accessions\n')
		for l in res:
			w.write(f'{l[0]}\t{l[1]}\t{l[2]}\t{l[3]}\n')

## Exec
if __name__ == "__main__":
	# Define arguments
	args = parser.parse_args()
	infolder = args.infolder
	taxid_list = args.taxid_list
	coverage_output = args.coverage_output
	outdir = args.outdir

	# Load class
	ncbi = NCBITaxa()

	# Read taxid
	taxid_list = load_taxid_list(taxid_list)

	# Estimate coverage
	pipeline(infolder, ncbi, taxid_list, coverage_output, outdir)
