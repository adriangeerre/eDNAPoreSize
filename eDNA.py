# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By    : Adrián Gómez Repollés
# Email         : adrian.gomez@mbg.au.dk
# Created Date  : 25/03/2022
# version       : '1.0'
# ---------------------------------------------------------------------------
""" Pipeline (gwf) to perform metagenomic analysis of paired-end eDNA illumi-
na reads. """ 
# ---------------------------------------------------------------------------
# Imports 
# ---------------------------------------------------------------------------

import os
import sys
from gwf import Workflow, AnonymousTarget

# Parameters
#-----------

account = "eDNA"
blastdb = "/project/eDNA/faststorage/blastdb/nt-old-30032020"
info_file = "samples.tsv" # Columns: Name Forward Reverse

# Workflow start
#---------------
gwf = Workflow(defaults={"account": account})

# Functions
#----------

# Merge paired-end reads
def merged_paired_reads(prefix, forward, reverse, merged, unmerged_forward, unmerged_reverse, max_mismatches, min_insert_size, overlap, memory):
    inputs = [forward, reverse]
    outputs = [merged, unmerged_forward, unmerged_reverse]
    options = {'cores': 4, 'memory': f'{memory}g', 'queue': 'short', 'walltime': '02:00:00'}

    spec = '''
    mkdir -p 01-Merge/{prefix}
	bbmerge.sh in1={forward} in2={reverse} out={merged} outu1={unmerged_forward} outu2={unmerged_reverse} maxmismatches={max_mismatches} mismatches={max_mismatches} mininsert={min_insert_size} minoverlap={overlap} > .logs/{prefix}.bbmerged.log
    '''.format(prefix=prefix, forward=forward, reverse=reverse, merged=merged, unmerged_forward=unmerged_forward, unmerged_reverse=unmerged_reverse, max_mismatches=max_mismatches, min_insert_size=min_insert_size, overlap=overlap, memory=memory)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Transform to fasta
def fastq2fasta(fqinput, faoutput, memory):
    inputs = [fqinput]
    outputs = [faoutput]
    options = {'cores': 4, 'memory': f'{memory}g', 'queue': 'short', 'walltime': '00:20:00'}

    spec = '''
    seqtk seq -a {fqinput} > {faoutput}
    '''.format(fqinput=fqinput, faoutput=faoutput)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Split the file
def split_file(prefix, num, in_route, out_route):
    inputs = [f'{in_route}/{prefix}.fasta']
    outputs = ["{}/{}_{}".format(out_route, prefix, "%04d" % i) for i in range(num)]
    options = {'cores': 1,'memory': '4g', 'queue': 'short', 'walltime': '00:30:00'}

    spec = '''
    mkdir -p {out_route}
    split -a 4 -d -l 10000 01-Merge/{prefix}/{prefix}.fasta '{prefix}_'
    mv {prefix}_* {out_route}
    '''.format(prefix = prefix, out_route=out_route, num = num)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Blast
def blast_file(fainput, in_route, out_route, dbs, blastdb, evalue, target_seqs, threads, memory):
    inputs = [f'{in_route}/{fainput}']
    outputs = [f'{out_route}/{fainput}.tsv']
    options = {'cores': 4,'memory': f'{memory}g', 'queue': 'normal', 'walltime': '24:00:00'}

    spec = '''
    mkdir -p {out_route}
    echo "Copying database $(date)"
    rsync -avr --exclude '*.tar.gz' --exclude '*.md5' --exclude 'seqs' {blastdb} /scratch/$SLURM_JOBID/
    echo "Database copied $(date)"
    blastn -query {in_route}/{fainput} -db /scratch/$SLURM_JOBID/nt-old-30032020/{dbs} -outfmt "6 std slen qlen sseq qseq qcovs qcovhsp qcovus staxids" -num_threads {threads} -evalue {evalue} -max_target_seqs {target_seqs} -out {out_route}/{fainput}.tsv
    '''.format(out_route=out_route, in_route=in_route, fainput=fainput, dbs=dbs, blastdb=blastdb, evalue=evalue, threads=threads, target_seqs=target_seqs)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Downstream analyses


# Targets
#--------
f = open(info_file,'r')

# Log folder
if not os.path.exists(".logs"): os.mkdir(".logs")

# Collect sample names
samples = []

for sample in f:
    
    if sample[0] == '#':
        continue

    # Prepare the data to be given to the functions
    sample = sample.strip().split("\t")
    forward = sample[0]	# forward 
    reverse = sample[1]	# reverse
    name = sample[2]	# name

    # Append sample names
    samples.append(name)

    # Execution
    #----------

    # Merge Paired-End Reads
    gwf.target_from_template(f'merge_{name}', merged_paired_reads(prefix=name, forward = forward, reverse = reverse, merged = f'01-Merge/{name}/{name}.merged.fastq', unmerged_forward=f'01-Merge/{name}/{name}_R1.unmerged.fastq', unmerged_reverse=f'01-Merge/{name}/{name}_R2.unmerged.fastq', max_mismatches = 0, min_insert_size = 100, overlap = 15, memory = 8))

    # Fastq to fasta
    gwf.target_from_template(f'fq2fa_{name}', fastq2fasta(fqinput = f'01-Merge/{name}/{name}.merged.fastq', faoutput = f'01-Merge/{name}/{name}.fasta', memory = 4))

    # Open file with number of splits
    lines_file = ".number_of_lines.tsv"
    d = {} # Dictionary number of lines
    if os.path.exists(lines_file):
        f = open(lines_file, "r")
        for row in f:
            (key,val) = row.split()
            d[key] = int(val)
        f.close()

	# Number of files
    div_lines = 5000 * 2
    if name not in d.keys():
	    if os.path.exists(f'01-Merge/{name}/{name}.fasta'):
	    	# Number of lines in fasta file
            num_lines = sum(1 for line in open(f'01-Merge/{name}/{name}.fasta'))
		    # Number of files created after split
            num_files = num_lines // div_lines
            if ((num_lines / div_lines) - num_files > 0):
		        num_files += 1
		    # Add number of lines
            f = open(lines_file, "a+")
            f.write(f'{name}\t{num_lines}\n')
            f.close()
    else:
		num_lines = d[name]
		# Number of files created after split
		num_files = num_lines // div_lines
		if ((num_lines / div_lines) - num_files > 0):
		    num_files += 1		

    # Split fasta file
    if "num_files" in locals():
    	gwf.target_from_template(f'split_{name}', split_file(prefix = name, num = num_files, in_route = f'01-Merge/{name}', out_route = f'02-Split/{name}'))
    
	    # Blast splitted files
        for j in range(num_files):
	        split_to_blast = "{}_{}".format(name, "%04d" % j)
            gwf.target_from_template(f'blast_{split_to_blast}', blast_file(fainput = split_to_blast, in_route = f'02-Split/{name}', out_route = f'03-Blast/{name}', dbs = 'nt', blastdb=blastdb, evalue=0.01, target_seqs=1000, threads=8, memory = 32))

    # Create Folder
    if not os.isdir(f"04-FilterBlast/{name}"): os.makedirs(f"04-FilterBlast/{name}")
    if not os.isdir("05-MergeBlast"): os.makedirs("05-MergeBlast")

    # Filter and Merge Blast results
    gwf.target(f'Filter_Blast_{name}', inputs=[f"03-Blast/{name}/{name}_{num}.tsv" for num in range(num_files)], outputs=[f"04-FilterBlast/{name}/{name}_filtered_{num}.tsv" for num in range(num_files)] + [f"05-MergeBlast/{name}.tsv"]) << """
        # Filter Blast
        for file in $(ls 03-Blast/{name}/ | cut -d '.' -f 1,2); do cut -f 1-14,17-18 $file | awk '($4 >= 100)' | awk '($3 >= 90)' | sort -k1,1 -k12,12nr -k11,11n > 04-FilterBlast/{name}/"$file"_filter.tsv; done

        # Merge Blast outputs
        cat 04-FilterBlast/{name}/*_filter.tsv > 05-MergeBlast/{name}.tsv
        """.format(name=name)

f.close()

# Install in ROBITools Environment 
# devtools::install_gitlab("obitools/ROBITaxonomy", host = "https://git.metabarcoding.org/")
# devtools::install_gitlab("obitools/ROBITools", host = "https://git.metabarcoding.org/")

# Obtain taxids
gwf.target('Obtain_taxids', inputs=[f"05-MergeBlast/{name}.tsv" for name in samples], outputs=["multiple.taxids.tsv","single.taxids.tsv"]) << """
# Obtain taxids
for file in $(ls 05-MergeBlast/)
do
	cut -f 16 $file | sort | uniq >> taxids.tmp.tsv
done
cat taxids.tmp.tsv | sort | uniq > taxids.tsv
rm taxids.tmp.tsv

# Divide taxids between single and multiple (semicolons separated)
grep ';' taxids.tsv > multiple.taxids.tsv
grep -v ';' taxids.tsv > single.taxids.tsv
"""

# MRCA for multiple taxids (TaxizeDB)
gwf.target('Clean_taxids_(MRCA)', inputs=["multiple.taxids.tsv"], outputs=["m2s.taxids.tsv"]) << """
Rscript mulipleTaxids2singleTaxid.lca.R
"""

# Get taxonomy (TaxizeDB)
gwf.target('Get_taxonomy', inputs=["single.taxids.tsv","m2s.taxids.tsv"], outputs=["taxid2taxonomy.tsv"]) << """
Rscript taxid2taxa.lca.R
"""

# Get MRCA from original Blast (TaxizeDB) -> Filtered for e-value < 0.0001, algn len >= 100 and pident >= 90
if not os.isdir("06-LCA"): os.makedirs("06-LCA")
gwf.target(f'Blast_MRCA_{name}', inputs=[f"05-MergeBlast/{name}.tsv" for name in samples] + ["taxid2taxonomy.tsv"], outputs=[f"06-LCA/{name}.lca.tsv" for name in samples]) << """
Rscript lca.R
"""

# Generate species-count matrix (TaxizeDB)
gwf.target('Count_matrix', inputs=[f"06-LCA/{name}.lca.tsv" for name in samples], outputs=["raw.taxa.count.tsv"]) << """
Rscript counts.lca.R
"""

# Rarefy (ROBITools)
gwf.target('Rarefy', inputs=["lca.taxa.count.tsv"], outputs=["Pore_size_table_for_rarefy.tsv","counts.lca.rarefy.tsv"]) << """
Rscript rarefyROBI.R
"""

# Obtain summary table and heatmap (Bioconductor)
gwf.target('Summary_table_and_heatmap', inputs=["counts.lca.rarefy.tsv"], outputs=["reads.lca.perc.phylum.tsv", "reads.lca.phylum.tsv"]) << """
Rscript summaryTable.lca.R
"""
