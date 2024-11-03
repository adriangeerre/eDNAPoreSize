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
    outputs = [f'{out_route}/{prefix}_{n:04}.fasta' for n in range(num)]
    options = {'cores': 1,'memory': '4g', 'queue': 'short', 'walltime': '00:30:00'}

    spec = '''
    mkdir -p {out_route}
    split -a 4 -d -l 30000 01-Merge/{prefix}/{prefix}.fasta '{prefix}_' --additional-suffix=.fasta
    mv {prefix}_* {out_route}
    '''.format(prefix = prefix, out_route=out_route, num = num)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Blast
def blast_file(fainput, in_route, out_route, dbs, blastdb, folderdb, evalue, target_seqs, threads, memory):
    inputs = [f'{in_route}/{fainput}.fasta']
    outputs = [f'{out_route}/{fainput}.tsv']
    options = {'cores': 4,'memory': f'{memory}g', 'queue': 'normal', 'walltime': '24:00:00'}

    spec = '''
    mkdir -p {out_route}
    echo "Copying database $(date)"
    rsync -avr --exclude '*.tar.gz' --exclude '*.md5' --exclude 'seqs' {blastdb} /scratch/$SLURM_JOBID/
    echo "Database copied $(date)"
    blastn -query {in_route}/{fainput}.fasta -db /scratch/$SLURM_JOBID/{folderdb}/{dbs} -outfmt "6 std slen qlen sseq qseq qcovs qcovhsp qcovus staxid" -num_threads {threads} -evalue {evalue} -max_target_seqs {target_seqs} -out {out_route}/{fainput}.tmp.tsv
    # Avoid gwf finish if time limit happens
    mv {out_route}/{fainput}.tmp.tsv {out_route}/{fainput}.tsv
    '''.format(out_route=out_route, in_route=in_route, fainput=fainput, dbs=dbs, blastdb=blastdb, folderdb=folderdb, evalue=evalue, threads=threads, target_seqs=target_seqs)

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
    # ERROR: CREATE .number_of_lines.tsv (Sample\tNumReads\n)
    lines_file = ".number_of_lines.tsv"
    d = {} # Dictionary number of lines
    if os.path.exists(lines_file):
        f = open(lines_file, "r")
        for row in f:
            (key,val) = row.split()
            d[key] = int(val)
        f.close()

	# Number of files
    div_lines = 30000
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
            split_to_blast = f'{name}_{j:04}'
            gwf.target_from_template(f'blast_{split_to_blast}', blast_file(fainput = split_to_blast, in_route = f'02-Split/{name}', out_route = f'03-Blast/{name}', dbs = 'nt', blastdb=blastdb, folderdb=os.path.basename(blastdb), evalue=0.01, target_seqs=500, threads=8, memory = 4))

    # Summarise Blast
    gwf.target(f'Summarise_Blast_{name}',
    inputs=[f"03-Blast/{name}/{name}_{num:04}.tsv" for num in range(num_files)], outputs=[f"07-Other/{name}.mean_pident_length_qcovs.per_read.tsv"], cores=1, memory='2g', walltime='2:00:00') << """
    mkdir -p 07-Other
    # Pident vs Length vs Qcovs (Per Read)
    Rscript scripts/summary.blast.R -s {name}
    """.format(name=name)

    # Create Folder
    if not os.path.isdir(f"04-FilterBlast/{name}") and os.path.isdir(f"03-Blast/{name}"): os.makedirs(f"04-FilterBlast/{name}")

    # Filter and Merge Blast results
    gwf.target(f'Filter_Blast_{name}',
     inputs=[f"03-Blast/{name}/{name}_{num:04}.tsv" for num in range(num_files)],
     outputs=[f"04-FilterBlast/{name}/{name}_{num:04}_filter.tsv" for num in range(num_files)],
     cores=1, 
     memory='4g', 
     walltime='8:00:00') << """
        # Filter Blast
        for file in $(ls 03-Blast/{name}/); do cut -f 1-14,17-20 03-Blast/{name}/$file | awk '($4 >= 100)' | awk '($3 >= 90)' | awk '($15 >= 90)' | sort -k1,1 -k12,12nr -k11,11n > 04-FilterBlast/{name}/$(echo $file | cut -d "." -f 1)_filter.tsv; done

        # Merge Blast outputs
        # cat 04-FilterBlast/{name}/*_filter.tsv > 05-MergeBlast/{name}.tsv
    """.format(name=name)

f.close()

# Obtain taxids
if os.path.isdir("04-FilterBlast/"):
    dirs = os.listdir("04-FilterBlast/")
    lens = {i:len(os.listdir(f"02-Split/{i}")) for i in dirs}
    inputs = [f"04-FilterBlast/{i}/{i}_{j:04}_filter.tsv" for i in dirs for j in range(lens[i])]
    gwf.target('Obtain_taxids', 
    inputs=inputs,
    outputs=["taxids.tsv"],
    cores=1, 
    memory='2g', 
    walltime='6:00:00') << """
    # Obtain taxids
    for file in $(ls -d 04-FilterBlast/*/*)
    do
        cut -f 18 $file | sort | uniq >> taxids.tmp.tsv
    done
    cat taxids.tmp.tsv | sort | uniq > taxids.tsv
    rm taxids.tmp.tsv
    """

    # Get taxonomy (TaxizeDB)
    gwf.target('Get_taxonomy',
    inputs=["taxids.tsv"],
    outputs=["taxid2taxonomy.tsv"],
    cores=1, 
    memory='2g', 
    walltime='6:00:00') << """
    Rscript scripts/taxid2taxa.lca.R
    """