# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By    : Adrián Gómez Repollés
# Email         : adrian.gomez@mbg.au.dk
# Created Date  : 30/09/2025
# version       : '1.0'
# ---------------------------------------------------------------------------
""" Pipeline (gwf) to perform blast merged reads against r18S sequences. """ 
# ---------------------------------------------------------------------------
# Imports 
# ---------------------------------------------------------------------------

import os
import sys
import glob
from gwf import Workflow, AnonymousTarget

# Parameters
#-----------

account = "eDNA"
infolder = "../10-Estimate18S/fasta_70"
blastdb = "blastdb"
outfolder = "blast_results"

# Workflow start
#---------------
gwf = Workflow(defaults={"account": account})

# Functions
#----------

# Blast
def blast_file(name, fasta, outdir, dbs, blastdb, folderdb, target_seqs, threads, memory):
	inputs = [fasta]
	outputs = [f'{outdir}/{name}_r18S.tsv']
	options = {'cores': 4,'memory': f'{memory}g', 'queue': 'normal', 'walltime': '24:00:00'}

	spec = '''
	echo "Copying database $(date)"
	rsync -avr {blastdb} /scratch/$SLURM_JOBID/
	echo "Database copied $(date)"
	blastn -query {fasta} -db /scratch/$SLURM_JOBID/{folderdb}/{dbs} -outfmt "6 std slen qlen sseq qseq qcovs qcovhsp qcovus" -num_threads {threads} -max_target_seqs {target_seqs} -out {outdir}/{name}_r18S.tmp.tsv
	# Avoid gwf finish if time limit happens
	mv {outdir}/{name}_r18S.tmp.tsv {outdir}/{name}_r18S.tsv
	'''.format(name=name, fasta=fasta, outdir=outdir, dbs=dbs, blastdb=blastdb, folderdb=folderdb, threads=threads, target_seqs=target_seqs)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Execution
#----------

# Define files
files = sorted(glob.glob(f"{infolder}/*.other_euka.fasta"))

# Loop files
for infile in files:
	# Define name
	name = os.path.basename(infile).replace(".other_euka.fasta","")
	# Create output folder
	if not os.path.isdir(outfolder): os.makedirs(outfolder)
	# Run target
	gwf.target_from_template(f'blast_r18S_{name}', blast_file(name=name, fasta=infile, outdir=outfolder, dbs="r18S", blastdb=blastdb, folderdb=os.path.basename(blastdb), target_seqs=50, threads=8, memory=4))