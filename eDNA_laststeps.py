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

# Execution
#----------

# Collect sample names
samples = []

# Open file again 
f = open(info_file,'r')

for sample in f:
    
    if sample[0] == '#':
        continue

    # Prepare the data to be given to the functions
    sample = sample.strip().split("\t")
    name = sample[2]	# name

    # Append sample names
    samples.append(name)

    # Lens
    dirs = os.listdir("04-FilterBlast/")
    lens = {i:len(os.listdir(f"02-Split/{i}")) for i in dirs}

    # Number of files
    num_files = lens[name]

    # Execution
    #----------

    # Get LCA from original Blast (TaxizeDB) -> Filtered for e-value < 0.0001, algn len >= 100, pident >= 90 and qcovs >= 90
    if not os.path.isdir(f"05-LCA/{name}") and os.path.isdir(f"04-FilterBlast/{name}"): os.makedirs(f"05-LCA/{name}")

    gwf.target(f'Blast_LCA_{name}',
    inputs= [f"04-FilterBlast/{name}/{name}_{num:04}_filter.tsv" for num in range(num_files)] + ["taxid2taxonomy.tsv"],
    outputs=[f"05-LCA/{name}/{name}_{num:04}_filter.lca.tsv" for num in range(num_files)],
    cores=4, 
    memory='16g', 
    walltime='48:00:00') << """
    Rscript scripts/lca.R -s {name} -t 4
    """.format(name=name)

    # Merge
    gwf.target(f'Merge_LCA_{name}',
    inputs=[f"05-LCA/{name}/{name}_{num:04}_filter.lca.tsv" for num in range(num_files)],
    outputs=[f"06-mergeLCA/{name}.lca.tsv"],
    cores=1, 
    memory='4g', 
    walltime='48:00:00') << """
    mkdir -p 06-mergeLCA/
    head -n 1 05-LCA/{name}/{name}_0000_filter.lca.tsv > 06-mergeLCA/{name}.lca.tsv
    cat 05-LCA/{name}/{name}_*_filter.lca.tsv | grep -wv "superkingdom" >> 06-mergeLCA/{name}.lca.tsv
    """.format(name=name)

f.close()

# Generate species-count matrix (TaxizeDB)
gwf.target('Count_matrix',
inputs=[f"06-mergeLCA/{name}.lca.tsv" for name in samples],
outputs=["lca.taxa.count.tsv"],
cores=1,
memory='4g',
walltime='2:00:00') << """
Rscript scripts/counts.lca.R
"""

# Rarefy (ROBITools)
gwf.target('Rarefy',
inputs=["lca.taxa.count.tsv"],
outputs=["Pore_size_table_for_rarefy.tsv","counts.lca.rarefy.tsv"],
cores=1,
memory='4g',
walltime='2:00:00') << """
Rscript scripts/rarefyROBI.R
"""

# Obtain summary table and heatmap (Bioconductor)
gwf.target('Summary_table_and_heatmap',
inputs=["counts.lca.rarefy.tsv"],
outputs=["reads.lca.perc.phylum.tsv", "reads.lca.phylum.tsv"],
cores=1,
memory='4g',
walltime='2:00:00') << """
Rscript scripts/summaryTable.lca.R
"""