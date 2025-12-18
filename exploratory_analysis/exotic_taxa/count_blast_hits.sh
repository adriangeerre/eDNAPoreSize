#!/bin/bash
#SBATCH --account eDNA
#SBATCH --partition normal
#SBATCH --mem 4G
#SBATCH -c 1
#SBATCH -t 2:00:00

CLASS=${1}
GENUS=${2}
TYPE=${3}

python scripts/ExploreExoticTaxa.py -i 04-FilterBlast -t taxid2taxonomy.tsv --superkingdom Eukaryota --class ${CLASS} --genus ${GENUS} -o 12-ExtractTaxaReads/${TYPE}