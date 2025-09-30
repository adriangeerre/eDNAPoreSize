#!/bin/bash
#SBATCH --account=eDNA
#SBATCH --time=24:00:00
#SBATCH --partition=normal
#SBATCH --mem=8G
#SBATCH -c 1

# Variables
infolder=${1}
outfolder=${2}

# Exec
python scripts/CountUniqueTaxids.py -i ${infolder} -o ${outfolder}