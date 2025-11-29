#!/bin/bash
#SBATCH --account=eDNA
#SBATCH --time=12:00:00
#SBATCH --partition=normal
#SBATCH --mem=8G
#SBATCH -c 1

# Variables
infolder=${1}
r18Sfolder=${2}
outfolder=${3}

# Exec
python scripts/FilterEukaryotaReadsThreshold.py -i ${infolder} -r ${r18Sfolder} -o ${outfolder}