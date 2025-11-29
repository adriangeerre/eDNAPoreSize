#!/bin/bash
#SBATCH --account=eDNA
#SBATCH --time=12:00:00
#SBATCH --partition=normal
#SBATCH --mem=8G
#SBATCH -c 1

# Variables
infolder=${1}
categories=${2}
outfolder=${3}

# Exec
python scripts/MatchReadsSuperkingdom.py -i ${infolder} -c ${categories} -o ${outfolder}