#!/bin/bash
#SBATCH --account eDNA
#SBATCH --partition normal
#SBATCH --mem 8G
#SBATCH -c 8
#SBATCH -t 4:00:00

# Variables
infolder=${1}
primers=${2}
outdir=${3}
threads=8

# Source conda to work with environments
source ~/programas/minconda3.9/etc/profile.d/conda.sh
conda activate Mothur

# trap ctrl-c and call ctrl_c()
trap ctrl_c INT
function ctrl_c() {
	echo "Execution halted by user."
	exit
}

# Create folders
mkdir -p ${outdir}

# Loop fastas
for fna in $(ls -d ${infolder}/*.fna); do
	# Define name
	name=$(basename ${fna} | sed 's/.fna//g')
	# Avoid recomputation
	if [ ! -f ${infolder}/${name}.pcr.fna ]; then
		# Run Mothur
		mothur "#pcr.seqs(fasta=${fna}, oligos=${primers}, pdiffs=3, rdiffs=3, keepdots=False, processors=${threads})"
		# Move output files
		mv ${infolder}/${name}.pcr.fna ${outdir}
		mv ${infolder}/${name}.bad.accnos ${outdir}
		mv ${infolder}/${name}.scrap.pcr.fna ${outdir}
	fi
done
