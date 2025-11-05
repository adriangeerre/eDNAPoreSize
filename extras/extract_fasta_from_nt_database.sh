#!/bin/bash
#SBATCH --account eDNA
#SBATCH --partition normal
#SBATCH --mem 12G
#SBATCH -c 1
#SBATCH -t 24:00:00

# Variables
accefolder=${1}
dbfolder=${2}
outdir=${3}

# Source conda to work with environments
source ~/programas/minconda3.9/etc/profile.d/conda.sh
conda activate eDNA

# trap ctrl-c and call ctrl_c()
trap ctrl_c INT
function ctrl_c() {
	echo "Execution halted by user."
	exit
}

# Create folders
mkdir -p ${outdir}

# Loop taxid files
for accefile in $(ls ${accefolder}); do
	# Define taxid
	taxid=$(echo ${accefile} | sed 's/.txt//g')
	echo ${taxid}
	# Loop accessions
	for acce in $(cat ${accefolder}/${accefile} | cut -f 2); do
		# Create or replace fasta file
		blastdbcmd -db ${dbfolder}/nt -entry ${acce} -outfmt "%f" >> ${outdir}/${taxid}.fna
	done
done