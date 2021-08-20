# Summary BLAST into species count (.tsv)

# Basic path
cd /datos/GenomeDK

# Obtain data
mkdir /datos/GenomeDK/blast/
cd blast/
rsync -avr agomez@login.genome.au.dk:/home/agomez/eDNA/faststorage/Adrian/Tesis/results/blast/* .

# Obtain BLAST output. Remove sseq, qseq y ssciname. Transform to tab-delimited file. Total: 16 columns.
for folder in $(ls)
do
	mkdir ../blast_1st/${folder}
 	cd $folder
 	for file in $(ls | cut -d '.' -f 1,2)
 	do
 		cat $file.csv | cut -d ',' -f 1-14,17-18 | sed 's/,/\t/g' | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > ../../blast_1st/${folder}/${file}_1st.tsv
	done
  cd ../
done
cd ../

# Merge BLAST per sample
mkdir blast_1st_merged
cd /datos/GenomeDK/blast_1st/
for folder in $(ls)
do
	cd ${folder}; cat * > $(echo ${folder} | cut -d '_' -f 1,2).tsv; mv $(echo ${folder} | cut -d '_' -f 1,2).tsv ../../blast_1st_merged/
	cd ../
done
cd ../

# Obtain taxids
cd blast_1st_merged
for file in $(ls)
do
	cut -f 16 ${file} | sort | uniq >> ../taxids.tmp.tsv
done
cd ../
cat taxids.tmp.tsv | sort | uniq > taxids.tsv
rm taxids.tmp.tsv

# Divide taxids between single and multiple (semicolons separated)
grep ';' taxids.tsv > multiple.taxids.tsv
grep -v ';' taxids.tsv > single.taxids.tsv

# MRCA for multiple taxids (Results: m2s.taxids.tsv)
/programas/miniconda3/envs/TaxizeDB/bin/Rscript mulipleTaxids2singleTaxid.R

# Get taxonomy (Results: taxid2taxonomy.tsv)
/programas/miniconda3/envs/TaxizeDB/bin/Rscript taxid2taxa.R

# Generate species-count matrix (Results: raw.taxa.count.tsv)
/programas/miniconda3/envs/TaxizeDB/bin/Rscript counts.raw.R

# Rarefy (Results: Pore_size_table_for_rarefy.tsv & counts.rarefy.tsv)
/programas/miniconda3/envs/ROBITools/bin/Rscript rarefyROBI.R

# Obtain summary table and heatmap (Results: reads.perc.phylum.tsv & reads.phylum.tsv)
/programas/miniconda3/envs/Bioconductor/bin/Rscript summaryTable.R