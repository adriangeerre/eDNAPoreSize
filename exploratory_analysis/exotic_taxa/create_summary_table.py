# Pre-run on Linux: ls 12-ExtractTaxaReads/*/*/*.tsv | xargs -I
# {} sh -c 'printf "{}\t$(cut -f 1 {} | grep -v qseqid | sort | uniq |
# wc -l)\n"' > 12-ExtractTaxaReads/read_counts.tsv &

# Imports
import pandas as pd

# Read count per file
infile = "12-ExtractTaxaReads/read_counts.tsv"

# Load file
df = pd.read_table(infile, header=None)

# Split path and adjust columns
df = pd.concat([df[0].str.split("/", expand=True), df[1]], axis=1)
df.columns = ['folder','type','filter','name','value']

# Pivot
final_df = df.pivot(columns=["type","name"], index="filter", values="value").T
final_df.to_csv("12-ExtractTaxaReads/final.read_counts.tsv", sep="\t", index=True)