# Libraries
library(tidyverse)
library(optparse)

# Parameters
# ----------
option_list = list(
    make_option(c("-s", "--sample"), action="store", default=NA, type='character', help="Sample name"))
opt = parse_args(OptionParser(option_list=option_list))

# Output
out_file <- paste("07-Other/", opt$sample, ".mean_pident_length_qcovs.per_read.tsv", sep="")

# Files
files <- list.files(paste("03-Blast/", opt$sample, sep=""), full.names=T)

# Loop
for (f in files) {
    df <- read.table(f, header=F, sep="\t")

    # Summarize
    df <- df %>% group_by(V1) %>% summarise(V3 = mean(V3), V4 = mean(V4), V17 = mean(V17)) %>% ungroup() %>% as.data.frame()

    # Save output
    if (paste(opt$sample, ".mean_pident_length_qcovs.per_read.tsv", sep="") %in% list.files("07-Other/")) {
        write.table(df, out_file, sep="\t", col.names=F, row.names=F, append=T, quote=F)
    } else {
        write.table(df, out_file, sep="\t", col.names=F, row.names=F,  append=F, quote=F)
    }
}
