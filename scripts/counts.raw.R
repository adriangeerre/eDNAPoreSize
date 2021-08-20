# Libraries
library(tidyverse)

# Data
f <- list.files("blast_1st_merged/")
dfs <- lapply(f, function(x) read.csv2(paste("blast_1st_merged/", x, sep=""), header=F, sep="\t") %>% mutate(file = str_replace(x,".tsv","")))
taxonomy <- read.delim("taxid2taxonomy.tsv")

# Add taxonomy
dfs <- lapply(dfs, function(x) merge(x, taxonomy, by.x="V16", by.y="taxid"))

# Transform data
dfs <- lapply(dfs, function(x) x[,c(18:25,17)])
dfs <- lapply(dfs, function(x) x %>% group_by(superkingdom,kingdom,phylum,order,class,family,genus,species,file) %>% count() %>% spread(file, n))

# Group data
df <- Reduce(function(...) merge(..., all=T), dfs)
df[9:23][is.na(df[9:23])] <- 0

# Summarize counts
df <- df %>% group_by(superkingdom,kingdom,phylum,class,order,family,genus,species) %>% summarise_if(is.numeric, sum)

# Arrange
df <- df[c(1:8,9,13,17,10,14,18,11,15,19,12,16,20,21,22,23)]

# Write table
write.table(df,file="raw.taxa.count.tsv",quote=F,sep="\t",row.names=F)


