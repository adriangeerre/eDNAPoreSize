# Libraries
library(taxizedb)
library(tidyverse)

# Data
single <- read.delim("single.taxids.tsv",header=F)
m2s <- read.delim("m2s.taxids.tsv",header=T)

# Group
ids <- c(single$V1,m2s$mrca)

# Taxonomy
tc <- c("taxid","superkingdom","kingdom","phylum","class","order","family","genus","species")
t <- classification(ids)
table(is.na(t))
t <- t[!is.na(t)]
t <- lapply(t, function(x) x %>% mutate(taxid = x$id[nrow(x)]) %>% filter(rank %in% tc) %>% select(-id) %>% unique() %>% spread(rank,name))
df <- Reduce(function(...) merge(..., all=T), t) # This takes a while (a lot!)
x <- df %>% arrange(superkingdom,kingdom,phylum,class,order,family,genus,species)

# Write table
write.table(x,file="taxid2taxonomy.tsv",quote=F,sep="\t",row.names=F)
