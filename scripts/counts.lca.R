# Libraries
library(tidyverse)

# Data
f <- list.files(path = "06-LCA/", pattern = "*.lca.tsv", full.names = TRUE)
dfs <- lapply(f, function(x) read.csv2(x, header=T, sep="\t") %>% mutate(file = str_replace(str_replace(x,".lca.tsv",""),"lca//","")))

# Transform data
dfs <- lapply(dfs, function(x) x %>% select(-X1) %>%  group_by(superkingdom,kingdom,phylum,order,class,family,genus,species,file) %>% count() %>% spread(file, n))

# Group data
df <- Reduce(function(...) merge(..., all=T), dfs)
df[9:23][is.na(df[9:23])] <- 0

# Arrange
df <- df[c(1:8,9,13,17,10,14,18,11,15,19,12,16,20,21,22,23)]

# Write table
write.table(df,file="lca.taxa.count.tsv",quote=F,sep="\t",row.names=F)
