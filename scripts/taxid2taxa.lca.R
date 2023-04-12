# Libraries
suppressMessages(library(taxizedb))
suppressMessages(library(tidyverse))
suppressMessages(library(plyr))

# Data
taxa <- read.delim("taxids.tsv",header=F)

# Group
ids <- taxa$V1

# Colnames
tc <- c("taxid","superkingdom","kingdom","phylum","class","order","family","genus","species")

# Taxonomy identification
t <- classification(ids)
t <- t[!is.na(t)]
t <- lapply(t, function(x) x %>% mutate(taxid = x$id[nrow(x)]) %>% filter(rank %in% tc) %>% group_by(rank) %>% filter(id == min(as.numeric(id))) %>% ungroup() %>% select(-id) %>% unique() %>% spread(rank,name))

# Merge classifications
vals <- data.frame(taxid=NA, class=NA, family=NA, genus=NA, kingdom=NA, order=NA, phylum=NA, species=NA, superkingdom=NA)
tf <- lapply(t, function(x) rbind.fill(x, vals)[1,]) # Faster than Reduce!
tf <- lapply(tf, function(x) x %>% select(all_of(tc))) # Sort to mix
df <- data.frame(matrix(unlist(tf), nrow=length(tf), byrow=TRUE),stringsAsFactors=FALSE) # Requires same number of columns in each data frame
colnames(df) <- tc
df <- df[!is.na(df$taxid),]
x <- df %>% arrange(superkingdom,kingdom,phylum,class,order,family,genus,species)

# Write table
write.table(x,file="taxid2taxonomy.tsv",quote=F,sep="\t",row.names=F)
