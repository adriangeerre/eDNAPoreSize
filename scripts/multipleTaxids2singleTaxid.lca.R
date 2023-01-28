# Libraries
library(taxizedb)
library(tidyverse)

# Data
df <- read.delim("multiple.taxids.tsv",header=F)

# Function
mrcaTaxids <- function(taxids) {
  taxs <- strsplit(taxids, split=";", fixed=T)[[1]]
  d <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(d) <- c("name", "rank", "id")
  for (t in taxs) {
    if (suppressMessages(is.na(classification(t)[[1]])) == TRUE) {
      next
    } else {
      d <- rbind(d, classification(t)[[1]])
    }
    colnames(d) <- c("name", "rank", "id")
  }
  d <- d %>% group_by(id) %>% mutate(count = n()) %>% ungroup() %>% filter(count >= max(count)) %>% unique() %>% slice(nrow(.)) %>% select(id)
  if (nrow(d) == 0) {
    return(NA)
  } else {
    return(d$id) 
  }
}

# Result
r <- c()
pb <- txtProgressBar(min = 1, max = length(df$V1), style = 3)
count = 1
for (taxids in df$V1) {
  setTxtProgressBar(pb, count)
  r <- append(r, mrcaTaxids(taxids))
  count = count + 1
}
close(pb)

# Save results
rdf <- df %>% mutate(mrca = r) %>% rename(taxids = V1)
write.table(rdf,file="m2s.taxids.tsv",quote=F,sep="\t",row.names=F)
