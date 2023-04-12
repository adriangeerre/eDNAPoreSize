# Libraries
suppressMessages(library(taxizedb))
suppressMessages(library(tidyverse))

# Data
d <- list.dirs("05-MergeBlast")
f <- lapply(d, function(x) list.files(x, pattern = "*.tsv", full.names = T))
f <- unlist(f)
l <- lapply(f, function(x) str_replace_all(str_split(x, "/")[[1]][[2]], ".tsv", "")) %>% unlist()
names(f) <- l

# Taxonomy
taxonomy <- read.delim("taxid2taxonomy.tsv") %>% unique()

taxid2taxa <- function(df, taxonomy) {
  df <- merge(df, taxonomy, by.x="V18", by.y="taxid")
  df <- df[,c(2:18,1,20:27,19)]
  return(df)
}

# MRCA
lca <- function(df) {
  # Get reads with only a unique hit
  df <- df %>% group_by(V1) %>% mutate(count = n()) %>% ungroup()
  onehit <- df[df$count==1,] %>% select(c(V1,superkingdom:species, count))
  multhits <- df[df$count>1,] 

  # Get MRCA per read
  lca <- multhits %>% select(c(V1,superkingdom:species, count))   # Select read name and classification columns
  lca <- split(lca, f=lca$V1)  # Split dataframe into a list with entry per read
  
  # Filter Superkingdom over 90% presence of Domain (Will help to clean reads with lots of hits)
  lca <- lapply(lca, function(x) x %>% group_by(superkingdom) %>% mutate(perc_spkm = n()/nrow(x)*100) %>% filter(perc_spkm >= 90) %>% select(c(superkingdom:species, count)))
  taxa <- lapply(lca, function(x) x %>% select(-count) %>% rownames_to_column %>% gather(taxonomy, value, -rowname) %>% select(-c(rowname)) %>% group_by(value) %>% unique() %>% ungroup() %>% group_by(taxonomy) %>% count() %>% filter(n == 1) %>% select(taxonomy) %>% ungroup())
  
  # Remove unclear branches 
  for (num in 1:length(lca)) {
    lca[[num]] <- lca[[num]] %>% select(c(taxa[[num]]$taxonomy, count)) %>% unique() %>% mutate(V1 = names(lca)[[num]])
  }

  # Merge the list into a dataframe
  lca <- Reduce(function(...) merge(..., all=T), lca)
  
  # Merge Read with 1 hit and reads with several hits
  lca <- bind_rows(onehit,lca)
  
  return(lca)
}


# Wrapper
execute <- function(files, taxonomy) {
  count = 1
  pb <- txtProgressBar(min = 1, max = length(files), style = 3)
  for (file in files) {
    sample <- str_replace_all(str_split(file, "/")[[1]][[2]], ".tsv", "")
    message(sample)
    setTxtProgressBar(pb, count)
    df <- read.table(file, header=F, sep="\t") %>% mutate(file = sample)
    df <- df %>% filter(V11 < 0.0001, V4 >= 100, V3 >= 90, V15 >= 90)
    df <- taxid2taxa(df, taxonomy)
    df <- suppressMessages(lca(df))
    count += 1
  }
  close(pb)
  return(df)
}

# Execute in loop
for (smpl in names(f)) {
  z <- execute(f[[smpl]], taxonomy)
  write.table(z,file=paste("06-LCA/",smpl,".lca.tsv", sep=""),quote=F,sep="\t",row.names=F)
  rm(z)
  message("")
}


