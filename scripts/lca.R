# Libraries
library(taxizedb)
library(tidyverse)

# Data
d <- list.dirs("05-MergeBlast/")
f <- lapply(d, function(x) list.files(x, pattern = "*.tsv", full.names = T))
l <- lapply(d, function(x) str_sub(str_split(x, "/")[[1]][3], 1, 4)) %>% unlist()
names(f) <- l

# Taxonomy
taxonomy <- read.delim("taxid2taxonomy.tsv") %>% unique()

taxid2taxa <- function(df, taxonomy) {
  df <- merge(df, taxonomy, by.x="X20", by.y="taxid")
  df <- df[,c(2:16,1,18:25,17)]
  return(df)
}

# MRCA
lca <- function(df) {
  # Get reads with only a unique hit
  df <- df %>% group_by(X1) %>% mutate(count = n()) %>% ungroup()
  onehit <- df[df$count==1,] %>% select(c(X1,superkingdom:species))
  multhits <- df[df$count>1,]

  # Get MRCA per read
  lca <- multhits %>% select(c(X1,superkingdom:species))   # Select read name and classification columns
  lca <- split(lca, f=lca$X1)  # Split dataframe into a list with entry per read
  
  # Filter Superkingdom over 90% presence of Domain (Will help to clean reads with lots of hits)
  lca <- lapply(lca, function(x) x %>% group_by(superkingdom) %>% mutate(perc_spkm = n()/nrow(x)*100) %>% filter(perc_spkm >= 90) %>% select(c(superkingdom:species)))
  taxa <- lapply(lca, function(x) x %>% rownames_to_column %>% gather(taxonomy, value, -rowname) %>% select(-c(rowname)) %>% group_by(value) %>% unique() %>% ungroup() %>% group_by(taxonomy) %>% count() %>% filter(n == 1) %>% select(taxonomy) %>% ungroup())
  
  # Remove unclear branches 
  for (num in 1:length(lca)) {
    lca[[num]] <- lca[[num]] %>% select(c(taxa[[num]]$taxonomy)) %>% unique() %>% mutate(X1 = names(lca)[[num]])
  }

  # Merge the list into a dataframe
  lca <- Reduce(function(...) merge(..., all=T), lca)
  
  # Merge Read with 1 hit and reads with several hits
  lca <- bind_rows(onehit,lca)
  
  return(lca)
}


# Wrapper
execute <- function(files) {
  message(str_sub(str_split(files[1], "/")[[1]][3], 1, 4))
  count = 1
  pb <- txtProgressBar(min = 1, max = length(files), style = 3)
  dtax <- data.frame()
  for (file in files) {
    setTxtProgressBar(pb, count)
    df <- read_csv(file, col_names=F, show_col_types=F, progress=F) %>% mutate(file = str_sub(str_split(file, "/")[[1]][4], 1, 4)) %>% select(-c(X15,X16))
    df <- df %>% filter(X11 < 0.0001, X14 >= 100, X3 >= 90)
    df <- taxid2taxa(df, taxonomy)
    df <- suppressMessages(lca(df))
    count = count + 1
    dtax <- bind_rows(dtax, df)
  }
  close(pb)
  return(dtax)
}

# Execute in loop
for (smpl in names(f)) {
  z <- execute(f[[smpl]])
  write.table(z,file=paste("06-LCA/",smpl,".lca.tsv", sep=""),quote=F,sep="\t",row.names=F)
  rm(z)
  message("")
}


