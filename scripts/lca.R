# Libraries
suppressMessages(library(taxizedb))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(foreach))

# Parameters
# ----------
option_list = list(
    make_option(c("-s", "--sample"), action="store", default=NA, type='character', help="Sample name"),
    make_option(c("-t", "--threads"), action="store", default=NA, type='numeric', help="Sample name"))
opt = parse_args(OptionParser(option_list=option_list))

# Data
f <- list.files(paste("04-FilterBlast", opt$sample, sep="/"), pattern=".tsv", full.names=T)

# Taxonomy
taxonomy <- read.delim("taxid2taxonomy.tsv") %>% unique()

taxid2taxa <- function(df, taxonomy) {
  df <- merge(df, taxonomy, by.x="V18", by.y="taxid")
  df <- df[,c(2:18,1,20:27,19)]
  return(df)
}

# MRCA
lca <- function(df, threads) {
  # Get reads with only a unique hit
  message("\tCount")
  df <- df %>% group_by(V1) %>% mutate(count = n()) %>% ungroup()
  onehit <- df[df$count==1,] %>% select(c(V1,superkingdom:species, count))
  multhits <- df[df$count>1,] 

  # Get MRCA per read
  lca <- multhits %>% select(c(V1,superkingdom:species, count)) # Select read name and classification columns
  message("\tSplit")
  lca <- split(lca, f=lca$V1) # Split dataframe into a list with entry per read
  
  # Filter Superkingdom over 90% presence of Domain (Will help to clean reads with lots of hits)
  message("\tFilter")
  lca <- lapply(lca, function(x) x %>% group_by(superkingdom) %>% mutate(perc_spkm = n()/nrow(x)*100) %>% filter(perc_spkm >= 90) %>% select(c(superkingdom:species, count)))

  # Taxonomy
  taxa <- lapply(lca, function(x) x %>% select(-count) %>% rownames_to_column %>% gather(taxonomy, value, -rowname) %>% select(-c(rowname)) %>% group_by(value) %>% unique() %>% ungroup() %>% group_by(taxonomy) %>% count() %>% filter(n == 1) %>% select(taxonomy) %>% ungroup())
  
  message("\tParrallel")
  # Remove unclear branches (Parallel)
  cl <- parallel::makeForkCluster(threads)
  doParallel::registerDoParallel(cl)
  lca <- foreach(num = seq_along(lca)) %dopar% {
    lca[[num]] %>% select(c(taxa[[num]]$taxonomy, count)) %>% unique() %>% mutate(V1 = names(lca)[[num]])
  }
  parallel::stopCluster(cl)

  # Merge the list into a dataframe
  lca <- Reduce(function(...) merge(..., all=T), lca)
  
  # Merge Read with 1 hit and reads with several hits
  lca <- bind_rows(onehit,lca)
  
  return(lca)
}

# Wrapper
execute <- function(file, taxonomy, sample, threads) {
  message("Read")
  df <- read.table(file, header=F, sep="\t") %>% mutate(file = sample)
  message("Filter")
  df <- df %>% filter(V11 < 0.0001, V4 >= 100, V3 >= 90, V15 >= 90)
  message("Taxid2taxa")
  df <- taxid2taxa(df, taxonomy)
  message("LCA")
  df <- suppressMessages(lca(df, threads))
  return(df)
}

# Execute in loop
for (tsv in f) {
  message(paste("Running file: ", tsv))
  name=str_remove(unlist(str_split(tsv, "/"))[3],".tsv")
  z <- execute(tsv, taxonomy, opt$sample, as.numeric(opt$threads))
  write.table(z, file=paste("05-LCA/", opt$sample, "/", name, ".lca.tsv", sep=""), quote=F, sep="\t", row.names=F)
}

