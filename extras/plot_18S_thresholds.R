# Libraries
library(tidyverse)
library(reshape2)

# Path
setwd("~/Documents/Computing/eDNA/New_Analysis/")

# Define total number of eukaryota reads
define_euka_reads_per_sample_and_threshold <- function() {
  files <- list.files("09-taxid_thresholds/read_counts", recursive = T, full.names = T)
  dfs <- lapply(files, function(x) read.table(x, header = T))
  df <- Reduce(function(...) rbind(...), dfs)
  total_euka_reads <- df %>% group_by(superkingdom,sample,threshold) %>% summarise(read_counts = sum(read_counts)) %>% filter(superkingdom == "E") %>% mutate(index = paste(sample,threshold, sep="_")) %>% ungroup() %>% select(index, read_counts)
  return(total_euka_reads)
}

total_euka_reads <- define_euka_reads_per_sample_and_threshold()

# List folders
files <- list.files("11-Blast_r18S/r18S_counts/", recursive = T, full.names = T)

# Read and add metadata
dfs <- lapply(files, function(x) read.table(x, header = T))

# Reduce
df <- Reduce(function(...) rbind(...), dfs)

# Add total eukaryota reads
df <- df %>% pivot_wider(id_cols = c("sample","threshold"), values_from = count, names_from = type) %>% mutate(index = paste(sample,threshold, sep="_"))
df <- merge(df, total_euka_reads, by = "index") %>% mutate("Non-18S" = read_counts - Accession - Blast) %>% select(-index, -read_counts) %>% melt(id.vars = c("sample", "threshold"))

# Rename domains
df <- df %>%
  mutate(sample = recode(sample, 
                         "GP_1" = "EN0.2A", "GP_2" = "EN0.2B", "GP_3" = "EN0.2C",
                         "A0_2" = "OP0.2A", "B0_2" = "OP0.2B", "C0_2" = "OP0.2C",
                         "A1_2" = "OP1.2A", "B1_2" = "OP1.2B", "C1_2" = "OP1.2C",
                         "A5_0" = "OP5.0A", "B5_0" = "OP5.0B", "C5_0" = "OP5.0C",
                         "A8_0" = "OP8.0A", "B8_0" = "OP8.0B", "C8_0" = "OP8.0C"))

# Add chosen threshold
df <- df %>% mutate(chosen = ifelse(threshold == 90, "Y", "N"))

# Count total 18S reads
eu18S <- df %>% filter(variable != "Non-18S") %>% group_by(sample, threshold) %>% summarise(r18S = sum(value))
pos18S <- df %>% group_by(sample, threshold) %>% summarise(total = sum(value)) %>% ungroup() %>% group_by(sample, threshold) %>% summarise(position = max(total)*0.35)
eu18S <- merge(eu18S, pos18S, by="sample")

# Plot
df$threshold <- factor(df$threshold, levels = c("70","75","80","85","90","95","99","100"))
df$variable <- factor(df$variable, levels = c("Non-18S", "Blast", "Accession"))

p <- df %>% 
  ggplot() + geom_col(aes(x=threshold, y=value, fill=variable, color=chosen, group=threshold)) +
  geom_text(data = eu18S, aes(x=threshold, y=r18S+10000, label = paste("n=",r18S,sep="")), color="darkred", size=2, angle=90, hjust = 0) +
  facet_wrap(sample~., nrow = 5, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c("gray90", "darkgoldenrod2","indianred2")) +
  scale_color_manual(values = c("gray35","black")) +
  scale_y_continuous(labels = scales::label_scientific()) +
  labs(fill = "Identification", x = "Similarity threshold", y = "No. of Eukaryota reads") +
  guides(color = F) +
  theme_minimal(base_size = 14)

pdf("18S_reads.pdf", width = 10, height = 6)
p
dev.off()
