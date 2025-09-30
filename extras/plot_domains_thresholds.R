# Libraries
library(tidyverse)

# Path
setwd("~/Documents/Computing/eDNA/New_Analysis/")

# List folders
files <- list.files("domains", recursive = T, full.names = T, pattern = "domains.tsv")

# Read and add metadata
dfs <- lapply(files, function(x) read.table(x, header = T) %>% mutate(similarity = basename(dirname(x)), sample = str_replace(basename(x), "\\.domains\\.tsv", "")))

# Reduce
df <- Reduce(function(...) rbind(...), dfs)

# Rename domains
df <- df %>%
  mutate(label = recode(domain,
  "A" = "Archaea",
  "B" = "Bacteria",
  "E" = "Eukaryota",
  "V" = "Viruses",
  "U" = "Unclassified",
  "O" = "Other",
  "M" = "Missing")) %>%
  mutate(sample = recode(sample, 
  "GP_1" = "EN0.2A", "GP_2" = "EN0.2B", "GP_3" = "EN0.2C",
  "A0_2" = "OP0.2A", "B0_2" = "OP0.2B", "C0_2" = "OP0.2C",
  "A1_2" = "OP1.2A", "B1_2" = "OP1.2B", "C1_2" = "OP1.2C",
  "A5_0" = "OP5.0A", "B5_0" = "OP5.0B", "C5_0" = "OP5.0C",
  "A8_0" = "OP8.0A", "B8_0" = "OP8.0B", "C8_0" = "OP8.0C"))

# Add chosen threshold
df <- df %>% mutate(chosen = ifelse(similarity == 90, "Y", "N"))

# Plot
df$similarity <- factor(df$similarity, levels = c("70","75","80","85","90","95","99","100"))
df$label <- factor(df$label, levels = c("Archaea", "Viruses", "Bacteria", "Eukaryota", "Other", "Unclassified", "Missing"))
p <- df %>% filter(domain %in% c("A","V","B","E")) %>% 
  ggplot() + geom_col(aes(x=similarity, y=Count, fill=label, color=chosen)) +
  facet_wrap(sample~., nrow = 5, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c("coral3","darkolivegreen3","darkcyan","orange3")) +
  scale_color_manual(values = c("gray35","black")) +
  labs(fill = "Superkingdom", x = "Similarity threshold", y = "No. of reads") +
  guides(color = F) +
  theme_minimal(base_size = 14)

pdf("similarity_thresholds.pdf", width = 10, height = 6)
p
dev.off()
  
