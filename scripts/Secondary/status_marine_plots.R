library(tidyverse)
library(reshape2)

# Data
status <- read.table("current_status.txt", header=F, sep="\t")
marine <- read.table("marine.tsv", header=F, sep="\t")
divers <- read.delim("arter.tsv", header=T, sep="\t")

## Status
#--------

# Transform data
status <- status %>% mutate(completed = (V3/V2) * 100)

# Plot
status %>% ggplot() + geom_col(aes(x=V1, y=completed, fill=completed)) +
  scale_fill_gradient2(low = "white", mid = "orange", high = "darkred", midpoint = 50) +
  labs(x = "samples", y = "Files Blasted") +
  geom_hline(yintercept = 100, linetype=2) +
  theme_classic()

## Database
#----------

# Full
df <- divers %>% select(Række, Klasse, Orden, Familie, Slægt, Arter, Marine)
dffull <- melt(df, id.vars = "Marine") %>% unique()
dffull <- dfm %>% group_by(variable) %>% summarise(count_full = n())
dffull$variable <- c("Kingdom", "Class", "Order", "Family", "Genus", "Species")
dffull$variable <- factor(dfmf$variable, levels = c("Kingdom", "Class", "Order", "Family", "Genus", "Species"))

# Marine
df <- divers %>% select(Række, Klasse, Orden, Familie, Slægt, Arter, Marine) %>% filter(Marine == "x")
dfmarine <- melt(df, id.vars = "Marine") %>% unique()
dfmarine <- dfmarine %>% group_by(variable) %>% summarise(count_marine = n())
dfmarine$variable <- c("Kingdom", "Class", "Order", "Family", "Genus", "Species")
dfmarine$variable <- factor(dfmf$variable, levels = c("Kingdom", "Class", "Order", "Family", "Genus", "Species"))

# Merge
dfm <- merge(dffull, dfmarine, by="variable")

# Plot
dfm %>% ggplot() + geom_col(aes(x=variable, y=count_full, fill=variable)) +
  geom_col(aes(x=variable, y=count_marine, fill=variable), color="black", ) +
  geom_text(aes(x=variable, y=count_marine, label=count_marine)) +
  labs(x="Taxonomic level", y="Count", fill="", color="Marine") +
  theme_classic()


