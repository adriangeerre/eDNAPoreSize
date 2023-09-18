# Libraries
library(tidyverse)
library(VennDiagram)

# --------------
## Metabarcoding
# --------------
metabar_taxa  <- read.delim("~/Documentos/GenomeDK_Metabar/classified.txt")
metabar_count <- read.delim("~/Documentos/GenomeDK_Metabar/DADA2_nochim.table") 

# Select controls
m_cne <- grepl("CNE", names(metabar_count))
m_cne[1] <- TRUE
m_ntc <- grepl("NTC", names(metabar_count))
m_ntc[1] <- TRUE

cne <- metabar_count[,m_cne]
ntc <- metabar_count[,m_ntc]

# Find contaminants
cne_cont <- cne %>% mutate(total = rowSums(cne[,-1])) %>% filter(total > 0)
ntc_cont <- ntc %>% mutate(total = rowSums(ntc[,-1])) %>% filter(total > 0)

# Select contaminant sequences
cont <- unique(c(cne_cont$X, ntc_cont$X))

# Remove contaminants and filter
metabarcoding <- metabar_taxa %>% filter(!qseqid %in% cont) %>% filter(pident.max.best > 90)

# Select contaminants
contaminants <- metabar_taxa %>% filter(qseqid %in% cont)

# Metabarcoding diversity
tx <- metabarcoding[,c(8:14)]
dk <- length(unique(tx$kingdom[!is.na(tx$kingdom)]))
dp <- length(unique(tx$phylum[!is.na(tx$phylum)]))
dc <- length(unique(tx$class[!is.na(tx$class)]))
do <- length(unique(tx$order[!is.na(tx$order)]))
df <- length(unique(tx$family[!is.na(tx$family)]))
dg <- length(unique(tx$genus[!is.na(tx$genus)]))
ds <- length(unique(tx$species[!is.na(tx$species)]))

message(paste("Metabarcoding diversity: kingdom: ", dk, ", phylum: ", dp, ", class: ", dc, ", order: ", do, ", family: ", df, ", genus: ", dg, " and species: ", ds, sep=""))

# --------
## Shotgun
# --------
shotgun_taxa_count <- read.delim("~/Documentos/GenomeDK_MRCA/counts.mrca.rarefy.tsv")

# Clean shotgun
shotgun <- shotgun_taxa_count %>% filter(!species %in% contaminants$species)

# Select Eukaryota
shotgun_eu <- shotgun %>% filter(superkingdom == "Eukaryota")

# -----------
## Comparison
# -----------
phyl_metabar <- unique(metabarcoding$phylum[!is.na(metabarcoding$phylum)])
phyl_shotgun <- unique(shotgun_eu$phylum[!is.na(shotgun_eu$phylum)])
taxa <- list(metabarcoding = phyl_metabar, shotgun_eu = phyl_shotgun)

# Function
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

# Venn Diagram (Saved manually)
display_venn(taxa,
             category.names = c("Shotgun", "Metabarcoding"),
             fill = c("royalblue", "gold"),
             cex = 2,
             cat.cex = 2.5, cat.fontface = "bold", 
             cat.dist = c(0.15, 0.15),
             scaled = F)

# Unique phyla
unq_metabar <- sort(taxa$metabarcoding[!taxa$metabarcoding %in% taxa$shotgun_eu])
unq_shotgun <- sort(taxa$shotgun_eu[!taxa$shotgun_eu %in% taxa$metabarcoding])

# Phylum table
df <- data.frame(phylum = unique(c(phyl_metabar, phyl_shotgun)))

df <- df %>% mutate(meta = ifelse(df$phylum %in% phyl_metabar, "metabar", NA)) %>%
  mutate(shot = ifelse(df$phylum %in% phyl_shotgun, "shotgun", NA)) %>% 
  mutate(metabarcoding = ifelse(meta == "metabar", phylum, NA)) %>% 
  mutate(shotgun = ifelse(shot == "shotgun", phylum, NA)) %>% 
  select(metabarcoding, shotgun) %>% 
  arrange(shotgun)

write.table(df, "~/Documentos/GenomeDK_Metabar/phylum_comparison_table.tsv", quote = F, row.names = F, sep="\t")
