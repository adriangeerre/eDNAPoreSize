# Libraries (Venn Diagram)
library(tidyverse)
library(VennDiagram)

# Libraries (NMDS)
library(ROBITools)
library(reshape2)
library(vegan)
library(ggalt)
library(ggpubr)
library(ggrepel)

# Libraries (Rarefaction)
library(vegan)

# --------------
## Metabarcoding
# --------------
# Data
metabar_taxa  <- read.delim("GenomeDK_Metabar/classified.txt")
metabar_count <- read.delim("GenomeDK_Metabar/DADA2_nochim.table")

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
cont <- unique(c(cne_cont$id, ntc_cont$id))

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

# Remove contaminant from counts
metabar_count_clean <- metabar_count %>% filter(!id %in% cont)
metabar_count_clean <- metabar_count_clean[,!grepl("NTC", names(metabar_count_clean))]
metabar_count_clean <- metabar_count_clean[,!grepl("CNE", names(metabar_count_clean))]

# Rename columns
cols <- colnames(metabar_count_clean)

clean_cols <- str_replace(str_replace(cols[-1], "sample.", ""),"u","")
clean_cols <- str_replace(clean_cols, "GP1", "EN0.2A")
clean_cols <- str_replace(clean_cols, "GP2", "EN0.2B")
clean_cols <- str_replace(clean_cols, "GP3", "EN0.2C")
clean_cols <- str_replace(clean_cols, "1_0.2", "OP0.2A")
clean_cols <- str_replace(clean_cols, "2_0.2", "OP0.2B")
clean_cols <- str_replace(clean_cols, "3_0.2", "OP0.2C")
clean_cols <- str_replace(clean_cols, "1_1.2", "OP1.2A")
clean_cols <- str_replace(clean_cols, "2_1.2", "OP1.2B")
clean_cols <- str_replace(clean_cols, "3_1.2", "OP1.2C")
clean_cols <- str_replace(clean_cols, "1_5.0", "OP5.0A")
clean_cols <- str_replace(clean_cols, "2_5.0", "OP5.0B")
clean_cols <- str_replace(clean_cols, "3_5.0", "OP5.0C")
clean_cols <- str_replace(clean_cols, "1_8.0", "OP8.0A")
clean_cols <- str_replace(clean_cols, "2_8.0", "OP8.0B")
clean_cols <- str_replace(clean_cols, "3_8.0", "OP8.0C")

colnames(metabar_count_clean)[-1] <-  clean_cols
cols <- colnames(metabar_count_clean)

## Raw rarefaction curves
# -----------------------

# GP
gp1 <- cols[grepl("EN0.2A", cols)]
gp2 <- cols[grepl("EN0.2B", cols)]
gp3 <- cols[grepl("EN0.2C", cols)]
# A02
a02 <- cols[grepl("OP0.2A", cols)]
b02 <- cols[grepl("OP0.2B", cols)]
c02 <- cols[grepl("OP0.2C", cols)]
# A02
a12 <- cols[grepl("OP1.2A", cols)]
b12 <- cols[grepl("OP1.2B", cols)]
c12 <- cols[grepl("OP1.2C", cols)]
# A02
a50 <- cols[grepl("OP5.0A", cols)]
b50 <- cols[grepl("OP5.0B", cols)]
c50 <- cols[grepl("OP5.0C", cols)]
# A02
a80 <- cols[grepl("OP8.0A", cols)]
b80 <- cols[grepl("OP8.0B", cols)]
c80 <- cols[grepl("OP8.0C", cols)]

# Colors
colors <- c("#F8766D", "#B79F00", "#00BF7D", "#00B0F6", "#E76BF3")

# Raw curve
cc <- t(metabar_count_clean)
colnames(cc) <- cc[1,]
cc <- cc[-1,]
cc <- apply(cc, 2, as.numeric)
rownames(cc) <- colnames(metabar_count_clean)[-1]

# GP
raremax <- min(rowSums(cc[gp1,]))
png("GenomeDK_Metabar/raw_curves/raw_curve.replicates.GP1.png", width=2000, height=2000, res=300)
rarecurve(cc[gp1,], step = 20, sample = raremax, col = colors[1], cex = 0.6)
dev.off()

raremax <- min(rowSums(cc[gp2,]))
png("GenomeDK_Metabar/raw_curves/raw_curve.replicates.GP2.png", width=2000, height=2000, res=300)
rarecurve(cc[gp2,], step = 20, sample = raremax, col = colors[1], cex = 0.6)
dev.off()

raremax <- min(rowSums(cc[gp3,]))
png("GenomeDK_Metabar/raw_curves/raw_curve.replicates.GP3.png", width=2000, height=2000, res=300)
rarecurve(cc[gp3,], step = 20, sample = raremax, col = colors[1], cex = 0.6)
dev.off()

# 0.2
raremax <- min(rowSums(cc[a02,]))
png("GenomeDK_Metabar/raw_curves/raw_curve.replicates.A02.png", width=2000, height=2000, res=300)
rarecurve(cc[a02,], step = 20, sample = raremax, col = colors[2], cex = 0.6)
dev.off()

raremax <- min(rowSums(cc[b02,]))
png("GenomeDK_Metabar/raw_curves/raw_curve.replicates.B02.png", width=2000, height=2000, res=300)
rarecurve(cc[b02,], step = 20, sample = raremax, col = colors[2], cex = 0.6)
dev.off()

raremax <- min(rowSums(cc[c02,]))
png("GenomeDK_Metabar/raw_curves/raw_curve.replicates.C02.png", width=2000, height=2000, res=300)
rarecurve(cc[c02,], step = 20, sample = raremax, col = colors[2], cex = 0.6)
dev.off()

# 1.2
raremax <- min(rowSums(cc[a12,]))
png("GenomeDK_Metabar/raw_curves/raw_curve.replicates.A12.png", width=2000, height=2000, res=300)
rarecurve(cc[a12,], step = 20, sample = raremax, col = colors[3], cex = 0.6)
dev.off()

raremax <- min(rowSums(cc[b12,]))
png("GenomeDK_Metabar/raw_curves/raw_curve.replicates.B12.png", width=2000, height=2000, res=300)
rarecurve(cc[b12,], step = 20, sample = raremax, col = colors[3], cex = 0.6)
dev.off()

raremax <- min(rowSums(cc[c12,]))
png("GenomeDK_Metabar/raw_curves/raw_curve.replicates.C12.png", width=2000, height=2000, res=300)
rarecurve(cc[c12,], step = 20, sample = raremax, col = colors[3], cex = 0.6)
dev.off()

# 5.0
raremax <- min(rowSums(cc[a50,]))
png("GenomeDK_Metabar/raw_curves/raw_curve.replicates.A50.png", width=2000, height=2000, res=300)
rarecurve(cc[a50,], step = 20, sample = raremax, col = colors[4], cex = 0.6)
dev.off()

raremax <- min(rowSums(cc[b50,]))
png("GenomeDK_Metabar/raw_curves/raw_curve.replicates.B50.png", width=2000, height=2000, res=300)
rarecurve(cc[b50,], step = 20, sample = raremax, col = colors[4], cex = 0.6)
dev.off()

raremax <- min(rowSums(cc[c50,]))
png("GenomeDK_Metabar/raw_curves/raw_curve.replicates.C50.png", width=2000, height=2000, res=300)
rarecurve(cc[c50,], step = 20, sample = raremax, col = colors[4], cex = 0.6)
dev.off()

# 8.0
raremax <- min(rowSums(cc[a80,]))
png("GenomeDK_Metabar/raw_curves/raw_curve.replicates.A80.png", width=2000, height=2000, res=300)
rarecurve(cc[a80,], step = 20, sample = raremax, col = colors[5], cex = 0.6)
dev.off()

raremax <- min(rowSums(cc[b80,]))
png("GenomeDK_Metabar/raw_curves/raw_curve.replicates.B80.png", width=2000, height=2000, res=300)
rarecurve(cc[b80,], step = 20, sample = raremax, col = colors[5], cex = 0.6)
dev.off()

raremax <- min(rowSums(cc[c80,]))
png("GenomeDK_Metabar/raw_curves/raw_curve.replicates.C80.png", width=2000, height=2000, res=300)
rarecurve(cc[c80,], step = 20, sample = raremax, col = colors[5], cex = 0.6)
dev.off()

## Accumulation plots
# -------------------



## Rarefaction replicates
# -----------------------

# Function: rarefy replicates
rarefy_reps <- function(df, columns, name) {
    cs2mr <- df[,c(columns,"id")]
    colnames(cs2mr)[-ncol(cs2mr)] <- paste("sample:", colnames(cs2mr)[-ncol(cs2mr)], sep="")
    # Save temportal file
    f = paste('GenomeDK_Metabar/temporal/', name, '_for_rep_rarefy.txt', sep="")
    write.table(cs2mr, file=f, quote=FALSE, sep='\t', col.names = NA, row.names=TRUE)
    # Load temporal file into ROBITools
    dfimp <- import.metabarcoding.data(f)
    # Compute median
    mdn <- median(rowSums(dfimp@reads))
    # Median Rarefy
    dfrar <- ROBITools::rarefy(dfimp, n = mdn, MARGIN="sample")
    # Return
    return(dfrar)
}

# Rarefied curves
# GP1
rc_gp1 <- rarefy_reps(metabar_count_clean, gp1, "EN0.2A")
png("GenomeDK_Metabar/rarefied_curves/rarefied_curve.replicates.GP1.png", width=2000, height=2000, res=300)
rarecurve(rc_gp1@reads, step = 20, col = colors[1], cex = 0.6)
dev.off()

# GP2
rc_gp2 <- rarefy_reps(metabar_count_clean, gp2, "EN0.2B")
png("GenomeDK_Metabar/rarefied_curves/rarefied_curve.replicates.GP2.png", width=2000, height=2000, res=300)
rarecurve(rc_gp2@reads, step = 20, col = colors[1], cex = 0.6)
dev.off()

# GP3
rc_gp3 <- rarefy_reps(metabar_count_clean, gp3, "EN0.2C")
png("GenomeDK_Metabar/rarefied_curves/rarefied_curve.replicates.GP3.png", width=2000, height=2000, res=300)
rarecurve(rc_gp3@reads, step = 20, col = colors[1], cex = 0.6)
dev.off()

# A02
rc_a02 <- rarefy_reps(metabar_count_clean, a02, "OP0.2A")
png("GenomeDK_Metabar/rarefied_curves/rarefied_curve.replicates.A02.png", width=2000, height=2000, res=300)
rarecurve(rc_a02@reads, step = 20, col = colors[2], cex = 0.6)
dev.off()

# B02
rc_b02 <- rarefy_reps(metabar_count_clean, b02, "OP0.2B")
png("GenomeDK_Metabar/rarefied_curves/rarefied_curve.replicates.B02.png", width=2000, height=2000, res=300)
rarecurve(rc_b02@reads, step = 20, col = colors[2], cex = 0.6)
dev.off()

# C02
rc_c02 <- rarefy_reps(metabar_count_clean, c02, "OP0.2C")
png("GenomeDK_Metabar/rarefied_curves/rarefied_curve.replicates.C02.png", width=2000, height=2000, res=300)
rarecurve(rc_c02@reads, step = 20, col = colors[2], cex = 0.6)
dev.off()

# A12
rc_a12 <- rarefy_reps(metabar_count_clean, a12, "OP1.2A")
png("GenomeDK_Metabar/rarefied_curves/rarefied_curve.replicates.A12.png", width=2000, height=2000, res=300)
rarecurve(rc_a12@reads, step = 20, col = colors[3], cex = 0.6)
dev.off()

# B12
rc_b12 <- rarefy_reps(metabar_count_clean, b12, "OP1.2B")
png("GenomeDK_Metabar/rarefied_curves/rarefied_curve.replicates.B12.png", width=2000, height=2000, res=300)
rarecurve(rc_b12@reads, step = 20, col = colors[3], cex = 0.6)
dev.off()

# C12
rc_c12 <- rarefy_reps(metabar_count_clean, c12, "OP1.2C")
png("GenomeDK_Metabar/rarefied_curves/rarefied_curve.replicates.C12.png", width=2000, height=2000, res=300)
rarecurve(rc_c12@reads, step = 20, col = colors[3], cex = 0.6)
dev.off()

# A50
rc_a50 <- rarefy_reps(metabar_count_clean, a50, "OP5.0A")
png("GenomeDK_Metabar/rarefied_curves/rarefied_curve.replicates.A50.png", width=2000, height=2000, res=300)
rarecurve(rc_a50@reads, step = 20, col = colors[4], cex = 0.6)
dev.off()

# B50
rc_b50 <- rarefy_reps(metabar_count_clean, b50, "OP5.0B")
png("GenomeDK_Metabar/rarefied_curves/rarefied_curve.replicates.B50.png", width=2000, height=2000, res=300)
rarecurve(rc_b50@reads, step = 20, col = colors[4], cex = 0.6)
dev.off()

# C50
rc_c50 <- rarefy_reps(metabar_count_clean, c50, "OP5.0C")
png("GenomeDK_Metabar/rarefied_curves/rarefied_curve.replicates.C50.png", width=2000, height=2000, res=300)
rarecurve(rc_c50@reads, step = 20, col = colors[4], cex = 0.6)
dev.off()

# A80
rc_a80 <- rarefy_reps(metabar_count_clean, a80, "OP8.0A")
png("GenomeDK_Metabar/rarefied_curves/rarefied_curve.replicates.A80.png", width=2000, height=2000, res=300)
rarecurve(rc_a80@reads, step = 20, col = colors[5], cex = 0.6)
dev.off()

# B80
rc_b80 <- rarefy_reps(metabar_count_clean, b80, "OP8.0B")
png("GenomeDK_Metabar/rarefied_curves/rarefied_curve.replicates.B80.png", width=2000, height=2000, res=300)
rarecurve(rc_b80@reads, step = 20, col = colors[5], cex = 0.6)
dev.off()

# C80
rc_c80 <- rarefy_reps(metabar_count_clean, c80, "OP8.0C")
png("GenomeDK_Metabar/rarefied_curves/rarefied_curve.replicates.C80.png", width=2000, height=2000, res=300)
rarecurve(rc_c80@reads, step = 20, col = colors[5], cex = 0.6)
dev.off()

## Rarefaction samples
# --------------------

# Function: Aggregate replicates
collapse_reps <- function(reads) {
    # Name and median
    sample <- unique(str_replace(rownames(reads), "_[0-9]", ""))
    values <- colSums(rc_l[[1]]@reads)
    # Table
    df <- t(as.data.frame(values))
    rownames(df) <- sample
    # Return
    return(df)
}

# Create table
rc_l <- list(rc_gp1, rc_gp2, rc_gp3,
     rc_a02, rc_b02, rc_c02,
     rc_a12, rc_b12, rc_c12,
     rc_a50, rc_b50, rc_c50,
     rc_a80, rc_b80, rc_c80)

rc <- lapply(rc_l, function(x) collapse_reps(x@reads))
rc <- do.call("rbind", rc)

## Shotgun
# --------

# --------
## Shotgun
# --------
shotgun_taxa_count <- read.delim("GenomeDK_LCA/counts.lca.rarefy.tsv")

# Clean shotgun
#shotgun <- shotgun_taxa_count %>% filter(!species %in% contaminants$species)

# Select Eukaryota
shotgun_eu <- shotgun_taxa_count %>% filter(superkingdom == "Eukaryota")

sc <- shotgun_taxa_count[,c(9:23)]
sc <- t(sc)
colnames(sc) <- paste("tax", 1:ncol(sc), sep="")

col <- c(rep("#F8766D",3), rep("#B79F00",3), rep("#00BF7D",3), rep("#00B0F6",3), rep("#E76BF3",3))

png("GenomeDK_LCA/rarefied_curve.all.png", width=2000, height=2000, res=300)
rarecurve(sc, step = 20, col = col, cex = 0.6)
dev.off()

## Comparison
# -----------

# Function
display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

# Phylum
phyl_metabar <- unique(metabarcoding$phylum[!is.na(metabarcoding$phylum)])
phyl_shotgun <- unique(shotgun_eu$phylum[!is.na(shotgun_eu$phylum)])
taxa <- list(metabarcoding = phyl_metabar, shotgun_eu = phyl_shotgun)

# Venn Diagram (Saved manually)
display_venn(taxa,
  category.names = c("Shotgun", "Metabarcoding"),
  fill = c("royalblue", "gold"),
  cex = 3,
  cat.cex = 2.5, cat.fontface = "bold", 
  cat.pos = c(10,-10),
  cat.dist = c(0.05, 0.05),
  scaled = F)

# Class
cla_metabar <- unique(metabarcoding$class[!is.na(metabarcoding$class)])
cla_shotgun <- unique(shotgun_eu$class[!is.na(shotgun_eu$class)])
taxa <- list(metabarcoding = cla_metabar, shotgun_eu = cla_shotgun)

# Venn Diagram (Saved manually)
display_venn(taxa,
  category.names = c("Shotgun", "Metabarcoding"),
  fill = c("royalblue", "gold"),
  cex = 3,
  cat.cex = 2.5, cat.fontface = "bold",
  cat.pos = c(-10,10),
  cat.dist = c(0.05, 0.05),
  scaled = F)

# Order
ord_metabar <- unique(metabarcoding$order[!is.na(metabarcoding$order)])
ord_shotgun <- unique(shotgun_eu$order[!is.na(shotgun_eu$order)])
taxa <- list(metabarcoding = ord_metabar, shotgun_eu = ord_shotgun)

# Venn Diagram (Saved manually)
display_venn(taxa,
  category.names = c("Shotgun", "Metabarcoding"),
  fill = c("royalblue", "gold"),
  cex = 3,
  cat.cex = 2.5, cat.fontface = "bold",
  cat.pos = c(-10,10),
  cat.dist = c(0.05, 0.05),
  scaled = F)

# Family
fam_metabar <- unique(metabarcoding$family[!is.na(metabarcoding$family)])
fam_shotgun <- unique(shotgun_eu$family[!is.na(shotgun_eu$family)])
taxa <- list(metabarcoding = fam_metabar, shotgun_eu = fam_shotgun)

# Venn Diagram (Saved manually)
display_venn(taxa,
  category.names = c("Shotgun", "Metabarcoding"),
  fill = c("royalblue", "gold"),
  cex = 3,
  cat.cex = 2.5, cat.fontface = "bold",
  cat.pos = c(-10,10),
  cat.dist = c(0.05, 0.05),
  scaled = F)

# Genus
gen_metabar <- unique(metabarcoding$genus[!is.na(metabarcoding$genus)])
gen_shotgun <- unique(shotgun_eu$genus[!is.na(shotgun_eu$genus)])
taxa <- list(metabarcoding = gen_metabar, shotgun_eu = gen_shotgun)

# Venn Diagram (Saved manually)
display_venn(taxa,
  category.names = c("Shotgun", "Metabarcoding"),
  fill = c("royalblue", "gold"),
  cex = 3,
  cat.cex = 2.5, cat.fontface = "bold",
  cat.pos = c(-10,10),
  cat.dist = c(0.05, 0.05),
  scaled = F)

# Phylum
phyl_metabar <- unique(metabarcoding$phylum[!is.na(metabarcoding$phylum)])
phyl_shotgun <- unique(shotgun_eu$phylum[!is.na(shotgun_eu$phylum)])
taxa <- list(metabarcoding = phyl_metabar, shotgun_eu = phyl_shotgun)

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

# Presence absence table
rownames(df) <- ifelse(!is.na(df$metabarcoding), df$metabarcoding, df$shotgun )
df$metabarcoding <- ifelse(!is.na(df$metabarcoding), 1, 0)
df$shotgun <- ifelse(!is.na(df$shotgun), 1, 0)

# Recover phylum
df$phylum <- rownames(df)

# Order columns
df <- df[,c(3,1,2)]

# Share
df <- df %>% mutate(belonging = ifelse(metabarcoding == 1 & shotgun == 1, "shared",
  ifelse(metabarcoding == 1 & shotgun == 0, "metabarcoding","shotgun")))

write.table(df, "phylum_comparison_table.eukaryota.tsv", quote = F, row.names = F, sep="\t")

# Heatmap
dfm <- melt(df)
dfm$variable <- ifelse(dfm$variable == "shotgun", "Shotgun", "Metabarcoding")
dfm$value <- ifelse(dfm$value == 0, "Absence", "Presence")
dfm$phylum <- factor(dfm$phylum, levels = rev(sort(unique(dfm$phylum))))
p <- dfm %>% ggplot() + geom_tile(aes(x=variable, y=phylum), fill="gray60", color="black", alpha=0.8) +
# Shotgun
geom_tile(data=subset(dfm, belonging == "shotgun" & variable == "Shotgun"),
          aes(x=variable, y=phylum), fill="royalblue", color="black") +
geom_tile(data=subset(dfm, belonging == "shotgun" & variable == "Metabarcoding"),
          aes(x=variable, y=phylum), fill="gray90", color="white") +
# Metabarcoding
geom_tile(data=subset(dfm, belonging == "metabarcoding" & variable == "Metabarcoding"),
          aes(x=variable, y=phylum), fill="gold", color="black") +
geom_tile(data=subset(dfm, belonging == "metabarcoding" & variable == "Shotgun"),
          aes(x=variable, y=phylum), fill="gray90", color="white") +
labs(x = "", y = "Phylum", fill="") +
theme_classic() %+replace% theme(axis.text.x = element_text(size=18),
  axis.text.y = element_text(size=11),
  axis.title = element_text(size=20),
  legend.title = element_text(size=18, face="bold"),
  legend.text = element_text(size=16))

png(file="comp.Presence_Absence.phylum.png", width=2300, height=2500, res=300)
p
dev.off()

## Validation
# -----------

arter <- read.table("arter.tsv", sep="\t", header=T)
colnames(arter) <- c("Phylum","Class","Order","Family","Genus","Species", "Marine")
marine <- arter[arter$Marine == "Yes",]

table(arter$Marine)

# Shotgun
ps <- shotgun_taxa_count %>% filter(superkingdom == "Eukaryota", kingdom == "Metazoa")
ps <- ps[,c(3,9:23)]
ps <- melt(ps) %>% group_by(phylum,variable) %>% summarise(value = sum(value)) %>%
mutate(origin = ifelse(phylum %in% marine$Phylum, "Danish", "Exotic"))

unique(ps$phylum[ps$phylum %in% marine$Phylum])
unique(ps$phylum[!ps$phylum %in% marine$Phylum])

newnames <- c('EN0.2A','EN0.2B','EN0.2C','OP0.2A','OP0.2B','OP0.2C','OP1.2A','OP1.2B','OP1.2C','OP5.0A','OP5.0B','OP5.0C','OP8.0A','OP8.0B','OP8.0C')
ps$variable <- factor(ps$variable, levels = newnames)
ps <- ps %>% mutate(type = ifelse(grepl("EN", variable), "Enclosed", "Open"))

p <- ps %>% drop_na() %>%
ggplot() + geom_col(aes(x=variable, y=value, fill=origin), color="black") +
scale_fill_manual(values = c("forestgreen", "firebrick")) +
labs(fill = "", x="", y="#Reads per phyla", title="Metazoa") +
facet_wrap(~type, scales="free_x") +
theme_classic() %+replace% theme(
  title = element_text(size=20, face="bold"),
  axis.text.x = element_text(angle = 90, size=18),
  axis.text.y = element_text(size=18),
  axis.title = element_text(size=20),
  strip.text.x = element_text(size = 16, face="bold"),
  strip.background = element_blank(),
  legend.title = element_text(size=18, face="bold"),
  legend.text = element_text(size=16))

gt = ggplot_gtable(ggplot_build(p))
gt$widths[5] = 0.28*gt$widths[5]

#png(file="07-Plots/bars.spkm.raw.png", width=2300, height=2000, res=300)
grid.draw(gt)
#dev.off()

# Plot
ps %>% drop_na() %>% group_by(variable, origin) %>% summarise(count = n()) %>% ungroup() %>%
filter(variable == "OP0.2A") %>% 
ggplot() + geom_col(aes(x=origin, y=count, fill=origin)) +
geom_text(aes(x=origin, y=count+1, label=count), size=6) +
scale_fill_manual(values = c("forestgreen", "firebrick")) +
labs(fill = "", x="", y="#Phylum") +
theme_classic() %+replace% theme(axis.text.x = element_text(size=18),
  axis.text.y = element_text(size=18),
  axis.title = element_text(size=20),
  legend.text = element_text(size=16))






