# Libraries
library(ROBITools)
library(tidyverse)
library(reshape2)
library(vegan)
library(ggalt)
library(ggpubr)
library(ggrepel)

# Data (Eukaryota less than expected because multiple taxids refer to the same taxa)
dfraw <- read.delim("lca.taxa.count.tsv", check.names=F)
colnames(dfraw) <- str_replace(colnames(dfraw),"06-mergeLCA//","")

# Rename Smaples
oldnames <- c('GP_1','GP_2','GP_3','A0_2','B0_2','C0_2','A1_2','B1_2','C1_2','A5_0','B5_0','C5_0','A8_0','B8_0','C8_0')
newnames <- c('EN0.2A','EN0.2B','EN0.2C','OP0.2A','OP0.2B','OP0.2C','OP1.2A','OP1.2B','OP1.2C','OP5.0A','OP5.0B','OP5.0C','OP8.0A','OP8.0B','OP8.0C')
names(newnames) <- oldnames

dfraw <- dfraw %>% rename_at(vars(oldnames), ~newnames)

## Plots Superkingdom
spkm_bars <- dfraw %>% filter(superkingdom %in% c("Bacteria","Archaea","Eukaryota","Viruses")) %>% select(c(1,9:23)) %>% melt() %>% group_by(superkingdom, variable) %>% summarise(count = sum(value)) %>% ungroup() %>% group_by(variable) %>% mutate(perc = count/sum(count) * 100)

spkm_bars$superkingdom <- factor(spkm_bars$superkingdom, levels = c("Archaea","Viruses","Bacteria","Eukaryota"))
spkm_bars$variable <- factor(spkm_bars$variable, levels = newnames)

p <- spkm_bars %>% ggplot() + geom_col(aes(x=variable,y=perc,fill=superkingdom)) +
 labs(x = "Sample", y = "Percentage", fill="Domain") +
 scale_fill_manual(values = c("coral3", "darkolivegreen3", "darkcyan", "orange3")) +
 theme_classic() %+replace% theme(axis.text.x = element_text(angle = 45))

ggsave("07-Plots/bars.spkm.raw.png", p, units="px", width=3500)

## Plots Phylum
phyl_bars <- dfraw %>% filter(superkingdom %in% c("Bacteria","Archaea","Eukaryota","Viruses")) %>% select(c(1,3,9:23)) %>% melt() %>% drop_na() %>% group_by(superkingdom, phylum, variable) %>% summarise(count = sum(value)) %>% ungroup() %>% group_by(variable) %>% mutate(perc = count/sum(count) * 100) %>% filter(count > 0)

phyl_bars$superkingdom <- factor(phyl_bars$superkingdom, levels = c("Archaea","Viruses","Bacteria","Eukaryota"))
phyl_bars$variable <- factor(phyl_bars$variable, levels = newnames)

p <- phyl_bars %>% ggplot() + geom_col(aes(x=variable,y=perc,fill=superkingdom), color="black") +
 labs(x = "Sample", y = "Percentage", fill="Domain") +
 guides(fill = "none") +
 scale_fill_manual(values = c("coral3", "darkolivegreen3", "darkcyan", "orange3")) +
 theme_classic() %+replace% theme(axis.text.x = element_text(angle = 45)) 

ggsave("07-Plots/bars.phyl.raw.png", p, units="px", width=3500)

# PCA
pcaplot <- function(x) {
  mat  = t(x)
  pca  = prcomp(mat, scale=F)
  summary(pca)
  col = c(rep("blue3",3), rep("red3",3), rep("green3",3), rep("orange3",3),rep("darkcyan",3))
  col.leg = c("blue3", "red3", "green3", "orange3","darkcyan")
  nam = newnames
  legend = c("0.2","1.2","5.0","8.0")
  aval  = pca$sdev^2
  cvar  = pca$rotation %*% diag(pca$sdev)
  #screeplot(pca, npcs=min(10, length(pca$sdev)), type="lines", col="red", main="Screeplot")
  pos = c(rep(1,length(nam)-1), 3)
  xs  = pca$x[,1]
  ys  = pca$x[,2]
  
  xlab = paste("PC1 (", summary(pca)$importance[2,1]*100, " %)", sep='')
  ylab = paste("PC2 (", summary(pca)$importance[2,2]*100, " %)", sep='')
  plot(pca$x[,1:2], pch=19, col=col, main="Two-dimensional PCA (Presence/Absence)", xlab=xlab, ylab=ylab)
  abline(h=0, v=0, lty=2, col=8)
  text(xs, ys, labels=nam, pos=pos)
}

dfpca <- dfraw[9:23]
dfpca[dfpca > 0] <- 1
pca <- pcaplot(dfpca)

# ggsave("07-Plots/PCA.raw.png", pca, units="px", width=3500)

# Relative Abundance PCA (nMDS) - Mads!

nMDS_plot <- function(reads, samples) {
  # Create scores
  nmds = metaMDS(reads, distance = "bray")
  nmds.sites.scores = as.data.frame(vegan::scores(nmds)$sites)
  #nmds.species.scores = as.data.frame(vegan::scores(nmds)$species)

  # Add metadata
  nmds.sites.scores$sample = samples$sample
  nmds.sites.scores$pore = samples$pore

  # Factorize
  nmds.sites.scores <- nmds.sites.scores %>%
    dplyr::mutate(pore=factor(pore, levels=c("GP", "0.2", "1.2","5.0", "8.0")))

  # Plot
  p <- ggplot(nmds.sites.scores, aes(x = NMDS1, y = NMDS2, label=sample)) + 
    geom_point(size = 4, aes(colour = pore)) +
    geom_encircle(aes(fill = pore), s_shape = 1, expand = 0,
                  alpha = 0.2, color = "black", show.legend = FALSE) +
    geom_text_repel() +
    labs(colour='Filter type') +
    theme_classic()
  
  return(p)
}

#-------------------------------------------------------------------#

# Prepare data for ROBITools
cols <- c("sample:OP0.2A","sample:OP0.2B","sample:OP0.2C","sample:OP1.2A","sample:OP1.2B","sample:OP1.2C","sample:OP5.0A","sample:OP5.0B","sample:OP5.0C","sample:OP8.0A","sample:OP8.0B","sample:OP8.0C","sample:EN0.2A","sample:EN0.2B","sample:EN0.2C","count","superkingdom_name","kingdom_name","phylum_name","class_name","order_name","family_name","genus_name","species_name","id")

dfraw$count <- rowSums(dfraw[,c(9:23)])
dfraw$id <- paste("tax",1:dim(dfraw)[1],sep="")
dfraw <- dfraw[,c(9:24,1:8,25)]
colnames(dfraw) <- cols
write.table(dfraw, file='Pore_size_table_for_rarefy.txt', quote=FALSE, sep='\t', col.names = NA,row.names=TRUE) # Write table to read as metabarcoding data


# Normalization (Rarefy)
# Function rarefy gives the expected species richness in random subsamples of size sample from the community. The function rarefy is based on Hurlbert's (1971) formulation, and the standard errors on Heck et al. (1975). In order to rarefy the data we need to create a count of species (rows) per sample (columns). The workflow is based in Eva's script including Mads modifications (Vandloeb_Summer_20200529.Rmd).
dfimp <- import.metabarcoding.data("Pore_size_table_for_rarefy.txt") # Read data with ROBITools

## nMDS plot (before rarefaction)
samples <- dfimp@samples
samples <- samples %>% mutate(pore = ifelse(sample %in% c("OP0.2A","OP0.2B","OP0.2C"), "0.2", ifelse(sample %in% c("OP1.2A","OP1.2B","OP1.2C"), "1.2", ifelse(sample %in% c("OP5.0A","OP5.0B","OP5.0C"), "5.0", ifelse(sample %in% c("OP8.0A","OP8.0B","OP8.0C"), "8.0", "GP")))))
reads <- dfimp@reads
p <- nMDS_plot(reads, samples)

ggsave("07-Plots/nMDS.raw.png", p, units="px", width=3500)

## Rarefy by sample
mdn <- summary(rowSums(dfimp@reads[dfimp@samples$sample,]))[[3]]
dfrar <- ROBITools::rarefy(dfimp, n = mdn, MARGIN="sample")

## Update TaxaCounts
dfrar@motus$count <- colSums(dfrar@reads) # Add count to rarefy species

## Clean data after rarefy
table(colSums(dfrar@reads)>0) # 2040 taxa equal to 0
dfrar <- dfrar[,colSums(dfrar@reads)>0] # Remove species with no hits
dfrar <- dfrar[rowSums(dfrar@reads)>0,] # Remove samples with no hits

## Transform data to final dataframe
dfend <- dfrar@reads %>% t() %>% as.data.frame() %>% mutate(id = colnames(dfrar@reads))
dfend <- merge(dfend,dfrar@motus, by="id") %>% select(-c(Var.17,id,count))
dfend <- dfend[,c(16:23,1:15)]
colnames(dfend)[1:8] <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")

write.table(dfend, file='counts.lca.rarefy.tsv', quote=FALSE, sep='\t', col.names=T ,row.names=F)

# Plots
spkm_bars <- dfend %>% filter(superkingdom %in% c("Bacteria","Archaea","Eukaryota","Viruses")) %>% select(c(1,9:23)) %>% melt() %>% group_by(superkingdom, variable) %>% summarise(count = sum(value)) %>% ungroup() %>% group_by(variable) %>% mutate(perc = count/sum(count) * 100)

spkm_bars$superkingdom <- factor(spkm_bars$superkingdom, levels = c("Archaea","Viruses","Bacteria","Eukaryota"))
spkm_bars$variable <- factor(spkm_bars$variable, levels = newnames) 

p <- spkm_bars %>% ggplot() + geom_col(aes(x=variable,y=perc,fill=superkingdom)) +
 labs(x = "Sample", y = "Percentage", fill="Domain") +
 scale_fill_manual(values = c("coral3", "darkolivegreen3", "darkcyan", "orange3")) +
 theme_classic() %+replace% theme(axis.text.x = element_text(angle = 45))

ggsave("07-Plots/bars.spkm.rarefy.png", p, units="px", width=3500)

# PCA
# dfpca <- dfend[9:23]
# dfpca[dfpca > 0] <- 1
# pcaplot(dfpca)

## nMDS plot (after rarefaction)
samples <- dfrar@samples
samples <- samples %>% mutate(pore = ifelse(sample %in% c("OP0.2A","OP0.2B","OP0.2C"), "0.2", ifelse(sample %in% c("OP1.2A","OP1.2B","OP1.2C"), "1.2", ifelse(sample %in% c("OP5.0A","OP5.0B","OP5.0C"), "5.0", ifelse(sample %in% c("OP8.0A","OP8.0B","OP8.0C"), "8.0", "GP")))))
reads <- dfrar@reads
p <- nMDS_plot(reads, samples)

ggsave("07-Plots/nMDS.rarefy.png", p, units="px", width=3500)