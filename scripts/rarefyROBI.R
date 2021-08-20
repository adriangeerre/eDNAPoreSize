# Libraries
library(ROBITools)
library(tidyverse)
library(reshape2)

# Data
dfraw <- read.delim("raw.taxa.count.tsv")

# Plots
p <- dfraw %>% filter(superkingdom %in% c("Bacteria","Archaea","Eukaryota","Viruses")) %>% select(c(1,9:23)) %>% melt() %>% group_by(superkingdom, variable) %>% summarise(count = sum(value)) %>% ungroup() %>% group_by(variable) %>% mutate(perc = count/sum(count) * 100)

p$superkingdom <- factor(p$superkingdom, levels = c("Archaea","Viruses","Bacteria","Eukaryota"))
p$variable <- factor(p$variable, levels = c('GP_1','GP_2','GP_3','A0_2','B0_2','C0_2','A1_2','B1_2','C1_2','A5_0','B5_0','C5_0','A8_0','B8_0','C8_0')) 

p %>% ggplot() + geom_col(aes(x=variable,y=perc,fill=superkingdom))

# PCA
pcaplot <- function(x) {
  mat  = t(x)
  pca  = prcomp(mat, scale=F)
  summary(pca)
  col = c(rep("blue3",3), rep("red3",3), rep("green3",3), rep("orange3",3),rep("darkcyan",3))
  col.leg = c("blue3", "red3", "green3", "orange3","darkcyan")
  nam = c("A0_2","B0_2","C0_2","A1_2","B1_2","C1_2","A5_0","B5_0","C5_0","A8_0","B8_0","C8_0","GP_1","GP_2","GP_3")
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
pcaplot(dfpca)

#-------------------------------------------------------------------#

# Prepare data for ROBITools
cols <- c("sample:A02","sample:B02","sample:C02","sample:A12","sample:B12","sample:C12","sample:A50","sample:B50","sample:C50","sample:A80","sample:B80","sample:C80","sample:GP1","sample:GP2","sample:GP3","count","superkingdom_name","kingdom_name","phylum_name","class_name","order_name","family_name","genus_name","species_name","id")

dfraw$count <- rowSums(dfraw[,c(9:23)])
dfraw$id <- paste("tax",1:dim(dfraw)[1],sep="")
dfraw <- dfraw[,c(9:24,1:8,25)]
colnames(dfraw) <- cols
write.table(dfraw, file='Pore_size_table_for_rarefy.txt', quote=FALSE, sep='\t', col.names = NA,row.names=TRUE) # Write table to read as metabarcoding data


# Normalization (Rarefy)
# Function rarefy gives the expected species richness in random subsamples of size sample from the community. The function rarefy is based on Hurlbert's (1971) formulation, and the standard errors on Heck et al. (1975). In order to rarefy the data we need to create a count of species (rows) per sample (columns). The workflow is based in Eva's script including Mads modifications (Vandloeb_Summer_20200529.Rmd).
dfimp <- import.metabarcoding.data("Pore_size_table_for_rarefy.txt") # Read data with ROBITools

## Rarefy by sample
mdn <- summary(rowSums(dfimp@reads[dfimp@samples$sample,]))[[3]]
dfrar <- ROBITools::rarefy(dfimp, n = mdn, MARGIN="sample")

## Update TaxaCounts
dfrar@motus$count <- colSums(dfrar@reads) # Add count to rarefy species

## Clean data after rarefy
table(colSums(dfrar@reads)>0) # 1445 taxa equal to 0
dfrar <- dfrar[,colSums(dfrar@reads)>0] # Remove species with no hits
dfrar <- dfrar[rowSums(dfrar@reads)>0,] # Remove samples with no hits

## Transform data to final dataframe
dfend <- dfrar@reads %>% t() %>% as.data.frame() %>% mutate(id = colnames(dfrar@reads))
dfend <- merge(dfend,dfrar@motus, by="id") %>% select(-c(Var.17,id,count))
dfend <- dfend[,c(16:23,1:15)]
colnames(dfend) <- c("superkingdom","kingdom","phylum","class","order","family","genus","species","A0_2","B0_2","C0_2","A1_2","B1_2","C1_2","A5_0","B5_0","C5_0","A8_0","B8_0","C8_0","GP_1","GP_2","GP_3")

write.table(dfend, file='counts.rarefy.tsv', quote=FALSE, sep='\t', col.names=T ,row.names=F)

# Plots
p <- dfend %>% filter(superkingdom %in% c("Bacteria","Archaea","Eukaryota","Viruses")) %>% select(c(1,9:23)) %>% melt() %>% group_by(superkingdom, variable) %>% summarise(count = sum(value)) %>% ungroup() %>% group_by(variable) %>% mutate(perc = count/sum(count) * 100)

p$superkingdom <- factor(p$superkingdom, levels = c("Archaea","Viruses","Bacteria","Eukaryota"))
p$variable <- factor(p$variable, levels = c('GP_1','GP_2','GP_3','A0_2','B0_2','C0_2','A1_2','B1_2','C1_2','A5_0','B5_0','C5_0','A8_0','B8_0','C8_0')) 

p %>% ggplot() + geom_col(aes(x=variable,y=perc,fill=superkingdom))

# PCA
dfpca <- dfend[9:23]
dfpca[dfpca > 0] <- 1
pcaplot(dfpca)

