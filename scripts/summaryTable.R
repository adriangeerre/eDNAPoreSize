# Libraries
library(tidyverse)
library(reshape2)
library(ComplexHeatmap)

# Data
df <- read.delim("counts.rarefy.tsv")
s <- c('A0_2','B0_2','C0_2','A1_2','B1_2','C1_2','A5_0','B5_0','C5_0','A8_0','B8_0','C8_0','GP_1','GP_2','GP_3')

# Phylum
dfeu <- df %>% filter(superkingdom == "Eukaryota")
ph <- dfeu %>% group_by(phylum) %>% drop_na("phylum") %>% select(s) %>% melt() %>% ungroup() %>% group_by(phylum, variable) %>% summarise(value = sum(value)) %>% spread(phylum, value) %>% rename(sample = variable) %>% as.data.frame()
rownames(ph) <- ph$sample
ph <- ph %>% select(-sample) %>% t() %>% as.data.frame() %>% arrange(desc(A0_2))
write.table(ph, file='reads.phylum.tsv', quote=FALSE, sep='\t', col.names=T ,row.names=F)

phperc <- apply(ph,2,function(x){round(x/sum(x)*100,2)})
write.table(phperc, file='reads.perc.phylum.tsv', quote=FALSE, sep='\t', col.names=T ,row.names=F)

# Sort
labels <- (df[,c(1,2,3)] %>% filter(superkingdom == "Eukaryota") %>% unique() %>% select(kingdom, phylum))
labels <- merge(labels, rownames(ph), by.x="phylum", by.y="y") %>% arrange(kingdom)
labels$kingdom[is.na(labels$kingdom)] <- "Undefined"

# Color
lab <- colorRampPalette(c("lightblue","lightyellow","orange","red","darkred"),  space = "Lab")

# Scale Rows
scale_heatmap_row <- function(x) {
  x <- sweep(x, 1L, rowMeans(x, na.rm = F), check.margin = FALSE)
  sx <- apply(x, 1L, sd, na.rm = F)
  x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
  return(x)
}

phs <- scale_heatmap_row(ph)
phs <- phs[labels$phylum,]

# Reorder GPs
phs <- phs[,c(13:15,1:12)]

# Heatmap
annot <- rowAnnotation(Kingdom = labels$kingdom, show_annotation_name = F, col = list(Kingdom = c("Fungi"="coral3", "Metazoa"="darkcyan", "Viridiplantae"="darkolivegreen3", "Undefined"="rosybrown3")), border = T)
Heatmap(as.matrix(phs), name = "Row Scaled \n Counts", cluster_rows = F, cluster_columns = F, row_order = labels$phylum, col = lab(16), left_annotation = annot, row_split = labels$kingdom, row_names_side = "left", row_names_gp = gpar(fontsize = 10), border = "gray60", column_split = c(rep("Enclosed",3), rep("Open",12)))


# NON-EUKARYOTA

ph <- df %>% group_by(phylum) %>% drop_na("phylum") %>% select(s) %>% melt() %>% ungroup() %>% group_by(phylum, variable) %>% summarise(value = sum(value)) %>% spread(phylum, value) %>% rename(sample = variable) %>% as.data.frame()
rownames(ph) <- ph$sample
ph <- ph %>% select(-sample) %>% t() %>% as.data.frame() %>% arrange(desc(A0_2))

labels <- (df[,c(1,2,3)] %>% unique() %>% select(superkingdom, phylum))
labels <- merge(labels, rownames(ph), by.x="phylum", by.y="y") %>% arrange(superkingdom)

phs <- scale_heatmap_row(ph)
phs <- phs[labels$phylum,]
phs <- phs[,c(13:15,1:12)]

annot <- rowAnnotation(Superkingdom = labels$superkingdom, show_annotation_name = F, col = list(Superkingdom = c("Archaea"="coral3", "Bacteria"="darkcyan", "Eukaryota"="darkolivegreen3", "Viruses"="rosybrown3")), border = T)
p <- Heatmap(as.matrix(phs), name = "Row Scaled \n Counts", cluster_rows = F, cluster_columns = F, row_order = labels$phylum, col = lab(16), left_annotation = annot, row_split = labels$superkingdom, row_names_side = "left", row_names_gp = gpar(fontsize = 4), border = "gray60", column_split = c(rep("Enclosed",3), rep("Open",12)))

png("heatmap.full.phylum.png", width = 3200, height = 2400, res=300)
p
dev.off()

