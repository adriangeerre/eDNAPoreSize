# Library
library(tidyverse)

# Inside the folders (03-Blast & 04-FilterBlast) we run:
# ls | grep -v ".tsv" | xargs -I {} sh -c "cat {}/* | cut -f 4 > lens.{}.tsv" & # Alig. Length
# ls | grep -v ".tsv" | xargs -I {} sh -c "cat {}/* | cut -f 3 > pident.{}.tsv" & # Perc. Identity
# ls | grep -v ".tsv" | xargs -I {} sh -c "cat {}/* | cut -f 17 > qcovs.{}.tsv" & # Query Coverage

## Lengths
#  -------
histo_lens <- function(sample, color) {
    # Rename sample
    oldnames <- c('GP_1','GP_2','GP_3','A0_2','B0_2','C0_2','A1_2','B1_2','C1_2','A5_0','B5_0','C5_0','A8_0','B8_0','C8_0')
    newnames <- c('EN0.2A','EN0.2B','EN0.2C','OP0.2A','OP0.2B','OP0.2C','OP1.2A','OP1.2B','OP1.2C','OP5.0A','OP5.0B','OP5.0C','OP8.0A','OP8.0B','OP8.0C')
    names(newnames) <- oldnames
    name <- newnames[sample]
    # Read files
    pre = read.table(paste("03-Blast/lens.", sample, ".tsv", sep=""))
    pst = read.table(paste("04-FilterBlast/lens.post_filter.", sample, ".tsv", sep=""))
    # Histogram
    p <- ggplot() + 
        geom_histogram(data=pre, aes(x=V1), bins=100, fill=color, color="lightgray") +
        geom_histogram(data=pst, aes(x=V1), bins=100, fill=color, color = "black") +
        labs(x="Alignment length", y="Count", title=paste("Length: ", name, sep="")) +
        theme_classic() %+replace% theme(axis.text = element_text(size=16),
                                 axis.title = element_text(size=18),
                                 title = element_text(size=22, face="bold"))
    # Return
    return(p)
}

# GP
p <- histo_lens("GP_1", "#F8766D")
# Save figure
png(file=paste("07-Plots/lens.GP_1.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_lens("GP_2", "#F8766D")
# Save figure
png(file=paste("07-Plots/lens.GP_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_lens("GP_3", "#F8766D")
# Save figure
png(file=paste("07-Plots/lens.GP_3.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_lens("A0_2", "#B79F00")
# Save figure
png(file=paste("07-Plots/lens.A0_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_lens("B0_2", "#B79F00")
# Save figure
png(file=paste("07-Plots/lens.B0_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_lens("C0_2", "#B79F00")
# Save figure
png(file=paste("07-Plots/lens.C0_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_lens("A1_2","#00BF7D")
# Save figure
png(file=paste("07-Plots/lens.A1_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_lens("B1_2","#00BF7D")
# Save figure
png(file=paste("07-Plots/lens.B1_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_lens("C1_2","#00BF7D")
# Save figure
png(file=paste("07-Plots/lens.C1_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_lens("A5_0","#00B0F6")
# Save figure
png(file=paste("07-Plots/lens.A5_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_lens("B5_0","#00B0F6")
# Save figure
png(file=paste("07-Plots/lens.B5_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_lens("C5_0","#00B0F6")
# Save figure
png(file=paste("07-Plots/lens.C5_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_lens("A8_0","#E76BF3")
# Save figure
png(file=paste("07-Plots/lens.A8_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_lens("B8_0","#E76BF3")
# Save figure
png(file=paste("07-Plots/lens.B8_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_lens("C8_0","#E76BF3")
# Save figure
png(file=paste("07-Plots/lens.C8_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

## Identity
#  --------

histo_pident <- function(sample, color) {
    # Rename sample
    oldnames <- c('GP_1','GP_2','GP_3','A0_2','B0_2','C0_2','A1_2','B1_2','C1_2','A5_0','B5_0','C5_0','A8_0','B8_0','C8_0')
    newnames <- c('EN0.2A','EN0.2B','EN0.2C','OP0.2A','OP0.2B','OP0.2C','OP1.2A','OP1.2B','OP1.2C','OP5.0A','OP5.0B','OP5.0C','OP8.0A','OP8.0B','OP8.0C')
    names(newnames) <- oldnames
    name <- newnames[sample]
    # Read files
    pre = read.table(paste("03-Blast/pident.", sample, ".tsv", sep=""))
    pst = read.table(paste("04-FilterBlast/pident.post_filter.", sample, ".tsv", sep=""))
    # Histogram
    p <- ggplot() + 
        geom_histogram(data=pre, aes(x=V1), bins=100, fill=color, color="lightgray") +
        geom_histogram(data=pst, aes(x=V1), bins=100, fill=color, color = "black") +
        labs(x="Perc. Identity", y="Count", title=paste("Identity: ", name, sep="")) +
        theme_classic() %+replace% theme(axis.text = element_text(size=16),
                                 axis.title = element_text(size=18),
                                 title = element_text(size=22, face="bold"))
    # Return
    return(p)
}

p <- histo_pident("GP_1", "#F8766D")
# Save figure
png(file=paste("07-Plots/pident.GP_1.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_pident("GP_2", "#F8766D")
# Save figure
png(file=paste("07-Plots/pident.GP_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_pident("GP_3", "#F8766D")
# Save figure
png(file=paste("07-Plots/pident.GP_3.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_pident("A0_2", "#B79F00")
# Save figure
png(file=paste("07-Plots/pident.A0_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_pident("B0_2", "#B79F00")
# Save figure
png(file=paste("07-Plots/pident.B0_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_pident("C0_2", "#B79F00")
# Save figure
png(file=paste("07-Plots/pident.C0_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_pident("A1_2","#00BF7D")
# Save figure
png(file=paste("07-Plots/pident.A1_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_pident("B1_2","#00BF7D")
# Save figure
png(file=paste("07-Plots/pident.B1_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_pident("C1_2","#00BF7D")
# Save figure
png(file=paste("07-Plots/pident.C1_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_pident("A5_0","#00B0F6")
# Save figure
png(file=paste("07-Plots/pident.A5_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_pident("B5_0","#00B0F6")
# Save figure
png(file=paste("07-Plots/pident.B5_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_pident("C5_0","#00B0F6")
# Save figure
png(file=paste("07-Plots/pident.C5_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_pident("A8_0","#E76BF3")
# Save figure
png(file=paste("07-Plots/pident.A8_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_pident("B8_0","#E76BF3")
# Save figure
png(file=paste("07-Plots/pident.B8_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_pident("C8_0","#E76BF3")
# Save figure
png(file=paste("07-Plots/pident.C8_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

## Coverage
#  --------

histo_qcovs <- function(sample, color) {
    # Rename sample
    oldnames <- c('GP_1','GP_2','GP_3','A0_2','B0_2','C0_2','A1_2','B1_2','C1_2','A5_0','B5_0','C5_0','A8_0','B8_0','C8_0')
    newnames <- c('EN0.2A','EN0.2B','EN0.2C','OP0.2A','OP0.2B','OP0.2C','OP1.2A','OP1.2B','OP1.2C','OP5.0A','OP5.0B','OP5.0C','OP8.0A','OP8.0B','OP8.0C')
    names(newnames) <- oldnames
    name <- newnames[sample]
    # Read files
    pre = read.table(paste("03-Blast/qcovs.", sample, ".tsv", sep=""))
    pst = read.table(paste("04-FilterBlast/qcovs.post_filter.", sample, ".tsv", sep=""))
    # Histogram
    p <- ggplot() + 
        geom_histogram(data=pre, aes(x=V1), bins=100, fill=color, color="lightgray") +
        geom_histogram(data=pst, aes(x=V1), bins=100, fill=color, color = "black") +
        labs(x="Query coverage", y="Count", title=paste("Coverage: ", name, sep="")) +
        xlim(0,100) +
        theme_classic() %+replace% theme(axis.text = element_text(size=16),
                                 axis.title = element_text(size=18),
                                 title = element_text(size=22, face="bold"))
    # Return
    return(p)
}

p <- histo_qcovs("GP_1", "#F8766D")
# Save figure
png(file=paste("07-Plots/qcovs.GP_1.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_qcovs("GP_2", "#F8766D")
# Save figure
png(file=paste("07-Plots/qcovs.GP_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_qcovs("GP_3", "#F8766D")
# Save figure
png(file=paste("07-Plots/qcovs.GP_3.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_qcovs("A0_2", "#B79F00")
# Save figure
png(file=paste("07-Plots/qcovs.A0_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_qcovs("B0_2", "#B79F00")
# Save figure
png(file=paste("07-Plots/qcovs.B0_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_qcovs("C0_2", "#B79F00")
# Save figure
png(file=paste("07-Plots/qcovs.C0_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_qcovs("A1_2","#00BF7D")
# Save figure
png(file=paste("07-Plots/qcovs.A1_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_qcovs("B1_2","#00BF7D")
# Save figure
png(file=paste("07-Plots/qcovs.B1_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_qcovs("C1_2","#00BF7D")
# Save figure
png(file=paste("07-Plots/qcovs.C1_2.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_qcovs("A5_0","#00B0F6")
# Save figure
png(file=paste("07-Plots/qcovs.A5_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_qcovs("B5_0","#00B0F6")
# Save figure
png(file=paste("07-Plots/qcovs.B5_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_qcovs("C5_0","#00B0F6")
# Save figure
png(file=paste("07-Plots/qcovs.C5_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_qcovs("A8_0","#E76BF3")
# Save figure
png(file=paste("07-Plots/qcovs.A8_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_qcovs("B8_0","#E76BF3")
# Save figure
png(file=paste("07-Plots/qcovs.B8_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()

p <- histo_qcovs("C8_0","#E76BF3")
# Save figure
png(file=paste("07-Plots/qcovs.C8_0.png", sep=""), width=2300, height=2000, res=300)
p
dev.off()