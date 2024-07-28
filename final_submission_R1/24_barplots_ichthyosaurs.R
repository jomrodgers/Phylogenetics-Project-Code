#ichthyosaurs
#barplots
#Fig 6a right panel

setwd("C:/Users/jmr95/OneDrive/Documents/AP&P_data")

#load packages
library (phytools)
library (geiger)
library (picante)


#############ICHTHYOSAURS
#load tree
phy00 <- read.tree("ichthyosaur_tree.tre")
phy00$root.time <- 251.902
plot(phy00, show.tip.label=F)

#load and prep data
ichthyo.data00 <- read.csv("ichthyosaur_data.csv")
rownames(ichthyo.data00) <- ichthyo.data00$Taxon

#match data and tree
ichthyocomp <- match.phylo.data(phy00, ichthyo.data00)
itree <- ichthyocomp$phy
idata <- ichthyocomp$data
idata[,2:5] <- apply(idata[,2:5], 2, as.numeric)
idata <- idata[itree$tip.label,]


phy <- ladderize(itree)
phy$edge <- phy$edge[nrow(mat):1,]
#phy <- itree

size <- data.frame(taxon=as.character(idata[,1]), 
                   size=as.numeric(idata[,2]))
rownames(size) <- size$taxon

size <- size[phy$tip.label,]
head(size)


################

#data for barplot
trait <- as.vector(size$size)
names(trait) <- phy$tip.label


pdf("barplots_ichthyosaurs.pdf", height=5, width=2)
barplot(rev(trait), horiz=T, axisnames=F, border=NA)
dev.off()


################

plotTree.barplot(itree, trait)
