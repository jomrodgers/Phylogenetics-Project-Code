#cetaceans
#barplots
#Fig 6b right panel

setwd("C:/Users/Lars Schmitz/Dropbox/@JIM2_manuscript_PCM/final_submission_R1")

#load packages
library (phytools)
library (strap)
library (geiger)
library (paleotree)
library (picante)
library (bayou)

#############CETACEANS
#load tree
cphy00 <- read.nexus("cetacean_tree_0529.nex")
cphy00$root.time <- 56

#load and prep data
cet.data00 <- read.csv("cetacea_data.csv")
head(cet.data00)
rownames(cet.data00) <- cet.data00$Taxon

#match data and tree
ccomp <- match.phylo.data(cphy00, cet.data00)
ctree <- ccomp$phy
cdata <- ccomp$data
cdata[,2:5] <- apply(cdata[,2:5], 2, as.numeric)
cdata <- cdata[ctree$tip.label,]


phy <- ladderize(ctree)
#phy$edge <- phy$edge[nrow(mat):1,]
#phy <- itree

size <- data.frame(taxon=as.character(cdata[,1]), 
                   size=as.numeric(cdata[,2]))
rownames(size) <- size$taxon

size <- size[phy$tip.label,]
head(size)


################

#data for barplot
trait <- as.vector(size$size*10) #convert to mm!
names(trait) <- phy$tip.label


pdf("barplots_cetaceans.pdf", height=5, width=2)
barplot(rev(trait), horiz=T, axisnames=F, border=NA)
dev.off()

###

plotTree.barplot(ctree, trait)
