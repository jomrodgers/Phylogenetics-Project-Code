#Timetrees for ichthyosaurs and cetaceans
#Needed for Figure 3 (composed in Illustrator)

rm(list=ls())

setwd("C:/Users/jmr95/OneDrive/Documents/AP&P_data")

#load packages
library (phytools)
library (strap)
library (geiger)
library (paleotree)
library (picante)

#############ICHTHYOSAURS
#load tree
phy00 <- read.tree("C:/Users/jmr95/OneDrive/Documents/AP&P_data/ichthyosaur_tree.tre")
phy00$root.time <- 251.902
plot(phy00, show.tip.label=F)

#load and prep data
ichthyo.data00 <- read.csv("ichthyosaur_data.csv")
rownames(ichthyo.data00) <- ichthyo.data00$Taxon

#match data and tree
ichthyocomp <- match.phylo.data(phy00, ichthyo.data00)
itre0 <- ichthyocomp$phy
idat0 <- ichthyocomp$data
idat0[,2:5] <- apply(idat0[,2:5], 2, as.numeric)
idat0 <- idat0[itre0$tip.label,]
iages <- cbind(FAD=idat0$FAD, LAD=idat0$LAD)
rownames(iages) <- idat0$Taxon

### time tree ichthyosaurs
itree <- DatePhylo(itre0, iages, method="equal", rlen=1, add.terminal = T)

##########CETACEANS

#load tree
cphy00 <- read.nexus("science.abf5787_data_s6/final_submission_R1/cetacean_tree_0529.nex")
cphy00$root.time <- 56

#load and prep data
cet.data00 <- read.csv("cetacea_data.csv")
head(cet.data00)
rownames(cet.data00) <- cet.data00$Taxon

#match data and tree
ccomp <- match.phylo.data(cphy00, cet.data00)
cphy01 <- ccomp$phy
cdat01 <- ccomp$data
cdat01[,2:5] <- apply(cdat01[,2:5], 2, as.numeric)
cdat01 <- cdat01[cphy01$tip.label,]
cages <- cbind(FAD=cdat01$FAD, LAD=cdat01$LAD)
rownames(cages) <- cdat01$Taxon

### time tree cetaceans
ctree <- DatePhylo(cphy01, cages, method="equal", rlen=1, add.terminal = T)

#############PLOTTING TO PDF

pdf("Fig_2_ichthyosaur_tree_timescaled_09_18_2020.pdf", 8.5, 11)
geoscalePhylo(tree=itree, ages=iages,
              units=c("Period","Epoch"), boxes="Epoch",
              x.lim=c(0, 260),
              direction="rightwards",
              cex.tip=0.6,
              cex.ts=0.6, tick.scale="no", 
              label.offset=0, show.tip.label=T,
              quat.rm=TRUE)
dev.off()

pdf("Fig_2_cetacean_tree_timescaled_09_18_2020.pdf", 8.5, 11)
geoscalePhylo(tree=ctree, ages=cages, quat.rm = TRUE,
              units=c("Period","Epoch"), boxes="Epoch",
              x.lim=c(0, 260),
              direction="rightwards",
              cex.tip=0.6,
              cex.ts=0.6, tick.scale="no", 
              label.offset=0, show.tip.label=FALSE)
dev.off()
