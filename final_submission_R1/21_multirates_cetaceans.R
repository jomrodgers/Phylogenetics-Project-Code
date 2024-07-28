#cetaceans
#multirates - exploring rate heterogeneity
#Fig 5b

rm(list=ls())

setwd("C:/Users/Lars Schmitz/Dropbox/@JIM2_manuscript_PCM/final_submission_R1")

#load packages
library (phytools)
library (strap)
library (geiger)
library (paleotree)
library (picante)
source("multirate_plot.R")

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
cphy01 <- ccomp$phy
cdat4multi <- ccomp$data
cdat4multi[,2:5] <- apply(cdat4multi[,2:5], 2, as.numeric)
cdat4multi <- cdat4multi[cphy01$tip.label,]

cages <- cbind(FAD=cdat4multi$FAD, LAD=cdat4multi$LAD)
rownames(cages) <- cdat4multi$Taxon

ctree <- timePaleoPhy(cphy01, cages, type = "equal", vartime = 1, randres = T, 
                      ntrees=1, dateTreatment= "firstLast", add.term=T)

phy <- ladderize(ctree)
phy$edge.length <- phy$edge.length/max(nodeHeights(phy))
max(nodeHeights(phy))

size <- data.frame(taxon=as.character(cdat4multi[,1]), 
                   size=as.numeric(cdat4multi[,2]))
rownames(size) <- size$taxon

size <- size[phy$tip.label,]

head(size)


################

#data for phytools' multirateBM()
trait <- as.vector(log10(size$size))
names(trait) <- phy$tip.label

################

fit1.ml<-multirateBM(phy, trait, lambda=1)


plot(fit1.ml, ftype="off", digits=4, 
     mar=c(1.1,1.1,4.1,1.1))

plot.multirateBM_hcl(fit1.ml, ftype="off", digits=4, 
                     mar=c(1.1,1.1,4.1,1.1))

title(main=expression(paste("Estimated rates (",
                            sigma[BM]^2,"), for ",lambda," = 1")),adj=0.1)

###

fit2.ml<-multirateBM(phy, trait, lambda=0.1)


plot(fit2.ml, ftype="off", digits=4, 
     mar=c(1.1,1.1,4.1,1.1))

plot.multirateBM_hcl(fit2.ml, ftype="off", digits=4, 
                     mar=c(1.1,1.1,4.1,1.1))

title(main=expression(paste("Estimated rates (",
                            sigma[BM]^2,"), for ",lambda," = 0.1")),adj=0.1)

pdf("multirate_cetaceans_lambda_0_1.pdf")
plot.multirateBM_hcl(fit1.ml, ftype="off", digits=4, 
                     mar=c(1.1,1.1,4.1,1.1))

title(main=expression(paste("Estimated rates (",
                            sigma[BM]^2,"), for ",lambda," = 1")),adj=0.1)
dev.off()


pdf("multirate_cetaceans_lambda_0_1_hcl.pdf")

plot.multirateBM_hcl(fit2.ml, ftype="off", digits=4, 
                     mar=c(1.1,1.1,4.1,1.1))

title(main=expression(paste("Estimated rates (",
                            sigma[BM]^2,"), for ",lambda," = 0.1")),adj=0.1)
dev.off()


###


fit3.ml<-multirateBM(phy, trait, lambda=0.01)
pdf("multirate_cetaceans_lambda_0_01_hcl.pdf")

plot.multirateBM_hcl(fit3.ml, ftype="off", digits=4, 
                     mar=c(1.1,1.1,4.1,1.1))

title(main=expression(paste("Estimated rates (",
                            sigma[BM]^2,"), for ",lambda," = 0.01")),adj=0.1)
dev.off()


###


fit4.ml<-multirateBM(phy, trait, lambda=10)
pdf("multirate_cetaceans_lambda_0_01_hcl.pdf")

plot.multirateBM_hcl(fit4.ml, ftype="off", digits=4, 
                     mar=c(1.1,1.1,4.1,1.1))

title(main=expression(paste("Estimated rates (",
                            sigma[BM]^2,"), for ",lambda," = 10")),adj=0.1)
dev.off()

