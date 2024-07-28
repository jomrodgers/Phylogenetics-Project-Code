#ichthyosaurs
#multirates - exploring rate heterogeneity
#Fig 5a

rm(list=ls())

setwd("C:/Users/Lars Schmitz/Dropbox/@JIM2_manuscript_PCM/final_submission_R1")

#load packages
library (phytools)
library (strap)
library (geiger)
library (paleotree)
library (picante)
source("multirate_plot.R")

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
iphy4multi <- ichthyocomp$phy
idat4multi <- ichthyocomp$data
idat4multi[,2:5] <- apply(idat4multi[,2:5], 2, as.numeric)
idat4multi <- idat4multi[iphy4multi$tip.label,]

iages <- cbind(FAD=idat4multi$FAD, LAD=idat4multi$LAD)
rownames(iages) <- idat4multi$Taxon

itree <- timePaleoPhy(iphy4multi, iages, type = "equal", vartime = 1, randres = T, 
                      ntrees=1, dateTreatment= "firstLast", add.term=T)

phy <- ladderize(itree)
phy$edge.length <- phy$edge.length/max(nodeHeights(phy))
max(nodeHeights(phy))

size <- data.frame(taxon=as.character(idat4multi[,1]), 
                   size=as.numeric(idat4multi[,2]))
rownames(size) <- size$taxon

size <- size[phy$tip.label,]

head(size)


################

#data for phytools' multirateBM()
trait <- as.vector(log10(size$size/10)) #/10 to convert to cm (as whales)
names(trait) <- phy$tip.label

################

fit1.ml<-multirateBM(phy, trait, lambda=0.1)


plot(fit1.ml, ftype="off", digits=4, 
     mar=c(1.1,1.1,4.1,1.1))

plot.multirateBM_hcl(fit.ml, ftype="off", digits=4, 
                     mar=c(1.1,1.1,4.1,1.1))

title(main=expression(paste("Estimated rates (",
                            sigma[BM]^2,"), for ",lambda," = 0.1")),adj=0.1)


pdf("multirate_ichthyosaurs_lambda_0_1.pdf")
plot.multirateBM_hcl(fit1.ml, ftype="off", digits=4, 
                     mar=c(1.1,1.1,4.1,1.1))

title(main=expression(paste("Estimated rates (",
                            sigma[BM]^2,"), for ",lambda," = 0.1")),adj=0.1)
dev.off()


###

fit2.ml<-multirateBM(phy, trait, lambda=0.01)


plot(fit2.ml, ftype="off", digits=4, 
     mar=c(1.1,1.1,4.1,1.1))

plot.multirateBM_hcl(fit2.ml, ftype="off", digits=4, 
     mar=c(1.1,1.1,4.1,1.1))

title(main=expression(paste("Estimated rates (",
                            sigma[BM]^2,"), for ",lambda," = 0.01")),adj=0.1)

pdf("multirate_ichthyosaurs_lambda_0_0_1.pdf")
plot.multirateBM_hcl(fit1.ml, ftype="off", digits=4, 
                     mar=c(1.1,1.1,4.1,1.1))

title(main=expression(paste("Estimated rates (",
                            sigma[BM]^2,"), for ",lambda," = 0.01")),adj=0.1)
dev.off()


###


fit3.ml<-multirateBM(phy, trait, lambda=0.001)


plot(fit3.ml, ftype="off", digits=4, 
     mar=c(1.1,1.1,4.1,1.1))

plot.multirateBM_hcl(fit3.ml, ftype="off", digits=4, 
                     mar=c(1.1,1.1,4.1,1.1))

title(main=expression(paste("Estimated rates (",
                            sigma[BM]^2,"), for ",lambda," = 0.001")),adj=0.1)

###


fit4.ml<-multirateBM(phy, trait, lambda=1)


plot(fit4.ml, ftype="off", digits=4, 
     mar=c(1.1,1.1,4.1,1.1))

plot.multirateBM_hcl(fit4.ml, ftype="off", digits=4, 
                     mar=c(1.1,1.1,4.1,1.1))

title(main=expression(paste("Estimated rates (",
                            sigma[BM]^2,"), for ",lambda," = 1")),adj=0.1)


###


fit5.ml<-multirateBM(phy, trait, lambda=10)


plot(fit5.ml, ftype="off", digits=4, 
     mar=c(1.1,1.1,4.1,1.1))

plot.multirateBM_hcl(fit5.ml, ftype="off", digits=4, 
                     mar=c(1.1,1.1,4.1,1.1))

title(main=expression(paste("Estimated rates (",
                            sigma[BM]^2,"), for ",lambda," = 10")),adj=0.1)
