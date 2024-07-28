#ichthyosaurs
#traitgram
#labelled
#accounting for tree distribution
#log10 skull length
#multiple colour options; final color edited in Illustrator
#Figure S10

rm(list=ls())

setwd("C:/Users/Lars Schmitz/Dropbox/@JIM2_manuscript_PCM/final_submission_R1")

#load packages
library (ape)
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
ichthyo.data00 <- read.csv("C:/Users/jmr95/OneDrive/Documents/AP&P_data/ichthyosaur_data.csv")
rownames(ichthyo.data00) <- ichthyo.data00$Taxon

#match data and tree
ichthyocomp <- match.phylo.data(phy00, ichthyo.data00)
iphy4asr <- ichthyocomp$phy
idat4asr <- ichthyocomp$data
idat4asr[,2:5] <- apply(idat4asr[,2:5], 2, as.numeric)
idat4asr <- idat4asr[iphy4asr$tip.label,]

iages <- cbind(FAD=idat4asr$FAD, LAD=idat4asr$LAD)
rownames(iages) <- idat4asr$Taxon

iskull <- setNames(idat4asr[,2],rownames(idat4asr))
iskull <- log10(iskull)
#iskull <- (iskull-min(iskull)) / (max(iskull)-min(iskull)) #use when plotting normlaized data

itreeset <- timePaleoPhy(iphy4asr, iages, type = "equal", vartime = 1, randres = T, 
                         ntrees=20, dateTreatment= "minMax")

#h <- sapply(itreeset,function(x) max(nodeHeights(x)))
#ii <- order(h,decreasing=TRUE)
#itreeset <- itreeset[ii]
#h <- h[ii]
itreeset <- mapply(function(x,r) { x$root.edge <- r; x },
                   x=itreeset, r=0, SIMPLIFY=FALSE) #max(h)-h
itreeset <- lapply(itreeset, rootedge.to.singleton)


itree <- timePaleoPhy(iphy4asr, iages, type = "equal", vartime = 1, randres = T, 
                      ntrees=1, dateTreatment= "firstLast", add.term=T)
itree$root.edge <- 0
itree <- rootedge.to.singleton(itree)



########################################################################################
pdf("Fig_S10_ichthyosaur_traitgram_with_20trees_labelled_presubmission_red.pdf", 6, 8)

par(xaxt="n",yaxt="n",mar=c(5.1,5.1,2.1,1.1))

for(i in 1:length(itreeset)){  
#ichthyosaurs
  itreeset[[i]]$edge.length <- itreeset[[i]]$edge.length +1e-7
  itreeset.n <- itreeset
  itreeset.n[[i]]$edge.length <- itreeset.n[[1]]$edge.length/max(nodeHeights(itreeset.n[[i]]))
  ia <- anc.ML(collapse.singles(itreeset.n[[i]]), iskull, model="EB")
  ia <- setNames(c(ia$ace[1],ia$ace), 1:itreeset.n[[i]]$Nnode + Ntip(itreeset.n[[i]]))
  itreeset[[i]] <- paintSubTree(itreeset[[i]],node=Ntip(itreeset[[i]])+2,
                                state="1",anc.state="0")
  icols <- setNames(c("transparent", make.transparent(palette()[2], 1/length(itreeset))),0:1)
  phenogram(itreeset[[i]], c(iskull, ia), 
            spread.labels=TRUE, 
            fsize=0.3, ftype="i",
            color=icols, add=(i>1),
            lwd=4,
            quiet=TRUE,
            xlim=c(0,250),
            xlab="", ylab="")
}

par(xaxt="s",yaxt="s",font.lab=2)
axis(1)
title(xlab="time since the root in million years", cex.lab=1)
axis(2)
title(ylab="relative body size proxy", cex.lab=1)

dev.off()
########################################################################################



########################################################################################
pdf("Fig_S10_ichthyosaur_traitgram_with_20trees_labelled_presubmission_blue.pdf", 6, 8)

par(xaxt="n",yaxt="n",mar=c(5.1,5.1,2.1,1.1))

for(i in 1:length(itreeset)){  
  #ichthyosaurs
  itreeset[[i]]$edge.length <- itreeset[[i]]$edge.length +1e-7
  itreeset.n <- itreeset
  itreeset.n[[i]]$edge.length <- itreeset.n[[1]]$edge.length/max(nodeHeights(itreeset.n[[i]]))
  ia <- anc.ML(collapse.singles(itreeset.n[[i]]), iskull, model="EB")
  ia <- setNames(c(ia$ace[1],ia$ace), 1:itreeset.n[[i]]$Nnode + Ntip(itreeset.n[[i]]))
  itreeset[[i]] <- paintSubTree(itreeset[[i]],node=Ntip(itreeset[[i]])+2,
                                state="1",anc.state="0")
  icols <- setNames(c("transparent", make.transparent(palette()[4], 1/length(itreeset))),0:1)
  phenogram(itreeset[[i]], c(iskull, ia), 
            spread.labels=TRUE, 
            fsize=0.3, ftype="i",
            color=icols, add=(i>1),
            lwd=4,
            quiet=TRUE,
            xlim=c(0,250),
            xlab="", ylab="")
}

par(xaxt="s",yaxt="s",font.lab=2)
axis(1)
title(xlab="time since the root in million years", cex.lab=1)
axis(2)
title(ylab="relative body size proxy", cex.lab=1)

dev.off()

################
