#traitgrams of cetacean and ichthyosaur body size proxies (normalized for comparison)
#combined into a single diagram to enhance direct comparison
#Figure 4
#color and design edits in Illustrator

rm(list=ls())

setwd("C:/Users/jmr95/OneDrive/Documents/AP&P_data/science.abf5787_data_s6/final_submission_R1/")

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
iskull <- (iskull-min(iskull)) / (max(iskull)-min(iskull))
force.multiPhylo = TRUE
itreeset <- timePaleoPhy(iphy4asr, iages, type = "equal", vartime = 1, randres = T, 
                         ntrees=100, dateTreatment= "minMax")

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

##########CETACEANS


#load tree
cphy00 <- read.nexus("C:/Users/jmr95/OneDrive/Documents/AP&P_data/science.abf5787_data_s6/final_submission_R1/cetacean_tree_0529.nex")
cphy00$root.time <- 56

#load and prep data
cet.data00 <- read.csv("C:/Users/jmr95/OneDrive/Documents/AP&P_data/cetacea_data.csv")
head(cet.data00)
rownames(cet.data00) <- cet.data00$Taxon
plot(tree, show.node.label = TRUE)

#match data and tree
ccomp <- match.phylo.data(cphy00, cet.data00)
cphy01 <- ccomp$phy
cdat01 <- ccomp$data
cdat01[,2:5] <- apply(cdat01[,2:5], 2, as.numeric)
cdat01 <- cdat01[cphy01$tip.label,]
cages <- cbind(FAD=cdat01$FAD, LAD=cdat01$LAD)
rownames(cages) <- cdat01$Taxon

cskull <-setNames(cdat01[,2],rownames(cdat01))
cskull <- (cskull-min(cskull)) / (max(cskull)-min(cskull))

ctreeset <- timePaleoPhy(cphy01, cages, type = "equal", vartime = 1, randres = T, 
                         ntrees=10, dateTreatment= "minMax")
ch <- sapply(ctreeset, function(x) max(nodeHeights(x)))
ii <- order(ch, decreasing=TRUE)
ctreeset <- ctreeset[ii]
ch <- ch[ii]
ctreeset <- mapply(function(x, r) { x$root.edge <- r; x },
                   x=ctreeset, r=max(ch)-ch, SIMPLIFY=FALSE)
ctreeset <- lapply(ctreeset, rootedge.to.singleton)

ctree <- timePaleoPhy(cphy01, cages, type = "equal", vartime = 1, randres = T, 
                      ntrees=1, dateTreatment= "firstLast", add.term=T)
ctree$root.edge <- max(ch)-max(nodeHeights(ctree))
ctree <- rootedge.to.singleton(ctree)


########################################################################################
pdf("Fig_4_combined_traitgrams_with100trees.pdf", 6, 8)

par(xaxt="n",yaxt="n",mar=c(5.1,5.1,2.1,1.1))


for(i in 1:length(ctreeset)){ 
#cetaceans
    ca <- fastAnc(collapse.singles(ctreeset[[i]]), cskull)
    ca <- setNames(c(ca[1], ca), 1:ctreeset[[i]]$Nnode + Ntip(ctreeset[[i]]))
    ctreeset[[i]] <- paintSubTree(ctreeset[[i]],node=Ntip(ctreeset[[i]])+2,
                                  state="1",anc.state="0")
    ccols <- setNames(c("blue", make.transparent(palette()[4], 1/length(ctreeset))),0:1)
    phenogram(ctreeset[[i]], c(cskull, ca), 
              spread.labels=TRUE, 
              fsize=1, 
              color=ccols, add=(i>1),
              lwd=4,
              quiet=TRUE, 
              xlim=c(0,175),
              xlab="", ylab="")
  }
  
  b <- fastAnc(collapse.singles(ctree), cskull)
  b <- setNames(c(b[1],b), 1:ctree$Nnode + Ntip(ctree))
  phenogram(ctree, c(cskull, b), 
            spread.labels=TRUE, 
            fsize=1,
            color="#2297E6", add=TRUE,
            lwd=3, 
            quiet=TRUE, 
            xlim=c(0,175),
            xlab="", ylab="")
  phenogram(ctree, c(cskull, b), 
            spread.labels=TRUE, 
            fsize=1,
            color="white", add=TRUE,
            lwd=1,
            quiet=TRUE, 
            xlim=c(0,175),
            xlab="", ylab="")
  
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
            fsize=0,
            color=icols, add=TRUE,
            lwd=4,
            quiet=TRUE, 
            xlim=c(0,175),
            xlab="", ylab="")
}

itree$edge.length <- itree$edge.length +1e-7
itree.n <- itree
itree.n$edge.length <- itree.n$edge.length/max(nodeHeights(itree.n))
a <- anc.ML(collapse.singles(itree.n), iskull, model="EB")
a <- setNames(c(a$ace[1],a$ace), 1:itree.n$Nnode + Ntip(itree.n))
phenogram(itree, c(iskull, a), 
          spread.labels=TRUE, 
          fsize=0,
          color="#DF536B", add=TRUE,
          lwd=3, 
          quiet=TRUE, 
          xlim=c(0,175),
          xlab="", ylab="")
phenogram(itree, c(iskull, a), 
          spread.labels=TRUE, 
          fsize=0,
          color="white", add=TRUE,
          lwd=2,
          quiet=TRUE, 
          xlim=c(0,175),
          xlab="", ylab="")




par(xaxt="s",yaxt="s",font.lab=2)
axis(1)
title(xlab="time since the root in million years", cex.lab=1)
axis(2)
title(ylab="relative body size proxy", cex.lab=1)

dev.off()
########################################################################################

