#ichthyosaurs
#bayou
#Fig 6a
#bar plots generated separately and combined woth bayou regimes in Illustrator

rm(list=ls())

setwd("C:/Users/Lars Schmitz/Dropbox/@JIM2_manuscript_PCM/final_submission_R1")

#load packages
library (phytools)
library (strap)
library (geiger)
library (paleotree)
library (picante)
library (bayou)

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
iphy4bayou <- ichthyocomp$phy
idat4bayou <- ichthyocomp$data
idat4bayou[,2:5] <- apply(idat4bayou[,2:5], 2, as.numeric)
idat4bayou <- idat4bayou[iphy4bayou$tip.label,]

iages <- cbind(FAD=idat4bayou$FAD, LAD=idat4bayou$LAD)
rownames(iages) <- idat4bayou$Taxon

itree <- timePaleoPhy(iphy4bayou, iages, type = "equal", vartime = 1, randres = T, 
                      ntrees=1, dateTreatment= "firstLast", add.term=T)

phy <- ladderize(itree)

phy <- reorder(phy, order = "postorder")

size <- data.frame(taxon=as.character(idat4bayou[,1]), 
                   size=as.numeric(idat4bayou[,2]))
rownames(size) <- size$taxon

size <- size[phy$tip.label,]

head(size)


################

#data for bayOU
trait <- as.vector(log10(size$size))
names(trait) <- phy$tip.label

range(trait)

max(nodeHeights(phy))

#kmax=length(phy$tip.label)*2-2
#phy$edge.length <- phy$edge.length/10


#prior
prior <- make.prior(phy, 
                      dists=list(dalpha="dhalfcauchy", 
                                 dsig2="dhalfcauchy", 
                                 dk="cdpois", 
                                 dtheta="dnorm"),
                      param=list(dalpha=list(scale=0.1), 
                                 dsig2=list(scale=0.1),
                                 dk=list(lambda=10, kmax=50), 
                                 dsb=list(bmax=1, prob=1), 
                                 dtheta=list(mean=mean(trait), 
                                 sd=1.5*sd(trait)))
)

startpars <- priorSim(prior, phy, plot=TRUE)$pars[[1]]
prior(startpars)



#chain 1
fit1 <- bayou.makeMCMC(phy, trait, model="OU", prior=prior, file.dir="./chain1", plot.freq=NULL, ticker.freq=10000)
fit1$run(3000000)

chain <- fit1$load()
chain <- set.burnin(chain, 0.3)
#plot(chain, auto.layout=F)
out <- summary(chain)


plotSimmap.mcmc(chain, burnin = 0.3, edge.width=4, pp.cutoff = 0.3)
plotBranchHeatMap(phy, chain, "theta", edge.width=4, burnin = 0.3, pal = cm.colors)
phenogram.density(phy, trait, burnin = 0.3, chain, pp.cutoff = 0.3)

write.csv(out$statistics, "./ibayou_results/bayOU_a_stats.csv")
write.csv(out$branch.posteriors, "./ibayou_results/bayOU_a_posteriors.csv")



#plots

pdf("./ibayou_results/bayOU_regimes_a.pdf")
par(mfrow=c(1,1))
plotSimmap.mcmc(chain, burnin=0.3, pp.cutoff = 0.5, lwd=4, edge.type="regimes", pal=hcl.colors, show.tip.label=F, circle.col="black", cex=0.5)
dev.off()

pdf("./ibayou_results/bayOU_regimes_a_cutoff03.pdf")
par(mfrow=c(1,1))
plotSimmap.mcmc(chain, burnin=0.3, pp.cutoff = 0.3, lwd=4, edge.type="regimes", pal=hcl.colors, show.tip.label=F, circle.col="black", cex=0.5)
dev.off()

pdf("./ibayou_results/bayOU_regimes_a_cutoff04.pdf")
par(mfrow=c(1,1))
plotSimmap.mcmc(chain, burnin=0.3, pp.cutoff = 0.4, lwd=4, edge.type="regimes", pal=hcl.colors, show.tip.label=F, circle.col="black", cex=0.5)
dev.off()

pdf("./ibayou_results/bayOU_regimes_a_cutoff04_labels.pdf")
par(mfrow=c(1,1))
plotSimmap.mcmc(chain, burnin=0.3, pp.cutoff = 0.4, lwd=4, edge.type="regimes", pal=hcl.colors, show.tip.label=T, circle.col="black", cex=0.5)
dev.off()


pdf("./ibayou_results/bayOU_heatmap_a.pdf")
par(mfrow=c(1,1))
plotBranchHeatMap(phy, chain, "theta", edge.width=4, burnin = 0.3, pal = hcl.colors)
dev.off()




startpars <- priorSim(prior, phy, plot=TRUE)$pars[[1]]
prior(startpars)


#chain 2
fit2 <- bayou.makeMCMC(phy, trait, model="OU", prior=prior, file.dir="./chain2", plot.freq=NULL, ticker.freq=10000)
fit2$run(3000000)

chain2 <- fit2$load()
chain2 <- set.burnin(chain2, 0.3)
out2 <- summary(chain2)

plotSimmap.mcmc(chain2, burnin = 0.3, pp.cutoff = 0.3)
plotBranchHeatMap(phy, chain2, "theta", burnin = 0.3, edge.width=4, pal = cm.colors)
phenogram.density(phy, trait, burnin = 0.3, chain2, pp.cutoff = 0.3)


write.csv(out2$statistics, "./ibayou_results/bayOU_b_stats.csv")
write.csv(out2$branch.posteriors, "./ibayou_results/bayOU_b_posteriors.csv")


#plots

pdf("./ibayou_results/bayOU_regimes_b.pdf")
par(mfrow=c(1,1))
plotSimmap.mcmc(chain2, burnin=0.3, pp.cutoff = 0.5, lwd=4, edge.type="regimes", pal=hcl.colors, show.tip.label=F, circle.col="black", cex=0.5)
dev.off()

pdf("./ibayou_results/bayOU_regimes_b_cutoff03.pdf")
par(mfrow=c(1,1))
plotSimmap.mcmc(chain2, burnin=0.3, pp.cutoff = 0.3, lwd=4, edge.type="regimes", pal=hcl.colors, show.tip.label=F, circle.col="black", cex=0.5)
dev.off()

pdf("./ibayou_results/bayOU_regimes_b_cutoff04.pdf")
par(mfrow=c(1,1))
plotSimmap.mcmc(chain2, burnin=0.3, pp.cutoff = 0.4, lwd=4, edge.type="regimes", pal=hcl.colors, show.tip.label=F, circle.col="black", cex=0.5)
dev.off()

pdf("./ibayou_results/bayOU_regimes_b_cutoff04_labels.pdf")
par(mfrow=c(1,1))
plotSimmap.mcmc(chain2, burnin=0.3, pp.cutoff = 0.4, lwd=4, edge.type="regimes", pal=hcl.colors, show.tip.label=T, circle.col="black", cex=0.5)
dev.off()


pdf("./ibayou_results/bayOU_heatmap_b.pdf")
par(mfrow=c(1,1))
plotBranchHeatMap(phy, chain2, "theta", edge.width=4, burnin = 0.3, pal = hcl.colors)
dev.off()



#diagnostics

pdf("./ibayou_results/bayOU_gelman.pdf")
par(mfrow=c(3,1))
RlnL <- gelman.R("lnL", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
abline(v=1000000, lty=2, col="red")
Ralpha <- gelman.R("alpha", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
abline(v=1000000, lty=2, col="red")
Rsig2 <- gelman.R("sig2", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
abline(v=1000000, lty=2, col="red")
dev.off()

L1 <- Lposterior(chain, phy, burnin=0.3)
L2 <- Lposterior(chain2, phy, burnin=0.3)


pdf("./ibayou_results/bayOU_posterior_diagnostics.pdf", useDingbats=F)
par(mfrow=c(1,1))
plot(L1$pp,L2$pp, xlim=c(0,1), ylim=c(0,1), xlab="Chain 1", ylab="Chain 2")
curve(1*x, add=TRUE, lty=2)
dev.off()


pdf("./ibayou_results/bayOU_pp_prior_a.pdf")
par(mfrow=c(1,1))
plot(density(L1[,1]), main="PP compared to prior")
abline(v=dsb(1, ntips = length(phy$tip.label), bmax = 1, prob = 1, log = FALSE), col="red", lty=2)
dev.off()

pdf("./ibayou_results/bayOU_pp_prior_b.pdf")
par(mfrow=c(1,1))
plot(density(L2[,1]), main="PP compared to prior")
abline(v=dsb(1, ntips = length(phy$tip.label), bmax = 1, prob = 1, log = FALSE), col="red", lty=2)
dev.off()



