#Cetaceans
#Model fitting
#Fig S11 H

#Re-sampling approach no. 6 to more thoroughly investigate underlying patterns of trait evolution
#Cetaceans include a large radiation of mysticeti, of which there is no evidence in ichthyosaurs, even though filter feeding may have been present in upper Triassic taxa
#In addition, thje cetacean dataset is larger than the ichthyosaur dataset.

#Do results change when filter feeders are excluded and the data are downsampled?

rm(list=ls())

setwd("C:/Users/Lars Schmitz/Dropbox/@JIM2_manuscript_PCM/final_submission_R1")

#load packages
library (phytools)
library (geiger)
library (paleotree)
library (reshape2)
library (ggplot2)
library (picante)

#load tree
cphy00 <- read.nexus("cetacean_tree_0529.nex")
cphy00$root.time <- 56

#dropping Chaeomysticeti, the filter feeders
Chaeomysticeti <- findMRCA(cphy00, tips=c("Maiabalaena_nesbittae",	"Balaenoptera_physalus"),
                           type="node")
filterfeeders <- tips(cphy00, node=Chaeomysticeti) 
phy00_nofilter <- drop.tip(cphy00, filterfeeders)

#load and prep data
cet.data00 <- read.csv("cetacea_data.csv")
head(cet.data00)
rownames(cet.data00) <- cet.data00$Taxon

#match data and tree
comp <- match.phylo.data(phy00_nofilter, cet.data00)
phy02 <- comp$phy
dat02 <- comp$data
dat02[,2:5] <- apply(dat02[,2:5], 2, as.numeric)
dat02 <- dat02[phy02$tip.label,]

#prepare loop
iterations <- 100 #set to 100 for final run

#data frame for results
aicc.scores <- data.frame(matrix(NA, nrow = iterations, ncol = 5))
colnames(aicc.scores) <- c("BM", "OU", "EB","trend", "drift")

akaike.weights <- data.frame(matrix(NA, nrow = iterations, ncol = 5))
colnames(akaike.weights) <- c("BM", "OU", "EB","trend", "drift") 

BM_parameters <- data.frame(matrix(NA, nrow = iterations, ncol = 3))
colnames(BM_parameters) <- c("sigsq", "z0", "AICc")

OU_parameters <- data.frame(matrix(NA, nrow = iterations, ncol = 4))
colnames(OU_parameters) <- c("alpha", "sigsq", "z0", "AICc")

EB_parameters <- data.frame(matrix(NA, nrow = iterations, ncol = 4))
colnames(EB_parameters) <- c("a", "sigsq", "z0", "AICc")

trend_parameters <- data.frame(matrix(NA, nrow = iterations, ncol = 4))
colnames(trend_parameters) <- c("slope", "sigsq", "z0", "AICc")

drift_parameters <- data.frame(matrix(NA, nrow = iterations, ncol = 4))
colnames(drift_parameters) <- c("drift", "sigsq", "z0", "AICc")

#model fitting loop
#for each iterations, the following steps are performed
#resample data to match N of ichthyosaurs
#match data and phylogeny
#define trait vector (skull)
#time-calibrate phylogeny and normalize branch length (tree)
#fit models

for (i in 1: iterations)try({
  
  dat02_resampled <- dat02[sample(nrow(dat02), 48), ]
  
  #match data and tree
  comp <- match.phylo.data(phy02, dat02_resampled)
  phy03 <- comp$phy
  dat03 <- comp$data
  dat03[,2:5] <- apply(dat03[,2:5], 2, as.numeric)
  dat03 <- dat03[phy03$tip.label,]
  ages <- cbind(FAD=dat03$FAD, LAD=dat03$LAD)
  rownames(ages) <- dat03$Taxon
  
  #define trait vector
  skull <- log10(dat03$Skull.width..cm.)
  names(skull) <- rownames(dat03)
  
  #time calibrate a single tree and normalize branch length
  tree <- timePaleoPhy(phy03, ages, type = "equal", vartime = 1, randres = T, ntrees=1, dateTreatment= "minMax")
  tree$edge.length <- tree$edge.length/max(nodeHeights(tree))
  
  #fit models
  brownianModel <- fitContinuous(tree, skull, model="BM")
  OUModel <- fitContinuous(tree, skull, model="OU", bounds=list(alpha=c(0,exp(6)))) #increased bounds
  EBModel <- fitContinuous(tree, skull, model="EB", bounds=list(a=c(log(10^(-5))/max(node.depth.edgelength(tree)),-1e-16))) #increased bounds
  trendModel <- fitContinuous(tree, skull, model="trend", bounds=list(slope = c(min = -100, max = 2000))) #increased bounds
  driftModel <- fitContinuous(tree, skull, model="drift")
  
  # calculate AIC weights
  bmAICC <- brownianModel$opt$aicc
  ouAICC <- OUModel$opt$aicc
  ebAICC <- EBModel$opt$aicc
  trendAICC <- trendModel$opt$aicc
  driftAICC <- driftModel$opt$aicc
  
  aicc <- c(bmAICC, ouAICC, ebAICC, trendAICC, driftAICC)
  aiccD <- aicc - min(aicc)
  aw <- exp(-0.5 * aiccD)
  aiccW <- aw/sum(aw)
  
  aicc.scores$BM[i] <- aicc[1]
  aicc.scores$OU[i] <- aicc[2]
  aicc.scores$EB[i] <- aicc[3]
  aicc.scores$trend[i] <- aicc[4]
  aicc.scores$drift[i] <- aicc[5]
  
  akaike.weights$BM[i] <- aiccW[1]
  akaike.weights$OU[i] <- aiccW[2]
  akaike.weights$EB[i] <- aiccW[3]
  akaike.weights$trend[i] <- aiccW[4]
  akaike.weights$drift[i] <- aiccW[5]
  
  BM_parameters$sigsq[i] <- brownianModel$opt$sigsq
  BM_parameters$z0[i] <- brownianModel$opt$z0
  BM_parameters$AICc[i] <-bmAICC
  
  OU_parameters$alpha[i] <- OUModel$opt$alpha
  OU_parameters$sigsq[i] <- OUModel$opt$sigsq
  OU_parameters$z0[i] <- OUModel$opt$z0
  OU_parameters$AICc[i] <- OUModel$opt$aicc
  
  EB_parameters$a[i] <- EBModel$opt$a
  EB_parameters$sigsq[i] <- EBModel$opt$sigsq
  EB_parameters$z0[i] <- EBModel$opt$z0
  EB_parameters$AICc[i] <- EBModel$opt$aicc
  
  trend_parameters$slope[i] <- trendModel$opt$slope
  trend_parameters$sigsq[i] <- trendModel$opt$sigsq
  trend_parameters$z0[i] <- trendModel$opt$z0
  trend_parameters$AICc[i] <- trendModel$opt$aicc
  
  drift_parameters$drift[i] <- driftModel$opt$drift
  drift_parameters$sigsq[i] <- driftModel$opt$sigsq
  drift_parameters$z0[i] <- driftModel$opt$z0
  drift_parameters$AICc[i] <- driftModel$opt$aicc
  
  print(paste("iteration", i, "of", iterations, sep=" "))
  
})

write.csv(aicc.scores, "aicc_cetaceans_nofilter_resampled.csv")
write.csv(akaike.weights, "aw_cetaceans_nofilter_resampled.csv")
write.csv(BM_parameters, "bm_cetaceans_nofilter_resampled.csv")
write.csv(OU_parameters, "OU_cetaceans_nofilter_resampled.csv")
write.csv(EB_parameters, "EB_cetaceans_nofilter_resampled.csv")
write.csv(trend_parameters, "trend_cetaceans_nofilter_resampled.csv")
write.csv(drift_parameters, "drift_cetaceans_nofilter_resampled.csv")

#akaike.weights <- read.csv("aw_cetaceans_nofilter_resampled.csv") #in case one later returns to the model fitting output
#akaike.weights <- akaike.weights[,2:6] #only needed if one returns to the saved csv files later!

density.data <- melt(akaike.weights)
pdf("Fig_S11_H_Akaike_weights_cetaceans_nofilterfeeders_resampled_100iterations.pdf", 8, 6)
ggplot(density.data, aes(x=value, fill=variable)) + geom_boxplot(alpha=0.25)
dev.off()
