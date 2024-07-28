#Cetaceans
#Model fitting
#Fig_S11_J

#Re-sampling approach no. 8 to more thoroughly investigate underlying patterns of trait evolution
#While all show quatic adaptions, the degree of aquatic lifestyle may differ between early ichthyosaurs and early whales.

#Do results change when only the Pelagiceti are incluced in the analysis?

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

#load and prep data
cet.data00 <- read.csv("cetacea_data.csv")
head(cet.data00)
rownames(cet.data00) <- cet.data00$Taxon

#find node that defines Pelagiceti
Pelagiceti <- findMRCA(cphy00, tips=c("Durodon_atrox",	"Stenella_longirostris"))

#get the tree for Pelagiceti
phy00_pel <- extract.clade(cphy00, Pelagiceti) 

#match data and tree
comp <- match.phylo.data(phy00_pel, cet.data00)
phy02 <- comp$phy
dat02 <- comp$data
dat02[,2:5] <- apply(dat02[,2:5], 2, as.numeric)
dat02 <- dat02[phy02$tip.label,]

#age bins for time calibration
ages <- cbind(FAD=dat02$FAD, LAD=dat02$LAD)
rownames(ages) <- dat02$Taxon

#define vector with data for model fitting
skull <- log10(dat02$Skull.width..cm.)
names(skull) <- rownames(dat02)

#generate set of time calibrated trees and normalize branch lengths
ntrees <- 100 #set to 100 for final analysis
trees <- timePaleoPhy(phy02, ages, type = "equal", vartime = 1, randres = T, ntrees=ntrees, dateTreatment= "minMax")
for (i in 1:length(trees)) {trees[[i]]$edge.length <- trees[[i]]$edge.length/max(nodeHeights(trees[[1]]))}

#prepare loop
iterations <- length(trees)

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

#perform loop: fit models over entire set of trees, only with extinct taxa
#trait skull
#trees

for (i in 1: iterations)try({
  
  brownianModel <- fitContinuous(trees[[i]], skull, model="BM") #default bounds seem fine
  OUModel <- fitContinuous(trees[[i]], skull, model="OU")
  EBModel <- fitContinuous(trees[[i]], skull, model="EB")
  trendModel <- fitContinuous(trees[[i]], skull, model="trend")
  driftModel <- fitContinuous(trees[[i]], skull, model="drift")
  
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

write.csv(aicc.scores, "aicc_cetaceans_pelagiceti.csv")
write.csv(akaike.weights, "aw_cetaceans_pelagiceti.csv")
write.csv(BM_parameters, "bm_cetaceans_pelagiceti.csv")
write.csv(OU_parameters, "OU_cetaceans_pelagiceti.csv")
write.csv(EB_parameters, "EB_cetaceans_pelagiceti.csv")
write.csv(trend_parameters, "trend_cetaceans_pelagiceti.csv")
write.csv(drift_parameters, "drift_cetaceans_pelagiceti.csv")

#akaike.weights <- read.csv("aw_cetaceans_pelagiceti.csv") #in case one later returns to the model fitting output
#akaike.weights <- akaike.weights[,2:6] #only needed if one returns to the saved csv files later!

density.data <- melt(akaike.weights)
pdf("Fig_S11_J_Akaike_weights_cetaceans_100_pelagiceti.pdf", 8, 6)
ggplot(density.data, aes(x=value, fill=variable)) + geom_boxplot(alpha=0.25)
dev.off()



