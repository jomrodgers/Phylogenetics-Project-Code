#Cetaceans
#Model fitting
#Fig 4 lower insert
#Fig S11b

#Analysis of all data

rm(list=ls())

setwd("C:/Users/jmr95/OneDrive/Documents/AP&P_data")

#load packages
library (phytools)
library (geiger)
library (paleotree)
library (reshape2)
library (ggplot2)
library (picante)

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

#trait vector (skull)
skull <- log10(cdat01$Skull.width..cm.)
names(skull) <- rownames(cdat01)

##################
#generating a set of 1,000 time-calibrated trees

iterations <- 1000 #set to 1000 for final run
trees <- timePaleoPhy(cphy01, cages, type = "equal", vartime = 1, randres = T, ntrees=iterations, dateTreatment= "minMax")

#write.tree(trees, "cetacea_trees.tre")
#trees <- read.tree("cetacea_trees.tre")



##################
#normalize tree banches
for (i in 1:length(trees)) {trees[[i]]$edge.length <- trees[[i]]$edge.length/max(nodeHeights(trees[[1]]))}

#data frame for results (AKaike weights as a start)
aicc.scores <- data.frame(matrix(NA, nrow = length(trees), ncol = 5))
colnames(aicc.scores) <- c("BM", "OU", "EB","trend", "drift")

akaike.weights <- data.frame(matrix(NA, nrow = length(trees), ncol = 5))
colnames(akaike.weights) <- c("BM", "OU", "EB","trend", "drift") 

BM_parameters <- data.frame(matrix(NA, nrow = length(trees), ncol = 3))
colnames(BM_parameters) <- c("sigsq", "z0", "AICc")

OU_parameters <- data.frame(matrix(NA, nrow = length(trees), ncol = 4))
colnames(OU_parameters) <- c("alpha", "sigsq", "z0", "AICc")

EB_parameters <- data.frame(matrix(NA, nrow = length(trees), ncol = 4))
colnames(EB_parameters) <- c("a", "sigsq", "z0", "AICc")

trend_parameters <- data.frame(matrix(NA, nrow = length(trees), ncol = 4))
colnames(trend_parameters) <- c("slope", "sigsq", "z0", "AICc")

drift_parameters <- data.frame(matrix(NA, nrow = length(trees), ncol = 4))
colnames(drift_parameters) <- c("drift", "sigsq", "z0", "AICc")


for (i in 1: iterations)try({
  
  brownianModel <- fitContinuous(trees[[i]], skull, model="BM") #default bounds appear to be fine!
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

write.csv(aicc.scores, "aicc_cetaceans.csv")
write.csv(akaike.weights, "aw_cetaceans.csv")
write.csv(BM_parameters, "bm_cetaceans.csv")
write.csv(OU_parameters, "OU_cetaceans.csv")
write.csv(EB_parameters, "EB_cetaceans.csv")
write.csv(trend_parameters, "trend_cetaceans.csv")
write.csv(drift_parameters, "drift_cetaceans.csv")

#plotting Akaike weights
#next two lines read in results from previously saved runs!
# akaike.weights <- read.csv("aw_cetaceans.csv")
# akaike.weights <- akaike.weights[,2:6] #only needed when loading previous results

  mean(akaike.weights$EB); median(akaike.weights$EB)
  mean(akaike.weights$trend); median(akaike.weights$trend)
  mean(akaike.weights$drift); median(akaike.weights$drift)
  mean(akaike.weights$BM); median(akaike.weights$BM)
  mean(akaike.weights$OU); median(akaike.weights$OU)
  
density.data <- melt(akaike.weights)
pdf("Fig_4_Akaike_weights_cetaceans_1000iterations_final.pdf", 8, 6)
ggplot(density.data, aes(x=value, fill=variable)) + geom_boxplot(alpha=0.25)
dev.off()

#same figure also uswd for S11b

