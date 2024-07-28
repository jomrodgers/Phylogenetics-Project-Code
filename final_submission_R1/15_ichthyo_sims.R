#Simulate traits under BM and EB model for ichthyosaur tree
#Test whether the method is able to detect the model under which the traits were simulated
#Fig S11 K and L

#time-calibrate tree
#simulate trait
#fit models

rm(list=ls())

setwd("C:/Users/Lars Schmitz/Dropbox/@JIM2_manuscript_PCM/final_submission_R1")

#load packages
library (phytools)
library (geiger)
library (paleotree)
library (reshape2)
library (ggplot2)
library (mvMORPH)
library (picante)

#load tree
phy00 <- read.tree("ichthyosaur_tree.tre")
phy00$root.time <- 251.902

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

#load model parameter estimates from empirical runs
bm_p <- read.csv("BM_ichthyos.csv")
    bm_theta <- mean(bm_p$z0)
    bm_sigsq <- mean(bm_p$sigsq)
eb_p <- read.csv("EB_ichthyos.csv")
    eb_theta <- mean(eb_p$z0)
    eb_sigsq <- mean(eb_p$sigsq)
    eb_ratechange <- mean(eb_p$a)

#number of iterations
iterations <- 100 #set to 100 for final analysis

###BM

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

for (i in 1:iterations){
  
  #time calibrate tree
  tree <- timePaleoPhy(itre0, iages, type = "equal", vartime = 1, randres = T, 
                        ntrees=1, dateTreatment= "minMax")
  
  #simulate data with BM
  trait <- mvSIM(tree, nsim = 1, error = NULL, model = c("BM1"),
        param = list(theta = bm_theta, sigma = bm_sigsq))
  
  #normalize branch lengths
  tree$edge.length <- tree$edge.length/max(nodeHeights(tree))
  
  #fit models
  brownianModel <- fitContinuous(tree, trait, model="BM")
  OUModel <- fitContinuous(tree, trait, model="OU", bounds=list(alpha=c(0,exp(5)))) #increased bounds for max(alpha) after running up against the upper bound (problem remains)
  EBModel <- fitContinuous(tree, trait, model="EB")
  trendModel <- fitContinuous(tree, trait, model="trend")
  driftModel <- fitContinuous(tree, trait, model="drift")
  
  #calculate AIC weights
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
}

#save results
write.csv(aicc.scores, "aicc_sim_ichthyos_BM.csv")
write.csv(akaike.weights, "aw_sim_ichthyos_BM.csv")
write.csv(BM_parameters, "bm_sim_ichthyos_BM.csv")
write.csv(OU_parameters, "OU_sim_ichthyos_BM.csv")
write.csv(EB_parameters, "EB_sim_ichthyos_BM.csv")
write.csv(trend_parameters, "trend_sim_ichthyos_BM.csv")
write.csv(drift_parameters, "drift_sim_ichthyos_BM.csv")

#akaike.weights <- read.csv("aw_sim_ichthyos_BM.csv") #in case one later returns to the model fitting output
#akaike.weights <- akaike.weights[,2:6] #only needed when reading previous results from csv file

#graph results
density.data <- melt(akaike.weights)
pdf("Fig_S11_K_Akaike_weights_sim_ichthyos_BM.pdf", 8, 6)
ggplot(density.data, aes(x=value, fill=variable)) + geom_boxplot(alpha=0.25)
dev.off()


###########################

#EB

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

for (i in 1:iterations){
  
  #time calibrate tree
  tree <- timePaleoPhy(itre0, iages, type = "equal", vartime = 1, randres = T, 
                       ntrees=1, dateTreatment= "minMax")
  #simulate data with EB
  trait <- mvSIM(tree, nsim = 1, error = NULL, model = c("EB"),
                 param = list(theta = eb_theta, sigma = eb_sigsq, beta = eb_ratechange))
  
  #normalize branch lengths
  tree$edge.length <- tree$edge.length/max(nodeHeights(tree))
  
  #fit models
  brownianModel <- fitContinuous(tree, trait, model="BM")
  OUModel <- fitContinuous(tree, trait, model="OU", bounds=list(alpha=c(0,exp(5)))) #increased bounds for max(alpha) after running up against the upper bound (problem remains)
  EBModel <- fitContinuous(tree, trait, model="EB")
  trendModel <- fitContinuous(tree, trait, model="trend")
  driftModel <- fitContinuous(tree, trait, model="drift")
  
  #calculate AIC weights
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
}

#save results
write.csv(aicc.scores, "aicc_sim_ichthyos_EB.csv")
write.csv(akaike.weights, "aw_sim_ichthyos_EB.csv")
write.csv(BM_parameters, "bm_sim_ichthyos_EB.csv")
write.csv(OU_parameters, "OU_sim_ichthyos_EB.csv")
write.csv(EB_parameters, "EB_sim_ichthyos_EB.csv")
write.csv(trend_parameters, "trend_sim_ichthyos_EB.csv")
write.csv(drift_parameters, "drift_sim_ichthyos_EB.csv")

#akaike.weights <- read.csv("aw_sim_ichthyos_EB.csv") #in case one later returns to the model fitting output
#akaike.weights <- akaike.weights[,2:6] #only needed when reading previous results from csv file

#graph results
density.data <- melt(akaike.weights)
pdf("Fig_S11_L_Akaike_weights_sim_ichthyos_EB.pdf", 8, 6)
ggplot(density.data, aes(x=value, fill=variable)) + geom_boxplot(alpha=0.25)
dev.off()


###########################