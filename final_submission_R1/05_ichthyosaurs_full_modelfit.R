#ichthyosaurs
#model fitting
#Fig 4 top insert
#Fig S10
#Fig S11, panel A

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
phy00 <- read.tree("ichthyosaur_tree.tre")
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
head(idat0)

#trait vector (skull)
skull <- log10(idat0$Skull.length..mm./10) #divide by 10 to match cetacean data [cm]]
names(skull) <- rownames(idat0)


#time calibration
iterations <- 1000 #set to 1000 for complete run
trees <- timePaleoPhy(itre0, iages, type = "equal", vartime = 1, randres = T, 
                      ntrees=iterations , dateTreatment= "minMax") #adjust number of trees if desired, especially for quick checks!
write.tree(trees, "ichthyosaur_trees.tre")


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

iterations <- length(trees)

for (i in 1: length(trees))try({

    brownianModel <- fitContinuous(trees[[i]], skull, model="BM")
    OUModel <- fitContinuous(trees[[i]], skull, model="OU", bounds=list(alpha=c(0,exp(3)))) #increased bounds for max(alpha) after running up against the upper bound (problem remains)
    EBModel <- fitContinuous(trees[[i]], skull, model="EB")
    trendModel <- fitContinuous(trees[[i]], skull, model="trend")
    driftModel <- fitContinuous(trees[[i]], skull, model="drift")
    
    # retrieve AICc scores and calculate Akaike weights
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

#saving results
write.csv(aicc.scores, "aicc_ichthyos.csv")
write.csv(akaike.weights, "aw_ichthyos.csv")
write.csv(BM_parameters, "bm_ichthyos.csv")
write.csv(OU_parameters, "OU_ichthyos.csv")
write.csv(EB_parameters, "EB_ichthyos.csv")
write.csv(trend_parameters, "trend_ichthyos.csv")
write.csv(drift_parameters, "drift_ichthyos.csv")

#plotting Akaike weights
#next two lines read in results from previous runs!
  #akaike.weights <- read.csv("aw_ichthyos.csv")
  #akaike.weights <- akaike.weights[,2:6]

  mean(akaike.weights$EB); median(akaike.weights$EB)
  mean(akaike.weights$trend); median(akaike.weights$trend)
  mean(akaike.weights$EB)+mean(akaike.weights$trend)
  
density.data <- melt(akaike.weights)
pdf("Figure_4_Akaike_weights_ichthyosaurs_1000iterations.pdf", 8, 6)
ggplot(density.data, aes(x=value, fill=variable)) + geom_boxplot(alpha=0.25)
dev.off()

#same plot used for S10 and S11a