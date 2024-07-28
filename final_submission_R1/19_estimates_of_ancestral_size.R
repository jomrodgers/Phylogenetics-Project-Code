#ancestral state estimates for cetacean skull iwdth and ichthyosaurian skull length
#What were the starting points?

rm(list=ls())

setwd("C:/Users/Lars Schmitz/Dropbox/@JIM2_manuscript_PCM/final_submission_R1")

#load packages
library (phytools)
library (strap)
library (geiger)
library (paleotree)
library (picante)

#CETACEANS

#load tree
cphy00 <- read.nexus("C:/Users/jmr95/OneDrive/Documents/AP&P_data/science.abf5787_data_s6/final_submission_R1/cetacean_tree_0529.nex")
cphy00$root.time <- 56

#load and prep data
cet.data00 <- read.csv("C:/Users/jmr95/OneDrive/Documents/AP&P_data/cetacea_data.csv")
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
#generating a set of time-calibrated trees

iterations <- 100 #set to 100 for final run
trees <- timePaleoPhy(cphy01, cages, type = "equal", vartime = 1, randres = T, ntrees=iterations, dateTreatment= "minMax")

write.tree(trees, "C:/Users/jmr95/OneDrive/Documents/AP&P_data/science.abf5787_data_s6/final_submission_R1/cetacea_trees.tre")
trees <- read.tree("C:/Users/jmr95/OneDrive/Documents/AP&P_data/science.abf5787_data_s6/final_submission_R1/cetacea_trees.tre")



##################
#normalize tree banches
for (i in 1:length(trees)) {trees[[i]]$edge.length <- trees[[i]]$edge.length/max(nodeHeights(trees[[1]]))}

#data frame for results
anc_skull_width <- data.frame(matrix(NA, nrow = length(trees), ncol = 1))
colnames(anc_skull_width) <- c("ancestral_skull_width")


iterations <- 100 #set to 100 for final run

for (i in 1: iterations)try({
  
  brownianModel <- fitContinuous(trees[[i]], skull, model="BM", bounds=list(sigsq=c(min=exp(-700), max=exp(200)))) #default bounds appear to be fine but expand anyway
  anc_skull_width$ancestral_skull_width[i] <- brownianModel$opt$z0
  
  print(paste("iteration", i, "of", iterations, sep=" "))
  
})


#write.csv("cetacea_anc_skull_width.csv", anc_skull_width)
#quick check of what the distribution looks
hist(anc_skull_width[,1])


#Ichthyosaurs


#load tree
phy00 <- read.tree("C:/Users/jmr95/OneDrive/Documents/AP&P_data/ichthyosaur_tree.tre")
phy00$root.time <- 251.902
plot(phy00, show.tip.label=F)

#load and prep data
ichthyo.data00 <- read.csv("C:/Users/jmr95/OneDrive/Documents/AP&P_data/ichthyosaur_data.csv")
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
iterations <- 100 #set to 100 for complete run
trees <- timePaleoPhy(itre0, iages, type = "equal", vartime = 1, randres = T, 
                      ntrees=iterations , dateTreatment= "minMax") #adjust number of trees if desired, especially for quick checks!
write.tree(trees, "C:/Users/jmr95/OneDrive/Documents/AP&P_data/ichthyosaur_tree.tre")


##################

#normalize tree banches
for (i in 1:length(trees)) {trees[[i]]$edge.length <- trees[[i]]$edge.length/max(nodeHeights(trees[[1]]))}

#data frame for results (AKaike weights as a start)
anc_skull_length<- data.frame(matrix(NA, nrow = length(trees), ncol = 1))
colnames(anc_skull_length) <- c("ancestral_skull_length") 

iterations <- length(trees)

for (i in 1: length(trees))try({
  
  EBModel <- fitContinuous(trees[[i]], skull, model="EB", bounds=list(a=c(log(10^(-5))/max(node.depth.edgelength(trees[[i]])),-1e-16))) #increased bounds
  anc_skull_length$ancestral_skull_length[i] <- EBModel$opt$z0
 
  print(paste("iteration", i, "of", iterations, sep=" "))
  
})


write.csv("ichthyosaur_ancestral_skull_length.csv", anc_skull_length)

#quick check of what the distribution looks
hist(anc_skull_length[,1])

cetacean_median_skwidth_ase <- median(anc_skull_width[,1]) #; mean(anc_skull_width[,1])
ichthyo_median_sklength_ase <- median(anc_skull_length[,1]) #; mean(anc_skull_length[,1])

10^cetacean_median_skwidth_ase
10^ichthyo_median_sklength_ase

