#Cetaceans
#Exploring time calibration results in more depth
#Extracting node age estimates

rm(list=ls())

setwd("C:/Users/Lars Schmitz/Dropbox/@JIM2_manuscript_PCM/final_submission_R1")

library(phytools)
library(ape)
library(ggplot2)
library(reshape2)

#getting tree set (n=1,000)
#either run the line below or read in the treeset we supply with the script
#trees <- timePaleoPhy(phy02, ages, type = "equal", vartime = 1, randres = T, ntrees=1000, dateTreatment= "minMax")
#write.tree(trees, "cetacea_trees.tre")

#read in time calibrated trees generated previously
trees <- read.tree("cetacea_trees.tre")


ages <- data.frame(matrix(NA, nrow = length(trees), ncol = 16))
colnames(ages) <- c("Pelagiceti",
                    "Neoceti",
                    "Mysticeti",
                    "Chaeomysticeti",
                    "Crown_Mysticeti",
                    "Balaenidae",
                    "Plicogulae",
                    "Balaenopteridae",
                    "Eschrichtiinae",
                    "Odontoceti",
                    "Crown_Odontoceti",
                    "Pan_Physeterioidea",
                    "Ziphiidae",
                    "Delphinida", 
                    "Crown_Delphinida",
                    "Delphinidae")

#findMRCA

for (i in 1:length(trees)){
tree <- trees[[i]]
Pelagiceti <- max(nodeHeights(tree)) - findMRCA(tree, 
                       tips=c("Durodon_atrox",	"Stenella_longirostris"),
                       type="height")
Neoceti <- max(nodeHeights(tree)) - findMRCA(tree, 
                      tips=c("Coronodon_havensteini",	"Stenella_longirostris"),
                      type="height")
Mysticeti <- max(nodeHeights(tree)) - findMRCA(tree, 
                      tips=c("Coronodon_havensteini",	"Balaenoptera_physalus"),
                      type="height")
Chaeomysticeti <- max(nodeHeights(tree)) - findMRCA(tree, 
                      tips=c("Maiabalaena_nesbittae",	"Balaenoptera_physalus"),
                      type="height")
Crown_Mysticeti <- max(nodeHeights(tree)) - findMRCA(tree, 
                      tips=c("Cophocetus_oregonensis",	"Balaenoptera_physalus"),
                      type="height")
Balaenidae <- max(nodeHeights(tree)) - findMRCA(tree, 
                      tips=c("Peripolocetus_vexillifer",	"Balaena_mysticetus"),
                      type="height")
Plicogulae <- max(nodeHeights(tree)) - findMRCA(tree, 
                      tips=c("Isanacetus_laticephalus",	"Balaenoptera_physalus"),
                      type="height")
Balaenopteridae <- max(nodeHeights(tree)) - findMRCA(tree, 
                      tips=c("Norrisanima_miocaena",	"Balaenoptera_physalus"),
                      type="height")
Eschrichtiinae <- max(nodeHeights(tree)) - findMRCA(tree, 
                      tips=c("Gricetoides_aurorae",	"Eschrichtius_robustus"),
                      type="height")
Odontoceti <- max(nodeHeights(tree)) - findMRCA(tree, 
                      tips=c("Olympicetus_avitus",	"Stenella_longirostris"),
                      type="height")
Crown_Odontoceti <- max(nodeHeights(tree)) - findMRCA(tree, 
                      tips=c("Urkudelphis_chawpipacha",	"Stenella_longirostris"),
                      type="height")
Pan_Physeterioidea <- max(nodeHeights(tree)) - findMRCA(tree, 
                      tips=c("Eudelphis_mortezelensis",	"Kogia_sima"),
                      type="height")
Ziphiidae <- max(nodeHeights(tree)) - findMRCA(tree, 
                      tips=c("Chavinziphius_maxillocristatus",	"Mesoplodon_peruvianus"),
                      type="height")
Delphinida <- max(nodeHeights(tree)) - findMRCA(tree, 
                      tips=c("Delphinodon_dividum",	"Stenella_longirostris"),
                      type="height")
Crown_Delphinida <- max(nodeHeights(tree)) - findMRCA(tree, 
                      tips=c("Inia_geoffrensis",	"Stenella_longirostris"),
                      type="height")
Delphinidae <- max(nodeHeights(tree)) - findMRCA(tree, 
                      tips=c("Leucopleurus_acutus",	"Stenella_longirostris"),
                      type="height")

ages$Pelagiceti[i] <- Pelagiceti
ages$Neoceti[i] <- Neoceti
ages$Mysticeti[i] <- Mysticeti
ages$Chaeomysticeti[i] <- Chaeomysticeti
ages$Crown_Mysticeti[i] <- Crown_Mysticeti
ages$Balaenidae[i] <- Balaenidae
ages$Plicogulae[i] <- Plicogulae
ages$Balaenopteridae[i] <- Balaenopteridae
ages$Eschrichtiinae[i] <- Eschrichtiinae
ages$Odontoceti[i] <- Odontoceti
ages$Crown_Odontoceti[i] <- Crown_Odontoceti
ages$Pan_Physeterioidea[i] <- Pan_Physeterioidea
ages$Ziphiidae[i] <- Ziphiidae
ages$Delphinida[i] <- Delphinida 
ages$Crown_Delphinida[i] <- Crown_Delphinida
ages$Delphinidae[i] <- Delphinidae
}


write.csv(ages, "cetacean_node_ages.csv")

age.data <- melt(ages)
pdf("cetacean_node_ages.pdf", 11, 8.5)
ggplot(age.data, aes(x=-value, fill=variable)) + 
  geom_boxplot(alpha=0.25) + 
  coord_flip() +
  xlab("million years")
dev.off()

summarytable <- as.data.frame(summary(ages))
write.csv(summarytable, "cetacean_node_age_summary.csv")


