#Cetaceans
#Cetacean timetree with taxon labels, upright for best overall readability
#Final edits, incl. node names and design improvements in Illustrator by Jorge

rm(list=ls())

setwd("C:/Users/jmr95/OneDrive/Documents/AP&P_data/science.abf5787_data_s6/final_submission_R1")

#load packages
library (phytools)
library (strap)
library (geiger)
library (paleotree)
library (picante)

#load tree
cphy00 <- read.nexus("cetacean_tree_0529.nex")
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

### time tree cetaceans
phy <- DatePhylo(cphy01, cages, method="equal", rlen=1, add.terminal = T)

pdf("Fig_S8_upright_cet_time_tree_with_taxon_namessssss.pdf", 11, 8.5)
geoscalePhylo(tree=phy, ages=cages, units=c("Epoch", "Age"), boxes=c("Epoch", "Age"), 
              direction="upwards", cex.ts=1.2, tick.scale="yes", label.offset=0, show.tip.label=T)
dev.off()



#Note: uncertain why gray timezones do not extend all the way across?






######################### A few changes to make the pdf look better

rm(list=ls())

# Set working directory (Update the path to your specific directory)
setwd("C:/Users/jmr95/OneDrive/Documents/AP&P_data/science.abf5787_data_s6/final_submission_R1")

# Load packages
library(phytools)
library(strap)
library(geiger)
library(paleotree)
library(picante)

# Load tree
cphy00 <- read.nexus("cetacean_tree_0529.nex")
cphy00$root.time <- 56

# Load and prep data
cet.data00 <- read.csv("C:/Users/jmr95/OneDrive/Documents/AP&P_data/cetacea_data.csv")
head(cet.data00)
rownames(cet.data00) <- cet.data00$Taxon

# Match data and tree
ccomp <- match.phylo.data(cphy00, cet.data00)
cphy01 <- ccomp$phy
cdat01 <- ccomp$data
cdat01[,2:5] <- apply(cdat01[,2:5], 2, as.numeric)
cdat01 <- cdat01[cphy01$tip.label,]
cages <- cbind(FAD=cdat01$FAD, LAD=cdat01$LAD)
rownames(cages) <- cdat01$Taxon

# Time tree cetaceans
phy <- DatePhylo(cphy01, cages, method="equal", rlen=1, add.terminal = T)

# Define the PDF output file path
pdf_file_path <- "Fig_S8_upright_cet_time_tree_pakicetus_included.pdf"

# Check if the file can be created
if (file.access(dirname(pdf_file_path), 2) == 0) {
  # Open PDF device
  pdf(pdf_file_path, width = 11, height = 8.5)
  # Generate the plot with geological units and Ma
  geoscalePhylo(tree=phy, ages=cages, units=c("Period", "Epoch", "Age"), 
                boxes="auto", cex.age=0.8, cex.ts=1.2, tick.scale=1, label.offset=0, show.tip.label=T)
  # Close PDF device
  dev.off()
  cat("PDF file created successfully: ", pdf_file_path, "\n")
} else {
  cat("Cannot write to directory: ", dirname(pdf_file_path), "\n")
}

