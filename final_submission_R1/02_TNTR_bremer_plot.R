#Ichthyosaurs
#Reading in Bremer tree obtained from TNT (Nicole), converting to Newick file (script from Nick Matzke), and plot
#Final edits in Illustrator (Fig S7)

rm(list=ls())

setwd("C:/Users/Lars Schmitz/Dropbox/@JIM2_manuscript_PCM/final_submission_R1")

library(ape)
source("tnt_R_utils_v1.R") #Nick Matzke
source("aa_generics_v1.R") #Nick Matzke

tntfn <- "ichthyo_bremer.nex"
trfn <- tntfile2newick(tntfn, brlens=FALSE, options="text_convert", translate=FALSE, branchlabels="=")
moref(trfn)
tre <- read.tree(trfn)
tre$node.label
tre <- ladderize(tre)

pdf("Fig_S7_ichthyosaur_bremer_support_circle_09_19.pdf", 8.5, 11) #needs edits in Illustrator, cricles not everyone's favorite
plot(tre, lwd=2, cex=0.75)
nodelabels(tre$node.label, frame="circle", bg="white", cex=0.35)
dev.off()



