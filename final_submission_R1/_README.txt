This folder contains data and script to perform analyses and plots for the "Jim2-Project"
Sander et al.

Please install the latest version of R and install the following packages:
	ape
	bayou
	dplyr #used during initial analyses
	geiger
	googlesheets4 #used during initial analyses
	ggplot2
	mvMORPH
	paleotree
	phytools
	picante
	reshape2
	strap
	tidyverse #used during initial analyses

Please also make sure to set your working directory correctly.

Required files (must be in working directory)
	DATA
		ichthyosaur_data.csv
		cetacea_data.csv
		tooth_size-csv
	TREES
		ichthyosaur_tree.tre
		cetacean_tree_0529.nex
	FUNCTIONS
		aa_generics_v1.R
		tnt_R_utils_v1.R
		multirate_plot.R
	
Required subfolders in working directory
	ibayou_results #for ichthyosaurs
	cbayou_results #for cetaceans

Brief explanation of R files
	01 time calibration of ichthyosaur and cetacean tree [for Figure 3]
	02 ichthyosaur topology with Bremer values [for Figure S8]
	03 cetacean timetree with taxon labels, upright [for Figure S9]
	04 cetacean timetree age estimates [for comparison with other studies]
	05 ichthyosaur model fitting [produces panel for Figures 4, S10, and S11 panel A]
	06 cetacean model fitting [produces panel for Figures 4, and S11 panel B]
	07 cetacean model fitting, resampling 1 [produces panel for Figure S11 C]
	08 cetacean model fitting, resampling 2 [produces panel for Figure S11 D]
	09 cetacean model fitting, resampling 3 [produces panel for Figure S11 E]
	10 cetacean model fitting, resampling 4 [produces panel for Figure S11 F]
	11 cetacean model fitting, resampling 5 [produces panel for Figure S11 G]
	12 cetacean model fitting, resampling 6 [produces panel for Figure S11 H]
	13 cetacean model fitting, resampling 7 [produces panel for Figure S11 I]
	14 cetacean model fitting, resampling 8 [produces panel for Figure S11 J]
	15 ichthyosaur simulations [produces two panels for Figure S11, K and L]
	16 combined traitgram [for Figure 4]
	17 traitgram for ichthyosaurs with taxon labels [for Figure S10]
	18 ichthyosaur tooth size comparison, conventional stats [for Figure S5]
	19 ancestral state estimates using fitContinuous
	20 rate heterogeneity analysis, multirateBM (ichthyosaurs) [for Figure 5a]
	21 rate heterogeneity analysis, multirateBM (cetaceans) [for Figure 5b]
	22 dynamics of the adaptive landscape, bayou (ichthyosaurs) [for Figure 6a]
	23 dynamics of the adaptive landscape, bayou (cetaceans) [for Figure 6b]
	24 barplots of ichthyosaur skull length
	25 barplots of cetacean skull width

