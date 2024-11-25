#########################
## NORMALISATION.      ##
#########################
rm(list=ls())
## IMPORTATIONS========================================
library(vegan)
library(ape)
library(dendextend)

Gut_Microbiota_Composition <- read.csv("Gut_Microbiota_Composition.csv", row.names = 1)
Gut_Microbiota_Functionality <- read.csv("Gut_Microbiota_Functionality.csv", row.names = 1)
Physiological_data <- read.csv("Physiological_data.csv", row.names = 1)
Metabolites_data <- read.csv("Metabolites_data.csv", row.names = 1)
Dietary_Food_Intakes <- read.csv("Dietary_Food_Intakes.csv", row.names = 1)
Dietary_Nutrient_Intakes <- read.csv("Dietary_Nutrient_Intakes.csv", row.names = 1)
Species_labelled <- read.csv("Species_Phylo_Label.csv", row.names = 1)
read.tree(filt_tree,"Phylogenetic_Tree.txt")

## NORMALISATION====================================================
### Gut Microbiota Composition
Gut_Microbiota_Composition.TSS <- as.data.frame(t(apply(Gut_Microbiota_Composition, 1, function(x) x/sum(x))))
rm(Gut_Microbiota_Composition)
### Fecal Metabolites Concentration
Metabolites_data.NORMALIZED=as.data.frame(apply(Metabolites_data[,-1], 2, function(x) x*(1-(Metabolites_data$H2O..../100))))
rm(Metabolites_data)
### Dietary Nutrien Intakes
Dietary_Nutrient_Intakes.OLIGO.TSS <- as.data.frame(t(apply(Dietary_Nutrient_Intakes[,1:14],1, function(x) x/sum(x))))
Dietary_Nutrient_Intakes.MACRO.TSS <- as.data.frame(t(apply(Dietary_Nutrient_Intakes[,15:22],1, function(x) x/sum(x))))
Dietary_Nutrient_Intakes.VIT.TSS <- as.data.frame(t(apply(Dietary_Nutrient_Intakes[,23:31],1, function(x) x/sum(x))))
rm(Dietary_Nutrient_Intakes)
### Points depending on Experimental Groups and Donors
PCH = as.numeric(as.character(factor(Physiological_data$Groupe, c("CTL", "FOOT", "BIKE"), c(15, 16, 17))))
