##########################
## Constrained Analysis ##
##########################
## Explanatory variables: Nutrient Intakes
## Explained variables: Gut Microbiota Compositions

rm(list=ls())
## IMPORTATIONS========================================
library(vegan)

Gut_Microbiota_Composition <- read.csv("Gut_Microbiota_Composition.csv", row.names = 1)
Physiological_data <- read.csv("Physiological_data.csv", row.names = 1)
Dietary_Food_Intakes <- read.csv("Dietary_Food_Intakes.csv", row.names = 1)
Dietary_Nutrient_Intakes <- read.csv("Dietary_Nutrient_Intakes.csv", row.names = 1)

## NORMALISATION====================================================
### Gut Microbiota Composition
Gut_Microbiota_Composition.TSS <- as.data.frame(t(apply(Gut_Microbiota_Composition, 1, function(x) x/sum(x))))
rm(Gut_Microbiota_Composition)
### Oligo-trace Intakes
Dietary_Nutrient_Intakes.OLIGO.TSS <- as.data.frame(t(apply(Dietary_Nutrient_Intakes[,1:14],1, function(x) x/sum(x))))
### Points depending on Experimental Groups and Donors
PCH = as.numeric(as.character(factor(Physiological_data$Groupe, c("CTL", "FOOT", "BIKE"), c(15, 16, 17))))

## CAP =========================================================
par(mar=c(4,4,2,2))
### Oligo-Trace
res.cap <- capscale(Gut_Microbiota_Composition.TSS~., data = Dietary_Nutrient_Intakes.OLIGO.TSS)
res.R <- RsquareAdj(res.cap) ## Computation of the R squared adjusted. 
res.anova <- anova.cca(res.cap, step = 1000) ## Computation of the Anova test. 
plot(res.cap)
paste("R =",round(res.R$adj.r.squared, 2),", p =",round(res.anova$`Pr(>F)`[1],2))
legend("topright", expression(paste(R^2, "= 0.02, p = 0.42")), bty="n")
