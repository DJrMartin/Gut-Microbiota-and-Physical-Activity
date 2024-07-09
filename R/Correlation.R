#################
## CORRELATION ##
#################
rm(list=ls())
## IMPORTATIONS========================================
library(vegan)

Gut_Microbiota_Composition <- read.csv("~/Desktop/Exomic_project/data/Humans/Gut_Microbiota_Composition.csv", row.names = 1)
Metabolites_data <- read.csv("~/Desktop/Exomic_project/data/Humans/Metabolites_data.csv", row.names = 1)
Physiological_data <- read.csv("~/Desktop/Exomic_project/data/Humans/Physiological_data.csv", row.names = 1)

## NORMALISATION====================================================
### Gut Microbiota Composition
Gut_Microbiota_Composition.TSS <- as.data.frame(t(apply(Gut_Microbiota_Composition, 1, function(x) x/sum(x))))
rm(Gut_Microbiota_Composition)
### Fecal Metabolites Concentration
Metabolites_data.NORMALIZED=as.data.frame(apply(Metabolites_data[,-1], 2, function(x) x*(1-(Metabolites_data$H2O..../100))))
rm(Metabolites_data)
### Points depending on Experimental Groups and Donors
PCH = as.numeric(as.character(factor(Physiological_data$Groupe, c("CTL", "FOOT", "BIKE"), c(15, 16, 17))))

## INTRA-SAMPLES ANALYSIS===========================================
Richness <- apply(Gut_Microbiota_Composition.TSS, 1, function(x) sum(x>0))
Shannon <- diversity(Gut_Microbiota_Composition.TSS)

## Data frame with the physiological variables (A)
data_for_correlation <- scale(data.frame("VO2max"=Physiological_data$VO2_max, "EEE"= (Physiological_data$Lipides_SV1+Physiological_data$Glucides_SV1),
                           "FO"=(Physiological_data$Lipides_SV1/(Physiological_data$Lipides_SV1+Physiological_data$Glucides_SV1)),
                           "Fat Mass" = Physiological_data$Fat_Mass, Shannon, Richness))

## Data frame with the fecal SCFA content (B)
SCFA <- scale(Metabolites_data.NORMALIZED[,c(47:dim(Metabolites_data.NORMALIZED)[2])])
colnames(SCFA) = c("Acetate", "Proprionate", "Isobutyrate", "Butyrate", "Isovalerate", "Valerate")

## Creation of matrix where the results will record. 
p.value = r = matrix(NA, nrow = ncol(SCFA), ncol= ncol(data_for_correlation))
rownames(r)=rownames(p.value)=colnames(data_for_correlation)
colnames(r)=colnames(p.value)=colnames(SCFA)

## Estimation of the pearson'r, R squared and p value of the correlation between both data frame (A and B)
for (i in 1:ncol(SCFA)){
  for (j in 1:ncol(data_for_correlation)){
    if (shapiro.test(SCFA[,i])$p.value>0.05&shapiro.test(data_for_correlation[,j])$p.value>0.05){
      r[j,i] = cor.test(SCFA[,i], data_for_correlation[,j])$estimate ## Pearson'r estimation
      p.value[j,i] = cor.test(SCFA[,i], data_for_correlation[,j])$p.value ## p.value estimation
    }else{
      r[j,i] = cor.test(SCFA[,i], data_for_correlation[,j], method='spearman')$estimate ## Spearman'r estimation
      p.value[j,i] = cor.test(SCFA[,i], data_for_correlation[,j], method='spearman')$p.value ## p.value estimation
    }
  }
  print(paste("Done:", colnames(SCFA)[i]))
}

## Parameters of the plot 
### colors ramp of the r correlation direction
colfunc <- colorRampPalette(c('firebrick',"beige",'cornflowerblue'))(20)
colors_r = matrix("NA", nrow(r), ncol(r))
for (i in 1:20){
  colors_r[which(r>seq(-1, 1,by=0.1)[i]&r<seq(-1, 1,by=0.1)[i+1])]=colfunc[i]
}

### PLOT the correlogram
par(mar=c(17,10,3,1))
ord <- rep(1:ncol(SCFA), each=ncol(data_for_correlation))
abs <- rep(1:ncol(data_for_correlation),ncol(SCFA)) 
plot(abs, ord, 
     axes=F,xlab='', ylab="",xlim=c(0.5, ncol(data_for_correlation)+.5),
     cex = 3,ylim=c(0,ncol(SCFA)+1), pch=15, col=colors_r) 
      ## The weights of the points are proportionnal to the R squared.
#### Axis.
box()
axis(1, at= 1:ncol(data_for_correlation), labels = FALSE)
text(1:ncol(data_for_correlation),rep(-2, ncol(data_for_correlation)), colnames(data_for_correlation),
     srt = 30, xpd=NA, adj=1, cex = 0.8, col=c(rep('#88CCEE', 4), rep("#117733",3)))
axis(2, at= 1:6, labels = FALSE)
text(-.3, 1:ncol(SCFA),colnames(SCFA) ,
     xpd=NA, srt = 0,cex=0.8,adj=1, col="#DDCC77")
#### p value.
text(abs[which(p.value<0.05&p.value>0.01)], ord[which(p.value<0.05&p.value>0.01)],"*", cex=1, col='black')
text(abs[which(p.value<0.01)], ord[which(p.value<0.01)],"**", cex=1, col='black')
