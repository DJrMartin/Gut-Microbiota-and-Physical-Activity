##################
## Linear Model ##
##################
rm(list=ls())
## IMPORTATIONS========================================
library(vegan)
library(ape)

Gut_Microbiota_Composition <- read.csv("Gut_Microbiota_Composition.csv", row.names = 1)
Physiological_data <- read.csv("Physiological_data.csv", row.names = 1)

## NORMALISATION====================================================
### Gut Microbiota Composition
Gut_Microbiota_Composition.TSS <- as.data.frame(t(apply(Gut_Microbiota_Composition, 1, function(x) x/sum(x))))
rm(Gut_Microbiota_Composition)

### Coordinates from Gut Microbiota Pairwise Dissimalirity Matrix
Coordinates <- data.frame(pcoa(vegdist(Gut_Microbiota_Composition.TSS))$vectors)[,1:2]

### Data frame of physiological data
df <- data.frame(VO2 =Physiological_data$VO2_max, 
                 Fat.Ox = Physiological_data$Lipides_SV1/(Physiological_data$Lipides_SV1+Physiological_data$Glucides_SV1),
                 EEE = Physiological_data$Lipides_SV1+Physiological_data$Glucides_SV1,
                 FM = Physiological_data$Fat_Mass, BMI = Physiological_data$BMI)
### BOOSTRAP
p = ncol(df)

B = 1000
PC1_exp = PC2_exp = NULL
for (i in 1:p){
  coef_lm_PC1=coef_lm_PC2=NULL
  for (cnt in 1:B){
    boostrap=sample(1:nrow(df), nrow(df), replace = T)
    coef_lm_PC1 = c(coef_lm_PC1,cor.test(Coordinates$Axis.1[boostrap],df[boostrap,i], method="spearman")$estimate)
    coef_lm_PC2 = c(coef_lm_PC2,cor.test(Coordinates$Axis.2[boostrap],df[boostrap,i], method="spearman")$estimate)
  }
  PC1_exp = cbind(PC1_exp,coef_lm_PC1)
  PC2_exp = cbind(PC2_exp,coef_lm_PC2)
}

### PLOT
par(mar=c(17,2,2,1))
plot(NA, xlim=c(0.5,p+.5), ylim=c(-1,1),
     axes=F, xlab="", ylab="")
abline(h=0, col="black", lty='dashed')
box()
axis(2)
axis(1,1:p,labels=FALSE)
abline(v=3.5)
text(1:p,rep(-1.3, p),colnames(df), 
     srt = 45, xpd=NA, adj=c(1,1), cex = 0.8)
legend("topleft", c("r with PCo1","r with PCo2"), 
       bty = 'n', fill = c("#88CCEE","#88CCEE50"), cex=0.8, title.adj = 0)
#### PCo1
##### 95%-Confident Intervals of the correlation coef. of physiological variables/body compo with PCo1
inf. = apply(PC1_exp,2,quantile,probs=.025)
sup. = apply(PC1_exp,2,quantile,probs=.975)
medians = apply(PC1_exp,2,quantile,probs=.5)
##### Segments of the CI
X.Y=data.frame(x1=0.8:(p-.2), y1=inf.,
               x2=0.8:(p-.2), y2=sup.)
X.Y=X.Y[which(is.na(X.Y[,2])==F),]
segments(X.Y$x1, X.Y$y1, X.Y$x2, X.Y$y2, lwd=2)
##### Medianes
points(x=0.8:(p-.2),y=medians, lwd=3, cex=2 , pch=18, col=c(rep('#88CCEE',3), rep("#909090",2)))

#### PCo2
##### 95%-Confident Intervals of the correlation coef. of physiological variables/body compo with PCo2
inf. = apply(PC2_exp,2,quantile,probs=.025)
sup. = apply(PC2_exp,2,quantile,probs=.975)
medians = apply(PC2_exp,2,quantile,probs=.5)
##### Segments of the CI
X.Y=data.frame(x1=1.2:(p+.2), y1=inf.,
               x2=1.2:(p+.2), y2=sup.)
X.Y=X.Y[which(is.na(X.Y[,2])==F),]
segments(X.Y$x1, X.Y$y1, X.Y$x2, X.Y$y2, lwd=2)
##### Medianes
points(x=1.2:(p+.2), y=medians, lwd=3, cex=2 , pch=18, col=c(rep('#88CCEE50',3), rep("#90909050",2)))

