#########################
## EXPLORATORY METHODS ##
#########################
rm(list=ls())
## IMPORTATIONS=====================================================
library(vegan)
library(ape)
library(dendextend)

Gut_Microbiota_Composition <- read.csv("~/Desktop/Exomic_project/data/Humans/Gut_Microbiota_Composition.csv", row.names = 1)
Physiological_data <- read.csv("~/Desktop/Exomic_project/data/Humans/Physiological_data.csv", row.names = 1)
Species_labelled <- read.csv("~/Desktop/Exomic_project/data/Humans/Species_Phylo_Label.csv", row.names = 1)
read.tree(filt_tree,"~/Desktop/Exomic_project/data/Humans/Phylogenetic_Tree.txt")

## NORMALISATION====================================================
### Total Sum Scaling
Gut_Microbiota_Composition.TSS <- as.data.frame(t(apply(Gut_Microbiota_Composition, 1, function(x) x/sum(x))))
rm(Gut_Microbiota_Composition)
### Points depending on Experimental Groups and Donors
PCH = as.numeric(as.character(factor(Physiological_data$Groupe, c("NoA", "ESP", "EC"), c(15, 16, 17))))

## INTRA-SAMPLES ANALYSIS============================================
Richness = apply(Gut_Microbiota_Composition.TSS, 1,  function(x) length(which(x>0)))
Shannon = vegan::diversity(Gut_Microbiota_Composition.TSS, 'shannon')

## INTER-SAMPLES ANALYSIS============================================
### Gut Microbiota Composition depending on Genus
Gut_Microbiota_Composition.TSS.1 <- Gut_Microbiota_Composition.TSS
colnames(Gut_Microbiota_Composition.TSS.1) <- Species_labelled$Genus

## Sum of Species belonging to the same genus.
Gut_Microbiota_Composition.TSS.GENUS=NULL
for (i in unique(Species_labelled$Genus)){
  Gut_Microbiota_Composition.TSS.GENUS <- 
    cbind(Gut_Microbiota_Composition.TSS.GENUS,
          apply(as.matrix(Gut_Microbiota_Composition.TSS.1[,which(Species_labelled$Genus==i)]), 1, sum))
}
colnames(Gut_Microbiota_Composition.TSS.GENUS)=unique(Species_labelled$Genus)
Gut_Microbiota_Composition.TSS.GENUS <- data.frame(Gut_Microbiota_Composition.TSS.GENUS)

#### Color Palette of Prevotella/Bacteroides Ratio (can be changed)
RATIO <- Gut_Microbiota_Composition.TSS.GENUS$Prevotella.s/Gut_Microbiota_Composition.TSS.GENUS$Bacteroides.s
colors.ratio <- colorRampPalette(c('#44AA99','#AA4499'))(dim(Physiological_data)[1])[rank(RATIO)]

### Principal Coordinates from the pairwise Dissimalirity Matrix (Bray-Curtis)
Coordinates <- data.frame(pcoa(vegdist(Gut_Microbiota_Composition.TSS))$vectors)[,1:2]

### Figure of the Principale Coordinates Analysis
Eigenvalues <- round(pcoa(vegdist(Gut_Microbiota_Composition.TSS))$values[1:2,3]*100)
### Correlation circle of the bacteria of the Principales Coordinates Analysis. BIPLOT
layout(matrix(c(1,2), nrow=1))
#### First plot: Samples
plot(Coordinates, xlab=paste('PCoA 1 (',Eigenvalues[1],'%)'), ylab=paste('PCoA 2 (',Eigenvalues[2],'%)'),
     col=colors.ratio, cex=1, pch= PCH)

#### Second plot: Variables
n_axes = 2
r_PC = matrix(0, ncol = n_axes , nrow = ncol(Gut_Microbiota_Composition.TSS))
for (i in 1:ncol(Gut_Microbiota_Composition.TSS)){
  r_PC[i,1] <- cor(as.numeric(Gut_Microbiota_Composition.TSS[,i]),Coordinates$Axis.1, method = 'spearman')
  r_PC[i,2] <- cor(as.numeric(Gut_Microbiota_Composition.TSS[,i]),Coordinates$Axis.2, method = 'spearman')
}

plot(1,1, xlim=c(-1.2,1.2), ylim=c(-1,1), col="white")
arrows(x0 = rep(0, ncol(Gut_Microbiota_Composition.TSS)),
       y0 = rep(0, ncol(Gut_Microbiota_Composition.TSS)), 
       x1 = r_PC[,1] , y1 = r_PC[,2], lwd=2)
w.species <- which(r_PC[,1]<(-0.5)|r_PC[,2]<(-0.4))
text(r_PC[w.species,],colnames(Gut_Microbiota_Composition.TSS)[w.species], cex=0.5, adj=c(1,0))
w.species <- which(r_PC[,1]>(0.5)|r_PC[,2]>(0.4))
text(r_PC[w.species,],colnames(Gut_Microbiota_Composition.TSS)[w.species], cex=0.5, adj=c(0,1))

### DENDROGRAMME AND COMPOSITION PLOT
### Dendrogramme of the Dissimilarity Matrix
Hiearachical_ascendant_clustering <- hclust(vegdist(Gut_Microbiota_Composition.TSS), method = "ward.D")
dt=as.dendrogram(Hiearachical_ascendant_clustering)
labels(dt) <- rep("",50) #### Labels suppression

### Dendrogramme and Composition Plot
layout(matrix(c(1,1,2,3,3,3), nrow = 1))
par(mar=c(5,1,2,1))
dt %>% plot(horiz = T, main='', axes=F)
plot(x=rep(1, 50), y = seq(1,50,by=1),  col=colors.ratio[Hiearachical_ascendant_clustering$order], 
     axes=F, xlab="", cex=1.5,
     pch=PCH[Hiearachical_ascendant_clustering$order])

### Selection of the most important GENUS
w=which(apply(Gut_Microbiota_Composition.TSS.GENUS, 2, function(x) median(x)>0.009)==TRUE) #### Most Expressed Genus
COMPOSITION=t(Gut_Microbiota_Composition.TSS.GENUS)[w,]
Others=1-colSums(COMPOSITION)
COMPOSITION=rbind(COMPOSITION,Others)

### Colors of the GENUS
col_genus <- paste0(RColorBrewer::brewer.pal(12,'Set3'),"90")
col_genus[c(1,4)]=c('#44AA9990','#AA449990')

### Barplot of the gut microbiota compositions
barplot(COMPOSITION[,Hiearachical_ascendant_clustering$order], horiz = T, 
        col =col_genus,names.arg =rep('', 50), main='', border=NA)


                          