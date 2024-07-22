rm(list=ls())
library(forcats)
library(rpart)
#############################
##### DONORS DATA  ##########
#############################
Physiological_data <- read.csv("~/Desktop/Exomic_project/data/Humans/Physiological_data.csv", row.names = 1)
Gut_Microbiota_Composition <- read.csv("~/Desktop/Exomic_project/data/Humans/Gut_Microbiota_Composition.csv", row.names = 1)
donors_ID <- levels(as.factor(Physiological_data$DONOR_ID))[-9]
DONORS = NULL
for(i in donors_ID){
  DONORS = rbind(DONORS,Physiological_data[Physiological_data$DONOR_ID==i,])
}
DONORS_composition <- t(Gut_Microbiota_Composition[rownames(DONORS),])

D_GROUP =NULL
for (i in 1:8){
  D <- donors_ID[1:8][i]
  D_GROUP <- c(D_GROUP,Physiological_data$Groupe[Physiological_data$DONOR_ID==D])
}
Names.donors <- paste0(donors_ID[1:8], " (",D_GROUP, ")")

###############################
##### MICE DATA  ##############
###############################
Physiological_data.MICE <- readxl::read_excel("~/Desktop/Exomic_project/data/Mice_FMT/meta_data_mice.xls")
Gut_Microbiota_Functionality.MICE <- read.csv("~/Desktop/Exomic_project/data/Mice_FMT/Gut_Microbiota_Functionality_mice.csv", row.names=1)
Gut_Microbiota_Composition.MICE <- read.csv("~/Desktop/Exomic_project/data/Mice_FMT/Gut_Microbiota_Composition_mice.csv", row.names=1)
Metabolites_data.MICE <- read.csv("~/Desktop/Exomic_project/data/Mice_FMT/Metabolites_data_mice.csv", row.names=1)
Species_labelled.MICE <- read.csv("~/Desktop/Exomic_project/data/Mice_FMT/Species_Phylo_Label_mice.csv", row.names=1)
filt_tree_MICE <- ape::read.tree("~/Desktop/Exomic_project/data/Mice_FMT/Phylogenetic_Tree_mice.txt")

##########################################################################################
MAP_METACYC <- read.table("~/Dropbox/THESE_DM/EXOMIC 2022/HUMANN3/map_metacyc.csv",  sep=';')
labels=substr(MAP_METACYC$V2,1,12)
MAP_PATH<-tidyr::separate(MAP_METACYC, "V2",c("A","B","C","D"), sep='\\|' )
MAP_PATH <- MAP_PATH[MAP_PATH$A!="Super-Pathways",]
##########################################################################################

##### DATA PREPROCESS.
## FEATURES FOR THE PLOTS=======
CTL.mice <- grep("No_FMT",Physiological_data.MICE$DONORS) # Control mice
COL.mice <- factor(as.factor(Physiological_data.MICE$DONORS), # Colors of mice in plot
                   levels(as.factor(Physiological_data.MICE$DONORS)), c(khroma::color("muted")(9)))
PCH.mice <- factor(as.factor(Physiological_data.MICE$DONORS), # PCH of mice in plot
                   levels(as.factor(Physiological_data.MICE$DONORS)), c(1:8, 15))
Mice.D <- Physiological_data.MICE$DONORS[-CTL.mice]
rownames(DONORS) <- DONORS$DONOR_ID

## Physiological FEATURES=======
# of mice
R_EX.CAP <- Physiological_data.MICE$CMA_4/Physiological_data.MICE$CMA_1
R_GLYCO_MR <- as.numeric(Physiological_data.MICE$GASTRO_GLUCOSE)
GLYCEMIE <- as.numeric(Physiological_data.MICE$GLYCEMIE)
R_FM <- as.numeric(Physiological_data.MICE$FatMass)/as.numeric(Physiological_data.MICE$WEIGTH_final)
R_INSULIN <- as.numeric(Physiological_data.MICE$INSULIN)
R_HOMA = (GLYCEMIE/100*5.5)*(R_INSULIN*25)/22.5
FOOD <- as.numeric(Physiological_data.MICE$FOOD_intakes)
WEIGHT_GAIN <- as.numeric(Physiological_data.MICE$WEIGTH_final)-as.numeric(Physiological_data.MICE$WEIGHT_initial)

# of gut microbiota of mice
RECEVORS.W9.Composition <- data.frame(Gut_Microbiota_Composition.MICE[grep("W9", rownames(Gut_Microbiota_Composition.MICE)),])
RECEVORS.W3.Composition <- data.frame(Gut_Microbiota_Composition.MICE[grep("W3", rownames(Gut_Microbiota_Composition.MICE)),])

R_Shannon_Compo. <- vegan::diversity(RECEVORS.W9.Composition)
R_Richness_Compo. <- apply(RECEVORS.W9.Composition, 1, function(x) sum(x>0))

# of Donors
D_VO2 <- DONORS[Mice.D,]$VO2_max
D_EEE <- DONORS[Mice.D,]$Lipides_SV1 + DONORS[Mice.D,]$Glucides_SV1
D_FO <- DONORS[Mice.D,]$Lipides_SV1/D_EEE
D_FM <- DONORS[Mice.D,]$Fat_Mass

##### DF creation of only transfected mice with their own donors.
df.D_R <- data.frame(
  ### the four features from Mice transfected.
  "R_EC" = R_EX.CAP[-CTL.mice], "R_SM-Gly."= R_GLYCO_MR[-CTL.mice],
  "R_FM"=R_FM[-CTL.mice],'R_HOMA-IR'=R_HOMA[-CTL.mice],
  'FOOD Intakes' = FOOD[-CTL.mice], "Weight Gain" = WEIGHT_GAIN[-CTL.mice],
  ### The diversities index of receveirs
  R_Shannon_Compo.=R_Shannon_Compo.[-CTL.mice], R_Richness_Compo.=R_Richness_Compo.[-CTL.mice],
  ### the four features from Donors.
  "D_VO2max"= D_VO2, D_FO, D_EEE, D_FM) 

####################
### Figure S5A #####
####################
BM <- cbind(as.numeric(Physiological_data.MICE$DENSITY_W2),as.numeric(Physiological_data.MICE$DENSITY_W3),as.numeric(Physiological_data.MICE$DENSITY_W8))
par(mar=c(5,5,5,2))
boxplot(BM, at=c(1,3,4), ylim=c(0, 11),
     axes=F, xlab="", ylab="Bacterial density\n(µg DNA/g of feces)")
axis(1, at=c(1,3,4), labels=F)
text(c(1,3,4),rep(-1.4, 4), c("Pre-FMT", "After ABX", "Post-FMT"), 
     srt = 45, xpd=NA, adj=c(1,1), cex = 1)
axis(2)
segments(x0=c(1,1,3), x1=c(3,4,4), y0=c(13, 14, 12), lwd=2 , xpd=NA)
text(2,13.5, "***", cex=1.5, xpd=NA)
text(2.5,14.5, "**", cex=1.5, xpd=NA)
text(3.5,12.5, "***", cex=1.5, xpd=NA)
abline(v=c(1.5, 2.5), lty="dashed")
text(2,5,"ANTIBIOTICS", srt=90)

####################
### Figure S5B #####
####################
quality.transf = id.mice = BC_postFMT = BC_preFMT = quality.transf_abs = SER = NULL
for (i in 1:8){
  w.recevors <- which(Physiological_data.MICE$DONORS==rev(DONORS$DONOR_ID)[i])
  print(w.recevors)
  receve <- merge(t(RECEVORS.W3.Composition)[,w.recevors], t(RECEVORS.W9.Composition)[,w.recevors], by="row.names", all=T)
  rownames(receve) <- receve$Row.names
  receve <- as.matrix(merge(receve[,-1],DONORS_composition[,i]/sum(DONORS_composition[,i]), by="row.names", all=T))
  species_names <- receve[,1]
  receve <- receve[,-1]
  receve[which(is.na(receve))]=0
  receve = apply(receve,2, as.numeric)
  
  id.mice <- c(id.mice, colnames(t(RECEVORS.W9.Composition)[,w.recevors]))
  BC_postFMT <- c(BC_postFMT,as.matrix(vegan::vegdist(t(receve)))[7,4:6])
  
  for (j in 1:3){
    BC_preFMT <- c(BC_preFMT,as.matrix(vegan::vegdist(t(receve)))[j,j+3])
    # Shared across species
    across.H_M <- rowSums(receve[,c(j, j+3, 7 )]>0)=="3"
    across <- length(which(across.H_M==T))
    across_abs <- sum(receve[which(rowSums(receve[,c(j, j+3, 7 )]>0)=="3"),j+3])
    # Shared in mice
    intra.mice <- length(which(rowSums(receve[across.H_M==F, c(j,j+3)]>0)=="2"))
    intra.mice_abs <- sum(receve[across.H_M==F,][which(rowSums(receve[across.H_M==F,c(j,j+3)]>0)=="2"),j+3])
    # Due to FMT
    FMT <- length(which(rowSums(receve[across.H_M==F,c(j+3, 7 )]>0)=="2"))
    FMT_abs <- sum(receve[across.H_M==F,][which(rowSums(receve[across.H_M==F,c(j+3,7)]>0)=="2"),j+3])
    # REC
    quality.transf <- cbind(quality.transf,c(across,  FMT,intra.mice, length(which(receve[,c(j+3 )]>0))))
    quality.transf_abs <- cbind(quality.transf_abs,c(across_abs, FMT_abs, intra.mice_abs, 1-sum(across_abs,intra.mice_abs,FMT_abs)))
    
    # Species engraftment rate
    SER <- c(SER,(FMT+across)/sum(receve[,7]>0))
  }
}
colSums(quality.transf_abs)
colnames(quality.transf)=colnames(quality.transf_abs)=id.mice
rownames(quality.transf)=rownames(quality.transf_abs)=c("Across Species",  "FMT", "Intra species","Reste")
quality.transf[4,] <- quality.transf[4,]-apply(quality.transf[1:3,], 2, sum)
quality.transf <- apply(quality.transf, 2, function(x) x/sum(x))

par(mar=c(5,5,2,2))
plot((quality.transf[1,]+quality.transf[3,]),(quality.transf[1,]+quality.transf[2,]), xlim=c(0,1), ylim=c(0,1),
     xlab="Fraction of post-FMT species\n shared with pre-FMT", 
     ylab="Fraction of post-FMT species\n shared with donors", pch=as.numeric(PCH.mice[-CTL.mice]))
abline(0,1, lty='dashed')

par(mar=c(5,1,2,3.5))
plot(x=rep(1, 8), y = seq(1,8,by=1), axes=F, xlab="", 
     pch=1:8, cex=1)
text(x=rep(0.8, 8), y = seq(1,8,by=1), Names.donors)

layout(mat = matrix(c(1:24), ncol = 3, byrow = T))
par(mar=c(rep(0.5,4)))
col.barplot=RColorBrewer::brewer.pal(8,'Set2')
for(i in 1:24){
  barplot(as.matrix(quality.transf[,i]),axes=F, horiz =T, col=col.barplot[5:8])
}

layout(matrix(c(1,2), nrow=1))
par(mar=c(5,5,2,2))
boxplot(mean(SER), ylim=c(0,0.5), ylab='Species engraftment rate')
points(rep(1, 24), y=SER, pch=as.numeric(PCH.mice[-CTL.mice]))
boxplot(mean(1-BC_postFMT),mean(1-BC_preFMT), ylim=c(0,0.5), ylab='Bray-Curtis similarity')
points(rep(c(1,2), each=24), y=c(1-BC_postFMT,1-BC_preFMT), pch=as.numeric(PCH.mice[-CTL.mice]))
segments(x0=1,x1=2, y0=.45, lwd=2)
text(1.5, .48, "**", cex=1.5)
text(c(1,2), -.05, xpd=NA, labels = c("PostFMT/Donors", "PreFMT/PostFMT"), 
     srt=45, adj=1, cex=0.8)

###############################
###### Figure S5C-D ###########
###############################
dev.off()
layout(matrix(c(1:2), nrow=1))
### BC-W3 =============
ind.coord_w3 <- ape::pcoa(vegan::vegdist(RECEVORS.W3.Composition))$vectors
plot(ind.coord_w3, pch=as.numeric(as.character(PCH.mice)),col="black", 
     xlab="PCoA 1 (27%)", ylab="PCoA 2 (16%)", main = "Before FMT",cex=1)
### BC-W9 =============
ind.coord_w9 <- ape::pcoa(vegan::vegdist(RECEVORS.W9.Composition))$vectors
plot(ind.coord_w9, pch=as.numeric(as.character(PCH.mice)),col="black",
     xlab="PCoA 1 (20%)", ylab="PCoA 2 (12%)", cex=1, main = "After FMT")
dev.off()

###### LEGENDS ==================================
par(mar=c(0,0,0,0))
plot.new()
legend(0,0.8, legend=c("CTL", Names.donors), ncol = 3,pch=c(15,1:8),
       bty='n',title.adj=0, cex=1, title="Mice experiments")
legend(0,0.4, legend=c("PreFMT/PostFMT/Donors",  "PostFMT/Donors","PreFMT/PostFMT", "Others"), 
       ncol = 2,fill=col.barplot[5:8], bty='n',title = "Shared Species",title.adj=0, cex=1)
