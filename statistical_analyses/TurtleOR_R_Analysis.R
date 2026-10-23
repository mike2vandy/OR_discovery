library(phytools)
library(ggplot2)
library(phylolm)
library(ggpubr)
library(hablar)
library(ape)
library(plyr)
library(cowplot)

#global variables
turtledata <- read.csv("Turtlegroupings.csv", header = TRUE, row.names = 1) #groupings
turtles <- read.tree("Timescaled.Turtles.nwk")

#plot turtle tree
plot(turtles, edge.width = 3)
add.scale.bar(length=50, lwd = 3)

###CLASS I PROPORTION BY AQUATIC AFFINITY###
#PGLS
ClassI_by_aquatic <-data.frame(X = turtledata$Habitat, Y = turtledata$Class.I.Prop, row.names = rownames(turtledata))
ClassI_by_aquatic_PGLS <- phylolm(Y~X, data = ClassI_by_aquatic, phy=turtles, model = "BM")
summary(ClassI_by_aquatic_PGLS)
#plot and save results
ClassI.plot <- ggplot(ClassI_by_aquatic, aes(x= X, y = Y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "Aquatic affinity",
       y = "Class I proportion")
ggsave("ClassI_by_aquatic.pdf", ClassI.plot, 
       device = 'pdf', units = 'in', height = 5, width = 6.5)

#####################################

###PSEUDOGENE PROPORTION BY INTACT REPERTOIRE SIZE###
#PGLS
Pseudo_prop_by_intact <- data.frame(X = turtledata$Total.Intact, Y= turtledata$Pseudogene.proportion, row.names= rownames(turtledata))
Pseudo_prop_by_intact_PGLS <- phylolm(Y~X, data = Pseudo_prop_by_intact, phy= turtles, model = "BM")
summary(Pseudo_prop_by_intact_PGLS)
#plot and save results
Pseudo.plot <- ggplot(Pseudo_prop_by_intact, aes(x= X, y = Y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "Number of intact genes",
       y = "Proportion of Pseudogenes (pseudo/pseudo+intact)")
ggsave("Pseudo_prop_by_intact.pdf", Pseudo.plot, 
       device = 'pdf', units = 'in', height = 5, width = 5.5)

###################################

###CREATE FUNCTION FOR PCAS AND PLOTS###

#Creates the phylo PCA function for all of the PCAs
#Note: to change grouping factor for visualization, adjust "Primary.Habitat" to another discrete variable in all 3 lines where it appears
phyloPCA <- function(newick, norm, meta, loadScale) {
  phyPCA <- phyl.pca(newick, norm, method = 'BM', mode = 'cov')
  
  variance <- diag(phyPCA$Eval)/sum(phyPCA$Eval)*100
  xLab <- as.character(paste("PC1 (", round(variance[1],2),"% variance)", sep = "")) 
  yLab <- as.character(paste("PC2 (", round(variance[2],2), "% variance)", sep = ""))
  
  pcaPl <- as.data.frame(phyPCA$S)
  pcAndMeta <- merge(pcaPl, meta, by = 0)
  
  #geom_polygon things
  find_hulls <- function(pcAndMeta) pcAndMeta[chull(pcAndMeta$PC1, pcAndMeta$PC2),]
  hulls <- ddply(pcAndMeta, "Primary.Habitat", find_hulls)
  
  #loadings
  loadings <- as.data.frame(phyPCA$L)
  
  #main ggplot
  PCAplot <- ggplot() + 
    geom_segment(data = loadings, 
                 aes(x = 0, y = 0, xend = PC1*loadScale, yend = PC2*loadScale), 
                 arrow = arrow(length = unit(1/2, "picas")),
                 color = "#CC2B18", size = 1) +
    #terrestrial #semi-aquatic #aquatic # marine
    scale_color_manual(values = c("#C46E42","#80BB83", "#467076", "#9A3856")) +
    scale_fill_manual(values = c("#C46E42","#80BB83", "#467076", "#9A3856")) +
    geom_point(data = pcAndMeta, 
               aes(x = PC1, y = PC2, color = Primary.Habitat), size = 3) +
    geom_polygon(data = hulls, alpha = 0.3, 
                 aes(x = PC1, y = PC2, fill = Primary.Habitat)) +
    annotate("text", x = (loadings$PC1*loadScale), y = (loadings$PC2*loadScale),
             label = row.names(loadings), color = "#CC2B18") +
    geom_text(data = pcAndMeta, 
              aes(x = PC1, y = PC2, label = Species), 
              vjust = 1, hjust = -0.1) +
    xlab(xLab) +
    ylab(yLab) +
    #xlim(-0.25, 0.2) +
    #ylim(-0.1, 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(size = 14))
  
  return(PCAplot)
}

###############################

###INTACT OR COUNT ANALYSIS###

#read in the gene counts divided by subfamily, transform to matrix
turtleor.count <- read.csv("TurtleOR.count.csv", header = TRUE, row.names = 1) #raw data
turtleorcount.matrix <- as.matrix(turtleor.count)

#PCA
turtleorcount.pca <- phyl.pca(turtles, turtleorcount.matrix, method = 'BM', mode = 'cov')

#visualize and save
intactGraph <- phyloPCA(turtles, turtleorcount.matrix, turtledata, 200)
plot(intactGraph)
ggsave("intact_count.pdf", intactGraph, 
       device = 'pdf', units = 'in', height = 6, width = 11)

#Store scores
turtleorcount.scores <- turtleorcount.pca$S

#set up the DF for the PGLS... move over scores, and set species names and groups
turtleorcount.df <- as.data.frame(turtleorcount.scores)
turtleorcount.df$species <- rownames(turtleor.count)
turtleorcount.df$primary_habitat <- turtledata$Primary.Habitat
turtleorcount.df$diet <- turtledata$Diet
turtleorcount.df$aquatic_cont <- turtledata$Habitat

#save data frame as csv
write.csv(turtleorcount.df, file = "Raw_OR_Count_PCA_groups.csv", row.names = FALSE)

#PGLS

#AQUATIC AFFINITY
#PC1
PC1_by_aquatic <- data.frame(X = turtleorcount.df$aquatic_cont, Y = turtleorcount.df$PC1, row.names = turtleorcount.df$species)
PC1_by_aquatic_PGLS <- phylolm(Y~X, data = PC1_by_aquatic, phy=turtles, model = "BM")
summary(PC1_by_aquatic_PGLS)
#PC2
PC2_by_aquatic <- data.frame(X = turtleorcount.df$aquatic_cont, Y = turtleorcount.df$PC2, row.names = turtleorcount.df$species)
PC2_by_aquatic_PGLS <- phylolm(Y~X, data = PC2_by_aquatic, phy=turtles, model = "BM")
summary(PC2_by_aquatic_PGLS)

#GENERAL HABITAT
#PC1
PC1_by_hab <- data.frame(X = turtleorcount.df$primary_habitat, Y = turtleorcount.df$PC1, row.names = turtleorcount.df$species)
PC1_by_hab_PGLS <- phylolm(Y~X, data = PC1_by_hab, phy=turtles, model = "BM")
summary(PC1_by_hab_PGLS)
#PC2
PC2_by_hab <- data.frame(X = turtleorcount.df$primary_habitat, Y = turtleorcount.df$PC2, row.names = turtleorcount.df$species)
PC2_by_hab_PGLS <- phylolm(Y~X, data = PC2_by_hab, phy=turtles, model = "BM")
summary(PC2_by_hab_PGLS)

#DIET
#PC1
PC1_by_diet <- data.frame(X = turtleorcount.df$diet, Y = turtleorcount.df$PC1, row.names = turtleorcount.df$species)
PC1_by_diet_PGLS <- phylolm(Y~X, data = PC1_by_diet, phy=turtles, model = "BM")
summary(PC1_by_diet_PGLS)
#PC2
PC2_by_diet <- data.frame(X = turtleorcount.df$diet, Y = turtleorcount.df$PC2, row.names = turtleorcount.df$species)
PC2_by_diet_PGLS <- phylolm(Y~X, data = PC2_by_diet, phy=turtles, model = "BM")
summary(PC2_by_diet_PGLS)


##################################################

###INTACT OR PROPORTION ANALYSIS###

#Read in the intact OR proportions (defined as intact subfamily count/total intact gene count), transform to matrix
turtleor.proportions <- read.csv("TurtleOR.proportions.csv", header = TRUE, row.names = 1) #raw data
turtleorproportions.matrix <- data.matrix(turtleor.proportions)

#PCA
turtleorproportions.pca <- phyl.pca(turtles, turtleorproportions.matrix, method = 'BM', mode = 'cov')

#visualize and save
intactpropGraph <- phyloPCA(turtles, turtleorproportions.matrix, turtledata, 0.1)
plot(intactpropGraph)
ggsave("intact_prop_PCA.pdf", intactpropGraph, 
       device = 'pdf', units = 'in', height = 6, width = 11)

#set up DF for PGLS
turtleorprop.scores <- turtleorproportions.pca$S
turtleorprop.df <- as.data.frame(turtleorprop.scores)
turtleorprop.df$species <- rownames(turtleor.proportions)
turtleorprop.df$primary_habitat <- turtledata$Primary.Habitat
turtleorprop.df$diet <- turtledata$Diet
turtleorprop.df$aquatic_cont <- turtledata$Habitat

#save data frame as csv
write.csv(turtleorprop.df, file = "OR_proportions_PCA_groups.csv", row.names = FALSE)

#PGLS

#AQUATIC AFFINITY
#PC1
PC1_by_aquatic.prop <- data.frame(X = turtleorprop.df$aquatic_cont, Y = turtleorprop.df$PC1, row.names = turtleorprop.df$species)
PC1_by_aquatic_PGLS.prop <- phylolm(Y~X, data = PC1_by_aquatic.prop, phy=turtles, model = "BM")
summary(PC1_by_aquatic_PGLS.prop)
#PC2
PC2_by_aquatic.prop <- data.frame(X = turtleorprop.df$aquatic_cont, Y = turtleorprop.df$PC2, row.names = turtleorprop.df$species)
PC2_by_aquatic_PGLS.prop <- phylolm(Y~X, data = PC2_by_aquatic.prop, phy=turtles, model = "BM")
summary(PC2_by_aquatic_PGLS.prop)

#HABITAT
#PC1
PC1_by_hab.prop <- data.frame(X = turtleorprop.df$primary_habitat, Y = turtleorprop.df$PC1, row.names = turtleorprop.df$species)
PC1_by_hab_PGLS.prop <- phylolm(Y~X, data = PC1_by_hab.prop, phy=turtles, model = "BM")
summary(PC1_by_hab_PGLS.prop)
#PC2
PC2_by_hab.prop <- data.frame(X = turtleorprop.df$primary_habitat, Y = turtleorprop.df$PC2, row.names = turtleorprop.df$species)
PC2_by_hab_PGLS.prop <- phylolm(Y~X, data = PC2_by_hab.prop, phy=turtles, model = "BM")
summary(PC2_by_hab_PGLS.prop)

#DIET
#PC1
PC1_by_diet.prop <- data.frame(X = turtleorprop.df$diet, Y = turtleorprop.df$PC1, row.names = turtleorprop.df$species)
PC1_by_diet_PGLS.prop <- phylolm(Y~X, data = PC1_by_diet.prop, phy=turtles, model = "BM")
summary(PC1_by_diet_PGLS.prop)
#PC2
PC2_by_diet.prop <- data.frame(X = turtleorprop.df$diet, Y = turtleorprop.df$PC2, row.names = turtleorprop.df$species)
PC2_by_diet_PGLS.prop <- phylolm(Y~X, data = PC2_by_diet.prop, phy=turtles, model = "BM")
summary(PC2_by_diet_PGLS.prop)

#################################

###PSEUDOGENE PROPORTION ANALYSIS###

#Read in the pseudogene proportions (number of pseudogenes/number of intact + pseudogenes within each subfamily), transform to matrix
turtleor.pseudo.proportions <- read.csv("TurtleOR.pseudo.prop.no55.no1.csv", header = TRUE, row.names = 1) #raw data
turtleor.pseudo.prop.matrix <- data.matrix(turtleor.pseudo.proportions)

#PCA
turtleor.pseudo.prop.pca <- phyl.pca(turtles, turtleor.pseudo.prop.matrix, method = 'BM', mode = 'cov')

#Visualize and save
PseudoGraph <- phyloPCA(turtles, turtleor.pseudo.prop.matrix, turtledata, 1)
plot(PseudoGraph)
ggsave("pseudo.prop.no55.no1.diet.pdf", PseudoGraph, 
       device = 'pdf', units = 'in', height = 6, width = 11)

#set up DF
turtleor.pseudo.prop.scores <- turtleor.pseudo.prop.pca$S
turtleor.pseudo.prop.df <- as.data.frame(turtleor.pseudo.prop.scores)
turtleor.pseudo.prop.df$species <- rownames(turtleor.pseudo.proportions)
turtleor.pseudo.prop.df$primary_habitat <- turtledata$Primary.Habitat
turtleor.pseudo.prop.df$diet <- turtledata$Diet
turtleor.pseudo.prop.df$aquatic_cont <- turtledata$Habitat

#save data frame as csv
write.csv(turtleor.pseudo.prop.df, file = "OR_pseudoprop_no55_no1_PCA_groups.csv", row.names = FALSE)

#PGLS

#AQUATIC AFFINITY
#PC1
PC1_by_aquatic.pseudo <- data.frame(X = turtleor.pseudo.prop.df$aquatic_cont, Y = turtleor.pseudo.prop.df$PC1, row.names = turtleor.pseudo.prop.df$species)
PC1_by_aquatic_PGLS.pseudo <- phylolm(Y~X, data = PC1_by_aquatic.pseudo, phy=turtles, model = "BM")
summary(PC1_by_aquatic_PGLS.pseudo)
#PC2
PC2_by_aquatic.pseudo <- data.frame(X = turtleor.pseudo.prop.df$aquatic_cont, Y = turtleor.pseudo.prop.df$PC2, row.names = turtleor.pseudo.prop.df$species)
PC2_by_aquatic_PGLS.pseudo <- phylolm(Y~X, data = PC2_by_aquatic.pseudo, phy=turtles, model = "BM")
summary(PC2_by_aquatic_PGLS.pseudo)

#HABITAT
#PC1
PC1_by_hab.pseudo <- data.frame(X = turtleor.pseudo.prop.df$primary_habitat, Y = turtleor.pseudo.prop.df$PC1, row.names = turtleor.pseudo.prop.df$species)
PC1_by_hab_PGLS.pseudo <- phylolm(Y~X, data = PC1_by_hab.pseudo, phy=turtles, model = "BM")
summary(PC1_by_hab_PGLS.pseudo)
#PC2
PC2_by_hab.pseudo <- data.frame(X = turtleor.pseudo.prop.df$primary_habitat, Y = turtleor.pseudo.prop.df$PC2, row.names = turtleor.pseudo.prop.df$species)
PC2_by_hab_PGLS.pseudo <- phylolm(Y~X, data = PC2_by_hab.pseudo, phy=turtles, model = "BM")
summary(PC2_by_hab_PGLS.pseudo)

#DIET
#PC1
PC1_by_diet.pseudo <- data.frame(X = turtleor.pseudo.prop.df$diet, Y = turtleor.pseudo.prop.df$PC1, row.names = turtleor.pseudo.prop.df$species)
PC1_by_diet_PGLS.pseudo <- phylolm(Y~X, data = PC1_by_diet.pseudo, phy=turtles, model = "BM")
summary(PC1_by_diet_PGLS.pseudo)
#PC2
PC2_by_diet.pseudo <- data.frame(X = turtleor.pseudo.prop.df$diet, Y = turtleor.pseudo.prop.df$PC2, row.names = turtleor.pseudo.prop.df$species)
PC2_by_diet_PGLS.pseudo <- phylolm(Y~X, data = PC2_by_diet.pseudo, phy=turtles, model = "BM")
summary(PC2_by_diet_PGLS.pseudo)

#####################################

###INDIVIDUAL SUBFAMILY ANALYSIS###

#create df to store data
subfamily.analysis.df <- data.frame(matrix(ncol = 12, nrow = 13))
row.names(subfamily.analysis.df) = c('OR1', 'OR2', 'OR4', 'OR5.9', 'OR6', 'OR8', 'OR10', 'OR11', 'OR12', 'OR13', 'OR14', 'OR51', 'OR52')
colnames(subfamily.analysis.df) <- c('Count.Estimate', 'Count.STDERR', 'Count.Tval', 'Count.raw.pval',  'Count.BF.pval', 'Count.r.squared','Prop.Estimate', 'Prop.STDERR', 'Prop.Tval', 'Prop.raw.pval', 'Prop.BF.pval', 'Prop.r.squared')

#INTACT COUNTS

#create input dataframe
turtleorcount.subfam.df <- as.data.frame(turtleorcount.matrix) #turtleorcount.matrix created in INTACT OR COUNT ANALYSIS
turtleorcount.subfam.df$species <- rownames(turtleor.count)
turtleorcount.subfam.df$aquatic_cont <- turtledata$Habitat

#PGLS

#OR1
OR1_by_aquatic <- data.frame(X = turtleorcount.subfam.df$aquatic_cont, Y = turtleorcount.subfam.df$OR1, row.names = turtleorcount.subfam.df$species)
OR1_by_aquatic_PGLS <- phylolm(Y~X, data = OR1_by_aquatic, phy=turtles, model = "BM")
OR1<-summary(OR1_by_aquatic_PGLS)
subfamily.analysis.df[1,1] <-OR1$coefficients[2,1]
subfamily.analysis.df[1,2] <-OR1$coefficients[2,2]
subfamily.analysis.df[1,3] <-OR1$coefficients[2,3]
subfamily.analysis.df[1,4] <-OR1$coefficients[2,4]
subfamily.analysis.df[1,6] <-OR1$r.squared
#OR2
OR2_by_aquatic <- data.frame(X = turtleorcount.subfam.df$aquatic_cont, Y = turtleorcount.subfam.df$OR2, row.names = turtleorcount.subfam.df$species)
OR2_by_aquatic_PGLS <- phylolm(Y~X, data = OR2_by_aquatic, phy=turtles, model = "BM")
OR2<-summary(OR2_by_aquatic_PGLS)
subfamily.analysis.df[2,1] <-OR2$coefficients[2,1]
subfamily.analysis.df[2,2] <-OR2$coefficients[2,2]
subfamily.analysis.df[2,3] <-OR2$coefficients[2,3]
subfamily.analysis.df[2,4] <-OR2$coefficients[2,4]
subfamily.analysis.df[2,6] <-OR2$r.squared
#OR4
OR4_by_aquatic <- data.frame(X = turtleorcount.subfam.df$aquatic_cont, Y = turtleorcount.subfam.df$OR4, row.names = turtleorcount.subfam.df$species)
OR4_by_aquatic_PGLS <- phylolm(Y~X, data = OR4_by_aquatic, phy=turtles, model = "BM")
OR4<-summary(OR4_by_aquatic_PGLS)
subfamily.analysis.df[3,1] <-OR4$coefficients[2,1]
subfamily.analysis.df[3,2] <-OR4$coefficients[2,2]
subfamily.analysis.df[3,3] <-OR4$coefficients[2,3]
subfamily.analysis.df[3,4] <-OR4$coefficients[2,4]
subfamily.analysis.df[3,6] <-OR4$r.squared
#OR5/9
OR5.9_by_aquatic <- data.frame(X = turtleorcount.subfam.df$aquatic_cont, Y = turtleorcount.subfam.df$OR5_9, row.names = turtleorcount.subfam.df$species)
OR5.9_by_aquatic_PGLS <- phylolm(Y~X, data = OR5.9_by_aquatic, phy=turtles, model = "BM")
OR5.9<-summary(OR5.9_by_aquatic_PGLS)
subfamily.analysis.df[4,1] <-OR5.9$coefficients[2,1]
subfamily.analysis.df[4,2] <-OR5.9$coefficients[2,2]
subfamily.analysis.df[4,3] <-OR5.9$coefficients[2,3]
subfamily.analysis.df[4,4] <-OR5.9$coefficients[2,4]
subfamily.analysis.df[4,6] <-OR5.9$r.squared
#OR6
OR6_by_aquatic <- data.frame(X = turtleorcount.subfam.df$aquatic_cont, Y = turtleorcount.subfam.df$OR6, row.names = turtleorcount.subfam.df$species)
OR6_by_aquatic_PGLS <- phylolm(Y~X, data = OR6_by_aquatic, phy=turtles, model = "BM")
OR6<-summary(OR6_by_aquatic_PGLS)
subfamily.analysis.df[5,1] <-OR6$coefficients[2,1]
subfamily.analysis.df[5,2] <-OR6$coefficients[2,2]
subfamily.analysis.df[5,3] <-OR6$coefficients[2,3]
subfamily.analysis.df[5,4] <-OR6$coefficients[2,4]
subfamily.analysis.df[5,6] <-OR6$r.squared
#OR8
OR8_by_aquatic <- data.frame(X = turtleorcount.subfam.df$aquatic_cont, Y = turtleorcount.subfam.df$OR8, row.names = turtleorcount.subfam.df$species)
OR8_by_aquatic_PGLS <- phylolm(Y~X, data = OR8_by_aquatic, phy=turtles, model = "BM")
OR8<-summary(OR8_by_aquatic_PGLS)
subfamily.analysis.df[6,1] <-OR8$coefficients[2,1]
subfamily.analysis.df[6,2] <-OR8$coefficients[2,2]
subfamily.analysis.df[6,3] <-OR8$coefficients[2,3]
subfamily.analysis.df[6,4] <-OR8$coefficients[2,4]
subfamily.analysis.df[6,6] <-OR8$r.squared
#OR10
OR10_by_aquatic <- data.frame(X = turtleorcount.subfam.df$aquatic_cont, Y = turtleorcount.subfam.df$OR10, row.names = turtleorcount.subfam.df$species)
OR10_by_aquatic_PGLS <- phylolm(Y~X, data = OR10_by_aquatic, phy=turtles, model = "BM")
summary(OR10_by_aquatic_PGLS)
OR10<-summary(OR10_by_aquatic_PGLS)
subfamily.analysis.df[7,1] <-OR10$coefficients[2,1]
subfamily.analysis.df[7,2] <-OR10$coefficients[2,2]
subfamily.analysis.df[7,3] <-OR10$coefficients[2,3]
subfamily.analysis.df[7,4] <-OR10$coefficients[2,4]
subfamily.analysis.df[7,6] <-OR10$r.squared
#OR11
OR11_by_aquatic <- data.frame(X = turtleorcount.subfam.df$aquatic_cont, Y = turtleorcount.subfam.df$OR11, row.names = turtleorcount.subfam.df$species)
OR11_by_aquatic_PGLS <- phylolm(Y~X, data = OR11_by_aquatic, phy=turtles, model = "BM")
OR11<-summary(OR11_by_aquatic_PGLS)
subfamily.analysis.df[8,1] <-OR11$coefficients[2,1]
subfamily.analysis.df[8,2] <-OR11$coefficients[2,2]
subfamily.analysis.df[8,3] <-OR11$coefficients[2,3]
subfamily.analysis.df[8,4] <-OR11$coefficients[2,4]
subfamily.analysis.df[8,6] <-OR11$r.squared
#OR12
OR12_by_aquatic <- data.frame(X = turtleorcount.subfam.df$aquatic_cont, Y = turtleorcount.subfam.df$OR12, row.names = turtleorcount.subfam.df$species)
OR12_by_aquatic_PGLS <- phylolm(Y~X, data = OR12_by_aquatic, phy=turtles, model = "BM")
OR12<-summary(OR12_by_aquatic_PGLS)
subfamily.analysis.df[9,1] <-OR12$coefficients[2,1]
subfamily.analysis.df[9,2] <-OR12$coefficients[2,2]
subfamily.analysis.df[9,3] <-OR12$coefficients[2,3]
subfamily.analysis.df[9,4] <-OR12$coefficients[2,4]
subfamily.analysis.df[9,6] <-OR12$r.squared
#OR13
OR13_by_aquatic <- data.frame(X = turtleorcount.subfam.df$aquatic_cont, Y = turtleorcount.subfam.df$OR13, row.names = turtleorcount.subfam.df$species)
OR13_by_aquatic_PGLS <- phylolm(Y~X, data = OR13_by_aquatic, phy=turtles, model = "BM")
OR13<-summary(OR13_by_aquatic_PGLS)
subfamily.analysis.df[10,1] <-OR13$coefficients[2,1]
subfamily.analysis.df[10,2] <-OR13$coefficients[2,2]
subfamily.analysis.df[10,3] <-OR13$coefficients[2,3]
subfamily.analysis.df[10,4] <-OR13$coefficients[2,4]
subfamily.analysis.df[10,6] <-OR13$r.squared
#OR14
OR14_by_aquatic <- data.frame(X = turtleorcount.subfam.df$aquatic_cont, Y = turtleorcount.subfam.df$OR14, row.names = turtleorcount.subfam.df$species)
OR14_by_aquatic_PGLS <- phylolm(Y~X, data = OR14_by_aquatic, phy=turtles, model = "BM")
OR14<-summary(OR14_by_aquatic_PGLS)
subfamily.analysis.df[11,1] <-OR14$coefficients[2,1]
subfamily.analysis.df[11,2] <-OR14$coefficients[2,2]
subfamily.analysis.df[11,3] <-OR14$coefficients[2,3]
subfamily.analysis.df[11,4] <-OR14$coefficients[2,4]
subfamily.analysis.df[11,6] <-OR14$r.squared
#OR51
OR51_by_aquatic <- data.frame(X = turtleorcount.subfam.df$aquatic_cont, Y = turtleorcount.subfam.df$OR51, row.names = turtleorcount.subfam.df$species)
OR51_by_aquatic_PGLS <- phylolm(Y~X, data = OR51_by_aquatic, phy=turtles, model = "BM")
OR51<-summary(OR51_by_aquatic_PGLS)
subfamily.analysis.df[12,1] <-OR51$coefficients[2,1]
subfamily.analysis.df[12,2] <-OR51$coefficients[2,2]
subfamily.analysis.df[12,3] <-OR51$coefficients[2,3]
subfamily.analysis.df[12,4] <-OR51$coefficients[2,4]
subfamily.analysis.df[12,6] <-OR51$r.squared
#OR52
OR52_by_aquatic <- data.frame(X = turtleorcount.subfam.df$aquatic_cont, Y = turtleorcount.subfam.df$OR52, row.names = turtleorcount.subfam.df$species)
OR52_by_aquatic_PGLS <- phylolm(Y~X, data = OR52_by_aquatic, phy=turtles, model = "BM")
OR52<-summary(OR52_by_aquatic_PGLS)
subfamily.analysis.df[13,1] <-OR52$coefficients[2,1]
subfamily.analysis.df[13,2] <-OR52$coefficients[2,2]
subfamily.analysis.df[13,3] <-OR52$coefficients[2,3]
subfamily.analysis.df[13,4] <-OR52$coefficients[2,4]
subfamily.analysis.df[13,6] <-OR52$r.squared
#OR55 not analyzed due to gene count

#INTACT PROPORTIONS

#create input DF
turtleorprop.subfam.df <- as.data.frame(turtleorproportions.matrix) #turtleorproportions.matrix created in INTACT OR PROPORTION ANALYSIS
turtleorprop.subfam.df$species <- rownames(turtleor.proportions)
turtleorprop.subfam.df$aquatic_cont <- turtledata$Habitat

#PGLS

#OR1prop
OR1prop_by_aquatic <- data.frame(X = turtleorprop.subfam.df$aquatic_cont, Y = turtleorprop.subfam.df$OR1, row.names = turtleorprop.subfam.df$species)
OR1prop_by_aquatic_PGLS <- phylolm(Y~X, data = OR1prop_by_aquatic, phy=turtles, model = "BM")
OR1prop<-summary(OR1prop_by_aquatic_PGLS)
subfamily.analysis.df[1,7] <-OR1prop$coefficients[2,1]
subfamily.analysis.df[1,8] <-OR1prop$coefficients[2,2]
subfamily.analysis.df[1,9] <-OR1prop$coefficients[2,3]
subfamily.analysis.df[1,10] <-OR1prop$coefficients[2,4]
subfamily.analysis.df[1,12] <-OR1prop$r.squared
#OR2prop
OR2prop_by_aquatic <- data.frame(X = turtleorprop.subfam.df$aquatic_cont, Y = turtleorprop.subfam.df$OR2, row.names = turtleorprop.subfam.df$species)
OR2prop_by_aquatic_PGLS <- phylolm(Y~X, data = OR2prop_by_aquatic, phy=turtles, model = "BM")
OR2prop<-summary(OR2prop_by_aquatic_PGLS)
subfamily.analysis.df[2,7] <-OR2prop$coefficients[2,1]
subfamily.analysis.df[2,8] <-OR2prop$coefficients[2,2]
subfamily.analysis.df[2,9] <-OR2prop$coefficients[2,3]
subfamily.analysis.df[2,10] <-OR2prop$coefficients[2,4]
subfamily.analysis.df[2,12] <-OR2prop$r.squared
#OR4prop
OR4prop_by_aquatic <- data.frame(X = turtleorprop.subfam.df$aquatic_cont, Y = turtleorprop.subfam.df$OR4, row.names = turtleorprop.subfam.df$species)
OR4prop_by_aquatic_PGLS <- phylolm(Y~X, data = OR4prop_by_aquatic, phy=turtles, model = "BM")
OR4prop<-summary(OR4prop_by_aquatic_PGLS)
subfamily.analysis.df[3,7] <-OR4prop$coefficients[2,1]
subfamily.analysis.df[3,8] <-OR4prop$coefficients[2,2]
subfamily.analysis.df[3,9] <-OR4prop$coefficients[2,3]
subfamily.analysis.df[3,10] <-OR4prop$coefficients[2,4]
subfamily.analysis.df[3,12] <-OR4prop$r.squared
#OR5/9
OR5.9prop_by_aquatic <- data.frame(X = turtleorprop.subfam.df$aquatic_cont, Y = turtleorprop.subfam.df$OR5_9, row.names = turtleorprop.subfam.df$species)
OR5.9prop_by_aquatic_PGLS <- phylolm(Y~X, data = OR5.9prop_by_aquatic, phy=turtles, model = "BM")
OR5.9prop<-summary(OR5.9prop_by_aquatic_PGLS)
subfamily.analysis.df[4,7] <-OR5.9prop$coefficients[2,1]
subfamily.analysis.df[4,8] <-OR5.9prop$coefficients[2,2]
subfamily.analysis.df[4,9] <-OR5.9prop$coefficients[2,3]
subfamily.analysis.df[4,10] <-OR5.9prop$coefficients[2,4]
subfamily.analysis.df[4,12] <-OR5.9prop$r.squared
#OR6prop
OR6prop_by_aquatic <- data.frame(X = turtleorprop.subfam.df$aquatic_cont, Y = turtleorprop.subfam.df$OR6, row.names = turtleorprop.subfam.df$species)
OR6prop_by_aquatic_PGLS <- phylolm(Y~X, data = OR6prop_by_aquatic, phy=turtles, model = "BM")
OR6prop<-summary(OR6prop_by_aquatic_PGLS)
subfamily.analysis.df[5,7] <-OR6prop$coefficients[2,1]
subfamily.analysis.df[5,8] <-OR6prop$coefficients[2,2]
subfamily.analysis.df[5,9] <-OR6prop$coefficients[2,3]
subfamily.analysis.df[5,10] <-OR6prop$coefficients[2,4]
subfamily.analysis.df[5,12] <-OR6prop$r.squared
#OR8prop
OR8prop_by_aquatic <- data.frame(X = turtleorprop.subfam.df$aquatic_cont, Y = turtleorprop.subfam.df$OR8, row.names = turtleorprop.subfam.df$species)
OR8prop_by_aquatic_PGLS <- phylolm(Y~X, data = OR8prop_by_aquatic, phy=turtles, model = "BM")
OR8prop<-summary(OR8prop_by_aquatic_PGLS)
subfamily.analysis.df[6,7] <-OR8prop$coefficients[2,1]
subfamily.analysis.df[6,8] <-OR8prop$coefficients[2,2]
subfamily.analysis.df[6,9] <-OR8prop$coefficients[2,3]
subfamily.analysis.df[6,10] <-OR8prop$coefficients[2,4]
subfamily.analysis.df[6,12] <-OR8prop$r.squared
#OR10prop
OR10prop_by_aquatic <- data.frame(X = turtleorprop.subfam.df$aquatic_cont, Y = turtleorprop.subfam.df$OR10, row.names = turtleorprop.subfam.df$species)
OR10prop_by_aquatic_PGLS <- phylolm(Y~X, data = OR10prop_by_aquatic, phy=turtles, model = "BM")
OR10prop<-summary(OR10prop_by_aquatic_PGLS)
subfamily.analysis.df[7,7] <-OR10prop$coefficients[2,1]
subfamily.analysis.df[7,8] <-OR10prop$coefficients[2,2]
subfamily.analysis.df[7,9] <-OR10prop$coefficients[2,3]
subfamily.analysis.df[7,10] <-OR10prop$coefficients[2,4]
subfamily.analysis.df[7,12] <-OR10prop$r.squared
#OR11prop
OR11prop_by_aquatic <- data.frame(X = turtleorprop.subfam.df$aquatic_cont, Y = turtleorprop.subfam.df$OR11, row.names = turtleorprop.subfam.df$species)
OR11prop_by_aquatic_PGLS <- phylolm(Y~X, data = OR11prop_by_aquatic, phy=turtles, model = "BM")
OR11prop<-summary(OR11prop_by_aquatic_PGLS)
subfamily.analysis.df[8,7] <-OR11prop$coefficients[2,1]
subfamily.analysis.df[8,8] <-OR11prop$coefficients[2,2]
subfamily.analysis.df[8,9] <-OR11prop$coefficients[2,3]
subfamily.analysis.df[8,10] <-OR11prop$coefficients[2,4]
subfamily.analysis.df[8,12] <-OR11prop$r.squared
#OR12prop
OR12prop_by_aquatic <- data.frame(X = turtleorprop.subfam.df$aquatic_cont, Y = turtleorprop.subfam.df$OR12, row.names = turtleorprop.subfam.df$species)
OR12prop_by_aquatic_PGLS <- phylolm(Y~X, data = OR12prop_by_aquatic, phy=turtles, model = "BM")
OR12prop<-summary(OR12prop_by_aquatic_PGLS)
subfamily.analysis.df[9,7] <-OR12prop$coefficients[2,1]
subfamily.analysis.df[9,8] <-OR12prop$coefficients[2,2]
subfamily.analysis.df[9,9] <-OR12prop$coefficients[2,3]
subfamily.analysis.df[9,10] <-OR12prop$coefficients[2,4]
subfamily.analysis.df[9,12] <-OR12prop$r.squared
#OR13prop
OR13prop_by_aquatic <- data.frame(X = turtleorprop.subfam.df$aquatic_cont, Y = turtleorprop.subfam.df$OR13, row.names = turtleorprop.subfam.df$species)
OR13prop_by_aquatic_PGLS <- phylolm(Y~X, data = OR13prop_by_aquatic, phy=turtles, model = "BM")
OR13prop<-summary(OR13prop_by_aquatic_PGLS)
subfamily.analysis.df[10,7] <-OR13prop$coefficients[2,1]
subfamily.analysis.df[10,8] <-OR13prop$coefficients[2,2]
subfamily.analysis.df[10,9] <-OR13prop$coefficients[2,3]
subfamily.analysis.df[10,10] <-OR13prop$coefficients[2,4]
subfamily.analysis.df[10,12] <-OR13prop$r.squared
#OR14prop
OR14prop_by_aquatic <- data.frame(X = turtleorprop.subfam.df$aquatic_cont, Y = turtleorprop.subfam.df$OR14, row.names = turtleorprop.subfam.df$species)
OR14prop_by_aquatic_PGLS <- phylolm(Y~X, data = OR14prop_by_aquatic, phy=turtles, model = "BM")
OR14prop<-summary(OR14prop_by_aquatic_PGLS)
subfamily.analysis.df[11,7] <-OR14prop$coefficients[2,1]
subfamily.analysis.df[11,8] <-OR14prop$coefficients[2,2]
subfamily.analysis.df[11,9] <-OR14prop$coefficients[2,3]
subfamily.analysis.df[11,10] <-OR14prop$coefficients[2,4]
subfamily.analysis.df[11,12] <-OR14prop$r.squared
#OR51prop
OR51prop_by_aquatic <- data.frame(X = turtleorprop.subfam.df$aquatic_cont, Y = turtleorprop.subfam.df$OR51, row.names = turtleorprop.subfam.df$species)
OR51prop_by_aquatic_PGLS <- phylolm(Y~X, data = OR51prop_by_aquatic, phy=turtles, model = "BM")
OR51prop<-summary(OR51prop_by_aquatic_PGLS)
subfamily.analysis.df[12,7] <-OR51prop$coefficients[2,1]
subfamily.analysis.df[12,8] <-OR51prop$coefficients[2,2]
subfamily.analysis.df[12,9] <-OR51prop$coefficients[2,3]
subfamily.analysis.df[12,10] <-OR51prop$coefficients[2,4]
subfamily.analysis.df[12,12] <-OR51prop$r.squared
#OR52prop
OR52prop_by_aquatic <- data.frame(X = turtleorprop.subfam.df$aquatic_cont, Y = turtleorprop.subfam.df$OR52, row.names = turtleorprop.subfam.df$species)
OR52prop_by_aquatic_PGLS <- phylolm(Y~X, data = OR52prop_by_aquatic, phy=turtles, model = "BM")
OR52prop<-summary(OR52prop_by_aquatic_PGLS)
subfamily.analysis.df[13,7] <-OR52prop$coefficients[2,1]
subfamily.analysis.df[13,8] <-OR52prop$coefficients[2,2]
subfamily.analysis.df[13,9] <-OR52prop$coefficients[2,3]
subfamily.analysis.df[13,10] <-OR52prop$coefficients[2,4]
subfamily.analysis.df[13,12] <-OR52prop$r.squared
#OR55 not analyzed due to gene count

#BF corrections
#count
ORcount.rawpval.vector <- as.vector(subfamily.analysis.df$Count.raw.pval)
ORcount.bf.pval <- p.adjust(ORcount.rawpval.vector, method= "bonferroni", n= length(ORcount.rawpval.vector))
subfamily.analysis.df$Count.BF.pval <-ORcount.bf.pval
#proportion
ORprop.rawpval.vector <- as.vector(subfamily.analysis.df$Prop.raw.pval)
ORprop.bf.pval <- p.adjust(ORprop.rawpval.vector, method= "bonferroni", n= length(ORprop.rawpval.vector))
subfamily.analysis.df$Prop.BF.pval <- ORprop.bf.pval

#save data frame as csv
write.csv(subfamily.analysis.df, file = "OR_subfamily_table.csv", row.names = TRUE)

#Plots of significant subfamilies
OR6.plot <- ggplot(OR6prop_by_aquatic, aes(x= X, y = Y)) +
  geom_point() +
  labs(title = "OR6",
       x = "Aquatic Affinity",
       y = "Intact Repertoire Proportion") +
  ylim(0, 0.5) +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("OR6.prop.pdf", OR6.plot, 
       device = 'pdf', units = 'in', height = 4, width = 6.5)

OR11.plot <- ggplot(OR11prop_by_aquatic, aes(x= X, y = Y)) +
  geom_point() +
  labs(title = "OR11",
       x = "Aquatic Affinity",
       y = "Intact Repertoire Proportion") +
  ylim(0, 0.5) +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("OR11.prop.pdf", OR11.plot, 
       device = 'pdf', units = 'in', height = 4, width = 6.5)

OR52.plot <- ggplot(OR52prop_by_aquatic, aes(x= X, y = Y)) +
  geom_point() +
  labs(title = "OR52",
       x = "Aquatic Affinity",
       y = "Intact Repertoire Proportion") +
  ylim(0, 0.5) +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("OR52.prop.pdf", OR52.plot, 
       device = 'pdf', units = 'in', height = 5, width = 6.5)

