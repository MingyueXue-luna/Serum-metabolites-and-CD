rm(list = ls())

require(knitr)
require(tidyverse)
require(RColorBrewer)
require(cowplot)
require(reshape2)
require(ggdendro)
require(vegan)
require(ape)
library(compositions)
library(PCAtools)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("PCAtools")


#help(PCAtools)

setwd("C:/Users/74511/Desktop/Projects data/Serum_metabomic/Income new batch/PCA_all/")
#rm(list=ls())

metadata_my <- read.csv("Meta_PCA_4181S1_3819S1.csv",row.names = 1)
Species <- read.csv("Serum_pc_use.csv", row.names = 1)

#metadata_my <- read.csv("Meta_Olink.csv",row.names = 1)
#Species <- read.csv("Olink.csv", row.names = 1)

str(metadata_my)

# filter the expression data to match the samples in our pdata
Species <- Species[,which(colnames(Species) %in% rownames(metadata_my))]

all(colnames(Species) == rownames(metadata_my))

#一些转换
#Species1 <- t(Species)
#Species1 <- log(Species+1)
#Species1 <- clr(Species)
#Species1 <- as.matrix(Species1)
#求最佳PC
p <- pca(Species, metadata = metadata_my, removeVar = 0)

screeplot(p, 
          ylim = c(0, 50),
          title = "SCREE plot of principal component 1-10",
          components = getComponents(p, c(1:10)),
          axisLabSize = 16, titleLabSize = 18)

a <- p$rotated[1:10]
biplot(p, showLoadings = TRUE, lab = NULL)
write.csv(a,file = "serum_PCs_2023_11_26.csv")

horn <- parallelPCA(Species)
horn$n

elbow <- findElbowPoint(p$variance)
elbow


#PCoA
metabo_all_dist <- dist(t(Species))
metabo_all_pcoa <- data.frame(pcoa(metabo_all_dist)$vectors) #pca
write.csv(metabo_all_pcoa,file = "Serum_PCoAs_11_26.csv")
eigen_meta1 <- pcoa(metabo_all_dist)



eigen_meta <- pcoa(metabo_all_dist)$values$Eigenvalues #pca
percent_var_meta <- signif(eigen_meta/sum(eigen_meta), 4)*100
write.csv(percent_var_meta,file = "percent_var_SERUM.csv")

#metabo_all_pcoa <- rownames_to_column(metabo_all_pcoa, var = "X.SampleID")

#Species1 <- as.matrix(Species1)

metabo_all_dist1 <- vegdist(t(Species), method = "bray")


##PCOA
# 计算加权bray-curtis距离
#specise_dist <- vegdist(Species1, method="bray", binary=F)
#specise_pcoa <- cmdscale(specise_dist, k=10, eig=T)
#specise_pcoa_points <- as.data.frame(specise_pcoa$points)
#sum_eig <- sum(specise_pcoa$eig)
#eig_percent <- round(specise_pcoa$eig/sum_eig*100,1)



specise_pcoa_result <- cbind(metabo_all_pcoa, metadata_my)
head(specise_pcoa_result)

perm <- with(metadata_my, how(nperm = 9999, blocks = Nestedgroup))


species.div <- adonis(metabo_all_dist ~ Outcome, data = metadata_my,permutations = perm)
#species.div <- adonis2(Species1 ~ Outcome, data = metadata_my, permutations = 10000)
#specise_adonis <- paste0("Species level adonis R2: ",round(species.div$R2,3), "; P-value: ", round(species.div$`Pr(>F)`,4))

plot_p <- species.div$aov.tab$`Pr(>F)`[1]
specise_pcoa_result$Nestedgroup <- as.factor(specise_pcoa_result$Nestedgroup)

library(ggplot2)

#画图65%的置信椭圆
food_micro <- ggplot(specise_pcoa_result, aes(x = Axis.1, y = Axis.2, fill = Nestedgroup)) +
  geom_point(size = 4.2, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = Outcome)) +
  stat_ellipse(color = "light grey", type = "norm", linetype = 2, level = 0.75, alpha = 0.8,lwd = 0.75) +
  scale_color_manual(values = c("#f58231", "#5f86b7")) +
  theme_classic() +
  #geom_line(aes(x= V1, y=V2, group=Con_1), col = "darkgrey", alpha = 0.6) +
  #geom_line(aes(x= V1, y=V2, group=Con_2), col = "darkgrey", alpha = 0.6) +
  #geom_line(aes(x= V1, y=V2, group=Con_3), col = "darkgrey", alpha = 0.6) +
  #geom_line(aes(x= V1, y=V2, group=Con_4), col = "darkgrey", alpha = 0.6) +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 6),
        axis.title = element_text(size=12)) +
  guides(color = guide_legend(ncol = 5))+
  annotate("text", x = 3.5, y = -4, label = paste0("p-value=",round(plot_p,4)), size = 4.5) +
  xlab(paste0("PC 1 [",percent_var_meta[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var_meta[2],"%]"))
  
food_micro_leg <- get_legend(food_micro)
food_micro + theme(legend.position = "right",
                   legend.box.just = "left",
                   legend.box.background = element_rect(color="black", size=0.8),
                   legend.box.margin = margin(3,3,3,3))  
  

# 设置工作目录
setwd("C:/Users/74511/Desktop/Projects data/Serum_metabomic/Income new batch/Procrust/")
###Procrust
CD<- read.csv("data_Auto_CD(2).csv",row.names = 1)
MC1 <- read.csv("data_Auto_MC1-1.csv", row.names = 1)

CD_dist <- dist(CD)
MC1_dist <- dist(MC1)


pcoa_CD <- as.data.frame(pcoa(CD_dist)$vectors)
pcoa_MC1 <- as.data.frame(pcoa(MC1_dist)$vectors)


# procrustes
pro1 <- procrustes(pcoa_CD, pcoa_MC1)
pro_test1 <- protest(pcoa_CD, pcoa_MC1, perm = 9999)
eigen1 <- sqrt(pro1$svd$d)
percent_var1 <- signif(eigen1/sum(eigen1), 4)*100
beta_pro1 <- data.frame(pro1$X)
trans_pro1 <- data.frame(pro1$Yrot)
beta_pro1$UserName <- rownames(beta_pro1)
beta_pro1$type <- "CD"
trans_pro1$UserName <- rownames(trans_pro1)
trans_pro1$type <- "MC1"
colnames(trans_pro1) <- colnames(beta_pro1)
pval1 <- signif(pro_test1$signif, 1)
plot1 <- rbind(beta_pro1, trans_pro1)


CD_MC1 <- ggplot(plot1) +
  geom_point(size = 4, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = type)) +
  scale_color_manual(values = c("#f58231","#CDCD00")) +
  theme_classic() +
  geom_line(aes(x= Axis.1, y=Axis.2, group=UserName), col = "darkgrey", alpha = 0.6) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position = 'bottom',
        axis.text = element_text(size=4),
        axis.title = element_text(size=9),
        aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  annotate("text", x = 10, y = -20, label = paste0("p-value=",pval1), size = 4) +
  xlab(paste0("PC 1 [",percent_var1[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var1[2],"%]"))
CD_MC_leg1 <- get_legend(CD_MC1)
p1 <- CD_MC1 + theme(legend.position = "right",
                     legend.box.just = "left",
                     legend.box.background = element_rect(color="black", size=0.9),
                     legend.box.margin = margin(3,3,3,3))


###MC2
CD<- read.csv("data_Auto_CD(2).csv",row.names = 1)
MC2 <- read.csv("data_Auto_MC2-1.csv", row.names = 1)

CD_dist <- dist(CD)
MC2_dist <- dist(MC2)


pcoa_CD <- as.data.frame(pcoa(CD_dist)$vectors)
pcoa_MC2 <- as.data.frame(pcoa(MC2_dist)$vectors)


# procrustes
pro2 <- procrustes(pcoa_CD, pcoa_MC2)
pro_test2 <- protest(pcoa_CD, pcoa_MC2, perm = 9999)
eigen2 <- sqrt(pro2$svd$d)
percent_var2 <- signif(eigen2/sum(eigen2), 4)*100
beta_pro2 <- data.frame(pro2$X)
trans_pro2 <- data.frame(pro2$Yrot)
beta_pro2$UserName <- rownames(beta_pro2)
beta_pro2$type <- "CD"
trans_pro2$UserName <- rownames(trans_pro2)
trans_pro2$type <- "MC2"
colnames(trans_pro2) <- colnames(beta_pro2)
pval2 <- signif(pro_test2$signif, 1)
plot2 <- rbind(beta_pro2, trans_pro2)
#head(plot)
#"#bfef45"
CD_MC2 <- ggplot(plot2) +
  geom_point(size = 4, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = type)) +
  scale_color_manual(values = c("#f58231","#CD69C9")) +
  theme_classic() +
  geom_line(aes(x= Axis.1, y=Axis.2, group=UserName), col = "darkgrey", alpha = 0.6) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position = 'bottom',
        axis.text = element_text(size=4),
        axis.title = element_text(size=9),
        aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  annotate("text", x = 8, y = -20, label = paste0("p-value=",pval2), size = 4) +
  xlab(paste0("PC 1 [",percent_var2[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var2[2],"%]"))
CD_MC_leg2 <- get_legend(CD_MC2)
p2 <- CD_MC2 + theme(legend.position = "right")

#MC3
###Procrust
CD3<- read.csv("data_Auto_CD_MC3-1.csv",row.names = 1)
MC3 <- read.csv("data_Auto_MC3-1.csv", row.names = 1)

CD_dist3 <- dist(CD3)
MC3_dist <- dist(MC3)


pcoa_CD3 <- as.data.frame(pcoa(CD_dist3)$vectors)
pcoa_MC3 <- as.data.frame(pcoa(MC3_dist)$vectors)


pro3 <- procrustes(pcoa_CD3, pcoa_MC3)
pro_test3 <- protest(pcoa_CD3, pcoa_MC3, perm = 9999)
eigen3 <- sqrt(pro3$svd$d)
percent_var3 <- signif(eigen3/sum(eigen3), 4)*100
beta_pro3 <- data.frame(pro3$X)
trans_pro3 <- data.frame(pro3$Yrot)
beta_pro3$UserName <- rownames(beta_pro3)
beta_pro3$type <- "CD"
trans_pro3$UserName <- rownames(trans_pro3)
trans_pro3$type <- "MC3"
colnames(trans_pro3) <- colnames(beta_pro3)
pval3 <- signif(pro_test3$signif, 1)
plot3 <- rbind(beta_pro3, trans_pro3)
head(plot)
CD_MC3 <- ggplot(plot3) +
  geom_point(size = 4, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = type)) +
  scale_color_manual(values = c("#f58231", "#5a2071")) +
  theme_classic() +
  geom_line(aes(x= Axis.1, y=Axis.2, group=UserName), col = "darkgrey", alpha = 0.6) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position = 'bottom',
        axis.text = element_text(size=4),
        axis.title = element_text(size=9),
        aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  annotate("text", x = 10, y = -21, label = paste0("p-value=",pval3), size = 4) +
  xlab(paste0("PC 1 [",percent_var3[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var3[2],"%]"))
CD_MC_leg3 <- get_legend(CD_MC3)
p3 <- CD_MC3 + theme(legend.position = "right")

#MC4
###Procrust
CD4<- read.csv("data_Auto_CD_MC4-1.csv",row.names = 1)
MC4 <- read.csv("data_Auto_MC4-1.csv", row.names = 1)

CD_dist4 <- dist(CD4)
MC4_dist <- dist(MC4)


pcoa_CD4 <- as.data.frame(pcoa(CD_dist4)$vectors)
pcoa_MC4 <- as.data.frame(pcoa(MC4_dist)$vectors)


pro4 <- procrustes(pcoa_CD4, pcoa_MC4)
pro_test4 <- protest(pcoa_CD4, pcoa_MC4, perm = 9999)
eigen4 <- sqrt(pro4$svd$d)
percent_var4 <- signif(eigen4/sum(eigen4), 4)*100
beta_pro4 <- data.frame(pro4$X)
trans_pro4 <- data.frame(pro4$Yrot)
beta_pro4$UserName <- rownames(beta_pro4)
beta_pro4$type <- "CD"
trans_pro4$UserName <- rownames(trans_pro4)
trans_pro4$type <- "MC4"
colnames(trans_pro4) <- colnames(beta_pro4)
pval4 <- signif(pro_test4$signif, 1)
plot4 <- rbind(beta_pro4, trans_pro4)
#head(plot)
CD_MC4 <- ggplot(plot4) +
  geom_point(size = 4, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = type)) +
  scale_color_manual(values = c("#f58231", "#96CDCD")) +
  theme_classic() +
  geom_line(aes(x= Axis.1, y=Axis.2, group=UserName), col = "darkgrey", alpha = 0.6) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position = 'bottom',
        axis.text = element_text(size=4),
        axis.title = element_text(size=9),
        aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  annotate("text", x = 10, y = -20, label = paste0("p-value=",pval4), size = 4) +
  xlab(paste0("PC 1 [",percent_var4[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var4[2],"%]"))
CD_MC_leg4 <- get_legend(CD_MC4)
p4 <- CD_MC4 + theme(legend.position = "right")

#combind the pictures
library(cowplot)
library(patchwork)

p_all <- plot_grid(p1, p2, p3, p4, ncol = 2)
p_all
ggsave("plot2.pdf", plot = p_all, height = 9.5, width = 10)
