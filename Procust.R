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

setwd("C:/Users/74511/Desktop/Projects data/Serum_metabomic/Income new batch/Procrust/")
#rm(list=ls())

  
###Procrust
CD<- read.csv("Sample_CD.csv",row.names = 1)
MC1 <- read.csv("Sample_MC.csv", row.names = 1)

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


