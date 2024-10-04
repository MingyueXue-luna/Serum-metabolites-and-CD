library(readr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(ggsignif)

setwd("C:/Users/74511/Desktop/Projects data/Serum_metabomic/Income714/")


rm(list = ls())

top_10 <- read_csv("top 10 serum.csv")


top_10_long <- melt(top_10, id = c("Outcome","Outcome1"))         # Reshaping iris data
head(top_10_long)


top_10_long$Outcome1 <- as.factor(top_10_long$Outcome1)

ggplot(top_10_long) +  # ggplot function
  geom_boxplot(aes(x = variable, y = value, colour=Outcome),
               # 间距调整
               #position=position_dodge(1),
               # 间距调整：
               width = 0.8,
               notch = TRUE,
               #outlier.shape = NA
               )+
  #geom_rect(aes(xmin=0.0, xmax=11, ymin=0.9, ymax = 1.1), fill="#eaeae0")+
  # 灰色竖线：
  #geom_vline(xintercept = c(1.5, 2.5, 3.5,4.5, 5.5, 6.5, 7.5, 8.5, 9.5), color = "#bcbdbf", alpha = 0.8)+
  # x轴标签：
  scale_x_discrete(labels = rep(c("Pre-CD", "MC"), 20))+
  scale_y_continuous(expand = c(0,0.1))+
  # 颜色：
  scale_fill_manual(values = c("#bcbdbf",  "#3675b6"))+
  scale_color_manual(values = c("#bcbdbf", "#3675b6"))+
  # 主题：
  theme_bw()+
  theme(panel.grid = element_blank(), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))+
  # 标题：
  ggtitle("Top 10 Serum Metabolits")+
  facet_grid(.~variable,space = "fixed", scales = "free_x")



#ggsave("boxplot.pdf", height = 4, width = 7)
  
  

  



