library(survival)
packageVersion("survival")

rm(list = ls())
setwd("C:/Users/74511/Desktop/Projects data/Serum_metabomic/Income new batch/Cor_16s_Olink")
#setwd("C:/Users/74511/Desktop/Projects data/Serum_metabomic/Income new batch/PC_cor")
#Shotgun <- read.csv("Auto_nomalized_regression.csv")
#Shotgun <- read.csv("Sensitivity_regression_shorter_follow.csv")
Shotgun <- read.csv("63Raw_top20_16s.csv")



Shotgun$Multiplex <- as.factor(Shotgun$Multiplex)
Shotgun$Relation.To.Proband <- as.factor(Shotgun$Relation.To.Proband)
#Shotgun$AS <- as.factor(Shotgun$AS)
names(Shotgun)


c-index




results<-read.csv("../Outcome new batch/Regression_616.csv",header=T)
results<-read.csv("try.csv",header=T)
results$fdr<-(p.adjust(results$p, method="fdr", n=length(results$p)))
results$fdr1<-(p.adjust(results$p, method="bonferroni", n=length(results$p)))

write.csv(results, "try.csv")
