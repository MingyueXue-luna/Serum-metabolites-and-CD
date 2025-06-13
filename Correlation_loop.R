library(ppcor)
#packageVersion("ppcor")
rm(list = ls())

setwd("C:/Users/74511/Desktop/Projects data/Serum_metabomic/Income new batch")
#setwd("C:/Users/74511/Desktop/Projects data/Quik analysis")

Serum <- read.csv(file = "63Raw_shannon_cor.csv")
#Serum <- read.csv(file = "Mested_NTF_02.csv")
#Serum$AS <- factor(Serum$AS, levels = c('0', '1', '2', '3','4','5'))
names(Serum)




results<-read.csv("META_CRP.csv",header=T)
results$fdr<-(p.adjust(results$pvalue, method="fdr", n=length(results$pvalue)))


