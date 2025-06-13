rm(list = ls())
Metabo <- read.csv("data_range_cdprs.csv")


  
}

write.csv(results, "GEE_cdprs.csv")

results<-read.csv("GEE_cdprs.csv",header=T)

results$fdr<-(p.adjust(results$pvalue, method="fdr", n=length(results$pvalue)))


