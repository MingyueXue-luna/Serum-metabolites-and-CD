rm(list = ls())
Metabo <- read.csv("data_range_cdprs.csv")

Metabo$Multiplex <- as.factor(Metabo$Multiplex)
Metabo$Relation.To.Proband <- as.factor(Metabo$Relation.To.Proband)
str(Metabo)

library(geepack)
library(stringr)
#library(survival)

gee.glm3<-geeglm(Log_crp~ X100020913 + Multiplex + Relation.To.Proband, id=Nestedgroup,
                 family=gaussian(), corstr = "exchangeable",
                 dat=Metabo) #linear regression

ss=summary(gee.glm3)




####GEE regression loop####

for( i in 13:75){
  if(i==13){
    results = matrix(ncol=5, nrow=0)
    colnames(results) = c("taxa", "zero", "pvalue", "estimate", "se")
  }
  print(i)
  tab_intermed=cbind("taxa"=Metabo[,i],
                     Metabo[c("Nestedgroup","Multiplex","Relation.To.Proband","FCP","Log10.LMR","CD_prs","Log_crp","AS")]
  )
  #
  colnames(tab_intermed)
  
  zeroProp = length(tab_intermed$taxa[tab_intermed$taxa == 0])/nrow(tab_intermed)
  
  #Perform glm (logistic and linear)###
  # gee.glm<-geeglm(taxa ~ CD_prs + Multiplex + Relation.To.Proband, id=Nestedgroup,
  #                 family=binomial(link = "logit"), corstr = "exchangeable",
  #                 dat=tab_intermed)
  # 
  # ss=summary(gee.glm)
  
  gee.glm3<-geeglm(taxa~ CD_prs + Multiplex + Relation.To.Proband, id=Nestedgroup,
                   family=gaussian(), corstr = "exchangeable",
                   dat=tab_intermed) #linear regression

  ss=summary(gee.glm3)
  # 
  results=rbind(results,
                c(names(Metabo)[i],
                  zeroProp,
                  ss$coefficients[2,4], #p-val
                  ss$coefficients[2,1], #estimte
                  ss$coefficients[2,2] #std.err
                  )
  )
  
}

write.csv(results, "GEE_cdprs.csv")

results<-read.csv("GEE_cdprs.csv",header=T)

results$fdr<-(p.adjust(results$pvalue, method="fdr", n=length(results$pvalue)))


