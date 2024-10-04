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

###clogit function
clogit1 = clogit(CD ~ PC1+strata(Nestedgroup)+FCP,
                 data=Shotgun)

res1 <- summary(clogit1)
pnorm(summary(clogit1)$coefficients[1,4], lower.tail = F)


m <- glm(Outcome ~ X35, family = binomial(),data=Shotgun)

pnorm(summary(m)$coefficients[2,3], lower.tail = F)


res1$coefficients[1,5]
res1$coefficients[1,1]
res1$coefficients[1,3]
res1$concordance[1:1]   #to get the c-index and relevant SE
res1$concordance[2:2]

res1$coefficients[1,5]
res1$coefficients[1,2]
paste(
  res1$conf.int[1,3], 
  res1$conf.int[1,4], 
  sep="-")
res1$concordance[1]



clogit1 = lm(Outcome ~ X35+Multiplex+Relation.To.Proband,
                 data=Shotgun)
res <- summary(clogit1)


2*pt(-abs(coef(res)[, 3]), clogit1$df)

Shotgun1$Butyricimonas_synergistica1 = factor(ifelse(Shotgun1$Butyricimonas_synergistica == 0, 0, 1),ordered = F)


OR  p 
c-index



####clogit loop####
#zeroThresh = 0.1

#for(i in c(6:1031)){
for(i in c(70:89)){
  if(i==70){
    results = matrix(ncol=6, nrow=0)
    colnames(results) = c("taxa","zero", "pvalue", "OR", "OR_CI_95pc", "C_index_model"
                          )
  }
  print(i)
  #,"Multiplex","Relation.To.Proband"
  tab_intermed=cbind("taxa"=Shotgun[,i], 
                     Shotgun[c("Outcome", "Nestedgroup","Relation.To.Proband","Multiplex"
                               )])
  
  colnames(tab_intermed)
  
  zeroProp = length(tab_intermed$taxa[tab_intermed$taxa == 0])/nrow(tab_intermed)
  #if(zeroProp>0.85) next
  #if(zeroProp>zeroThresh){
  
  #tab_intermed$taxa = factor(ifelse(tab_intermed$taxa == 0, 0, 1),ordered = F)
  
  #}

  ####perform survival analysis
  #+Multiplex+Relation.To.Proband
  clogit1 = clogit(Outcome ~ scale(taxa)+strata(Nestedgroup)+Multiplex+Relation.To.Proband,
            data=tab_intermed)
  
  ss=summary(clogit1)
  
  results = rbind(results,
                  c(names(Shotgun)[i],
                    zeroProp,
                    ss$coefficients[1,5],
                    ss$coefficients[1,2],
                    paste(
                      ss$conf.int[1,3], 
                      ss$conf.int[1,4], 
                      sep="-"),
                    paste(round(ss$concordance[1],2), round(ss$concordance[2],2), sep="__SD_")
                    )
  )
  
}


results
write.csv(results, "top20_16s_logistic-1.csv")

results<-read.csv("../Outcome new batch/Regression_616.csv",header=T)
results<-read.csv("try.csv",header=T)
results$fdr<-(p.adjust(results$p, method="fdr", n=length(results$p)))
results$fdr1<-(p.adjust(results$p, method="bonferroni", n=length(results$p)))

write.csv(results, "try.csv")
