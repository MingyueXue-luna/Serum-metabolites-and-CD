library(ppcor)
#packageVersion("ppcor")
rm(list = ls())

setwd("C:/Users/74511/Desktop/Projects data/Serum_metabomic/Income new batch")
#setwd("C:/Users/74511/Desktop/Projects data/Quik analysis")

Serum <- read.csv(file = "63Raw_shannon_cor.csv")
#Serum <- read.csv(file = "Mested_NTF_02.csv")
#Serum$AS <- factor(Serum$AS, levels = c('0', '1', '2', '3','4','5'))
names(Serum)

for(i in c(6)){
  if(i==6){
    results = matrix(ncol=4, nrow=0)
    colnames(results) = c("Meta1","Meta2", "estimate", "pvalue")
  }
  print(i)
  
  for(j in c(7:69)){
    print(c(i,j))
    tab_intermed=cbind("Se_PC"=Serum[,i],
                       "Genes_PC"=Serum[,j],
                       Serum[c("Nestedgroup","Multiplex","Relation.To.Proband")])
    
    colnames(tab_intermed)
    
    p <- pcor.test(tab_intermed$Se_PC,tab_intermed$Genes_PC,tab_intermed[,c("Nestedgroup","Multiplex","Relation.To.Proband")],method="spearman")
    
    
    results = rbind(results,
                    c(names(Serum)[i],
                      names(Serum)[j],
                      p[1],
                      p[2]
                    )
    )
    
  }
}

write.csv(results, "Shannon_meta_6_19.csv")
#write.csv(results, "Meta_olink_cor.csv")


results<-read.csv("Shannon_meta_6_19.csv",header=T)
results$fdr<-(p.adjust(results$pvalue, method="fdr", n=length(results$pvalue)))

###################

#setwd("C:/Users/74511/Desktop/Projects data/Serum_metabomic/Income new batch/")
setwd("C:/Users/74511/Desktop/Projects data/Quik analysis/")

Serum <- read.csv(file = "Mested_NTF_02.csv")
#Serum$AS <- factor(Serum$AS, levels = c('0', '1', '2', '3','4','5'))
names(Serum)

for(i in c(5:417)){
  if(i==5){
    results = matrix(ncol=3, nrow=0)
    colnames(results) = c("Se_PC", "estimate", "pvalue")
  }
  print(i)
  
  tab_intermed=cbind("Se_PC"=Serum[,i],
                      Serum[c("Nested_cohort","TNF.R2")])
    
  colnames(tab_intermed)
    
  p <- pcor.test(tab_intermed$Se_PC,tab_intermed$TNF.R2,tab_intermed[,c("Nested_cohort")],method="spearman")
  
  results = rbind(results,
                  c(names(Serum)[i],
                    p[1],
                    p[2]
                  )
    )
    
}



write.csv(results, "TNF_16S.csv")


results<-read.csv("TNF_16S.csv",header=T)
results$fdr<-(p.adjust(results$pvalue, method="fdr", n=length(results$pvalue)))

########识别top 20 16s是相关危险还是保护因素#######
setwd("C:/Users/74511/Desktop/My writing/Serum metabolites/Metabolomic paper_WT_8_28/Metabolomic paper")


Serum <- read.csv(file = "qantiles.csv")
#Serum$AS <- factor(Serum$AS, levels = c('0', '1', '2', '3','4','5'))
str(Serum)

for(i in c(2:42)){
  if(i==2){
    results = matrix(ncol=3, nrow=0)
    colnames(results) = c("Se_PC", "estimate", "pvalue")
  }
  print(i)
  
  tab_intermed=cbind("Se_PC"=Serum[,i],
                     Serum[c("Taxon")])
  
  colnames(tab_intermed)
  
  p <- pcor.test(tab_intermed$Se_PC,tab_intermed$Taxon,tab_intermed[,c()],method="spearman")
  
  results = rbind(results,
                  c(names(Serum)[i],
                    p[1],
                    p[2]
                  )
  )
  
}



write.csv(results, "quitail_corre.csv")


results<-read.csv("META_CRP.csv",header=T)
results$fdr<-(p.adjust(results$pvalue, method="fdr", n=length(results$pvalue)))


