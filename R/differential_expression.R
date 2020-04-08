
#ggforce for ggsina
if(install){
  install.packages("BiocManager")
  BiocManager::install("VGAM")
#  BiocManager::install("ggplot2")
#  BiocManager::install("gplots")
  #BiocManager::install("jsonlite")
  #BiocManager::install("gridExtra")
  #BiocManager::install("GGally")
 # BiocManager::install("binom")
 BiocManager::install("biomaRt")
}
library(VGAM)
library(biomaRt)

.findFile<-function(path, file, exact = T){
  for(i in 1:length(path)){
    if(exact){
      res = paste(path[i],file,sep="/") 
      if(file.exists(res)) { 
        return(res)
      }
    }else{
      files = grep(file, dir(path[i]) , v=T)
      if(length(files)==1){
        return(paste(path[i], files,sep="/"))
      }
    } 
  }
}

src = c("../../R" , "~/github/npTranscript/R" )
source(.findFile(src, "diff_expr_functs.R"))

#library(binom)
#library(VGAM)
files = dir()


a = read.table(files[1], head=T)
b = read.table(files[2], head=T)
geneID = as.character(a$Geneid)
x = a[,7]
y = b[,7]
DE1 = DEgenes(geneID,x,y,log=F);
DE2 = DEgenes(geneID,y,x,log=F);

genenames1 = as.character(DE1[which(DE1[,1]<1e-10),]$genenames)
genenames2 = as.character(DE2[which(DE1[,2]<1e-10),]$genenames)

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "csabaeus_gene_ensembl", 
                   mirror = "asia")

desc1 = getDescr(genenames1, mart)
desc2 = getDescr(genenames2, mart)

#attr= listAttributes(mart);


