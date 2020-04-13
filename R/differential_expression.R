
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
chrs = unlist(lapply(a$Chr,function(x) strsplit(as.character(x),";")[[1]][1]))

df = data.frame(control = a[,7],infected = b[,7],chrs = chrs)
names(df) = c("control", "infected", "chrom")
dimnames(df)[[1]] = geneID
DE1 = DEgenes(df,log=F,);
DE2 = DEgenes(df,log=F,inds = c(2,1));

#CHECK DISTRIBUTION
.qqplot(DE1,nme="FDR")
.qqplot(DE2, nme="FDR")


mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "csabaeus_gene_ensembl", 
                   mirror = "uswest") #asia useast

DE1 = getDescr(DE1, mart,thresh = 1e-5)
DE2 = getDescr(DE2, mart,thresh = 1e-5)

#FOLLOWING COMMAND OFTEN FAILS!  NEED TO RETRY FEW TIMES with different mirrors
goObjs = getGoIDs(  unique( as.character(DE1$geneID)),mart)
#chromObjs = getChromIDs(  unique( as.character(DE1$geneID)),mart)

#goids2 = getGoIDs(DE2,mart,1e-10)

sigGo1 = list()
sigGo2 = list()
for(i in 1:length(goObjs)){
 sigGo1[[i]] = findSigGo_(goObjs[[i]],DE1, fdr_thresh = 1e-5, go_thresh = 1e-4);
 sigGo2[[i]] = findSigGo_(goObjs[[i]],DE2, fdr_thresh = 1e-5, go_thresh = 1e-4);
 
}

names(sigGo1) = names(goObjs)
names(sigGo2) = names(goObjs)

sigChr1 = findSigChrom(DE1,fdr_thresh = 1e-5, go_thresh = 1e-4)
sigChr2 = findSigChrom(DE2,fdr_thresh = 1e-5, go_thresh = 1e-4)

 


#attr= listAttributes(mart);

findGenesByChrom(DE2,"M", fdr_thresh = 1e-10)

findGenes(goObjs[[1]],DE2,"GO:0001968", fdr_thresh = 1e-10)
findGenes(goObjs[[1]],DE2,"GO:0005372", fdr_thresh = 1e-10)
findGenes(goObjs[[1]],DE2,"GO:0006412", fdr_thresh = 1e-10)
findGenes(goObjs[[1]],DE2,"GO:0002020", fdr_thresh = 1e-5)

findGenes(goObjs[[1]],DE1,"GO:0006412", fdr_thresh = 1e-5)

#findGenes(chromObjs,DE1,"X", fdr_thresh = 1e-5)


