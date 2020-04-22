
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
 BiocManager::install("edgeR")
}
library(VGAM)
library(biomaRt)
library(edgeR)
library(stats)

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

src = c( "../../R" , "~/github/npTranscript/R" )
source(.findFile(src, "diff_expr_functs.R"))
prefix = "ENSC"; #for vervet monkey

files = dir()
featureCounts = T
if(!featureCounts){
  infilesT = grep("transcripts", dir(), v=T)
  df=.readTranscriptsHost(infilesT)
}else{
files = grep("featurecount", dir(), v=T)
df = .readFeatureCounts(grep("featurecount", dir(), v=T))
}


control_inds = 1
infected_inds = 2
DE1 = DEgenes(df, control_inds, infected_inds,log=F, edgeR = F, reorder=F);
DE2 = DEgenes(df,control_inds, infected_inds, log=F,edgeR = F, reorder=T);

#CHECK DISTRIBUTION (OPTIONAL)
.qqplot(DE1$FDR)
.qqplot(DE2$FDR)


mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "csabaeus_gene_ensembl", 
                   mirror = "uswest") #asia useast

DE1 = getDescr(DE1, mart,thresh = 1e-10, prefix=prefix)
DE2 = getDescr(DE2, mart,thresh = 1e-10, prefix=prefix)

#FOLLOWING COMMAND OFTEN FAILS!  NEED TO RETRY FEW TIMES with different mirrors
geneNames = unique(grep(prefix ,DE1$geneID,v=T))
goObjs = getGoIDs( geneNames,mart)
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

findGenesByChrom(DE1,"MT", fdr_thresh = 1e-5)
findGenesByChrom(DE2,"MT", fdr_thresh = 1e-5)


findGenes(goObjs[[1]],DE2,"GO:0001968", fdr_thresh = 1e-10)
findGenes(goObjs[[1]],DE2,"GO:0005372", fdr_thresh = 1e-10)
findGenes(goObjs[[1]],DE2,"GO:0006412", fdr_thresh = 1e-10)
findGenes(goObjs[[1]],DE2,"GO:0002020", fdr_thresh = 1e-5)

findGenes(goObjs[[1]],DE1,"GO:0006412", fdr_thresh = 1e-5)

#findGenes(chromObjs,DE1,"X", fdr_thresh = 1e-5)
#ACE2  ENSCSAG00000014921

getGene<-function(nme,DE){
  i = which(DE$leftGene== nme)
  DE[i,]
}
genes = list(ACE2 = "ENSCSAG00000014921", 
             MMP9 = "ENSCSAG00000014722",
             TMPRSS2 ="ENSCSAG00000008229", cathepsin_L="ENSCSAG00000007387", BSG="ENSCSAG00000013009" )

lapply(genes, getGene, DE1)
lapply(genes, getGene, DE2)
