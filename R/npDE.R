
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
library(rhdf5)
#library(seqinr)



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
featureCounts = F
chroms = NULL

transcripts = readTranscriptHostAll(grep("transcripts.txt", dir(), v=T), 
                                combined_depth_thresh = 100)
 
 # dfs = list()
 
control_inds = 1
infected_inds = 2
DE1 = DEgenes(df, control_inds, infected_inds,edgeR = F, reorder=T);
.qqplot(DE1$pvals,min.p= 1e-200)

#this calculates pvalues for base-level error rates
depth=.readH5All(transcripts, chroms = attr(transcripts,"chroms"),depth_thresh = 1000)

##this is visualisation
pv_inds = grep("pv", names(depth))
.qqplot(depth$pv1,min.p=1e-100)
.qqplot(depth$pv2,min.p=1e-100)
hist(depth[depth$pv1<1e-5,]$base)
hist(depth[depth$pv2<1e-5,]$base)
#minp = apply(pv_a[,pv_inds],1,min,na.rm=T)
#pv_a[which(minp==min(minp)),]
.vis(depth,i=1)
.vis(depth,i=2)




#DE2 = DEgenes(df,control_inds, infected_inds, log=F,edgeR = F, reorder=T);

#CHECK DISTRIBUTION (OPTIONAL)
.qqplot(DE1$FDR)
#.qqplot(DE1$FDR2)


mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "csabaeus_gene_ensembl", 
                   mirror = "uswest") #asia useast

DE1 = getDescr(DE1, mart,thresh = 1e-10, prefix=prefix)

#FOLLOWING COMMAND OFTEN FAILS!  NEED TO RETRY FEW TIMES with different mirrors
geneNames = unique(grep(prefix ,DE1$geneID,v=T))
goObjs = getGoIDs( geneNames,mart)

sigGo1 = lapply(goObjs ,findSigGo_,DE1, fdr_thresh = 1e-5, go_thresh = 1e-4, lessThan=T)
sigGo2 = lapply(goObjs ,findSigGo_,DE1, fdr_thresh = 1e-5, go_thresh = 1e-4, lessThan=F)

names(sigGo1) = names(goObjs)
names(sigGo2) = names(goObjs)


sigChr1 = findSigChrom(DE1,fdr_thresh = 1e-5, go_thresh = 1e-4, lessThan =T)
sigChr2 = findSigChrom(DE1,fdr_thresh = 1e-5, go_thresh = 1e-4, lessThan =F)

DE1_na = DE1[is.na(DE1$type),]
.qqplot(DE1_na$FDR)
sigChr1 = findSigChrom(DE1_na,fdr_thresh = 1e-5, go_thresh = 1e-4, lessThan =T)
sigChr2 = findSigChrom(DE1_na,fdr_thresh = 1e-5, go_thresh = 1e-4, lessThan =F)
findGenesByChrom(DE1_na,"26", fdr_thresh = 1e-5)
findGenesByChrom(DE1_na,"AQIB01159108.1", fdr_thresh = 1e-5)


#attr= listAttributes(mart);

##EXPLORATION
findGenesByChrom(DE1,"MT", fdr_thresh = 1e-5)

#go_categories = list("GO:0001968","GO:0005372","GO:0006412","GO:0002020","GO:0051702")

go_categories2 = as.character(sigGo2[[1]][,1])
getGoGenes(as.list(go_categories2),goObjs, lessThan=NULL, fdr_thresh = 1e-5)

go_categories1 = as.character(sigGo1[[1]][,1])
getGoGenes(as.list(go_categories1),goObjs, lessThan=NULL, fdr_thresh = 1e-5)




#findGenes(chromObjs,DE1,"X", fdr_thresh = 1e-5)
#ACE2  ENSCSAG00000014921

getGene<-function(nme,DE){
  i = which(DE$geneID== nme)
  DE[i,]
}
genes = list(ACE2 = "ENSCSAG00000014921", 
             MMP9 = "ENSCSAG00000014722",
             TMPRSS2 ="ENSCSAG00000008229", cathepsin_L="ENSCSAG00000007387", BSG="ENSCSAG00000013009",
            FURIN, "ENSCSAG00000017046")

lapply(genes, getGene, DE1)
