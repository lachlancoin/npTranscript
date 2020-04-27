
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
#set up mart
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "csabaeus_gene_ensembl", 
                   mirror = "uswest") #asia useast

files = dir()
featureCounts = F
chroms = NULL
control_inds = 1
infected_inds = 2

##READ TRANSCRIPT DATA
transcripts = readTranscriptHostAll(grep("transcripts.txt", dir(), v=T),  combined_depth_thresh = 100)
                             
##find DE genes
DE1 = DEgenes(transcripts, control_inds, infected_inds,edgeR = F, reorder=F);
pos1M = apply(cbind(as.character(DE1$chrs), as.character(round(DE1$start/1e5))),1,paste,collapse=".")
DE2 = cbind(DE1,pos1M)
sigChrLT = findSigChrom(DE2,fdr_thresh = 1e-5, go_thresh = 1e-2, nme="pvals2",nme2="pos1M")
print(sigChrLT);
print(sigChrGT)
desciGT = .getDescEnrich(DE2,fdr_thresh = 1e-5, go_thresh=1e-2,nme="pvals1",nme2="pos1M")
desciLT = .getDescEnrich(DE2,fdr_thresh = 1e-5, go_thresh=1e-2,nme="pvals2",nme2="pos1M")




pdf("qq.pdf")
.qqplot(DE1$pvals1, min.p= 1e-200,main="both")
.qqplot(DE1$pvals2, min.p= 1e-200,main="infected less",add=T)
.vis(DE1,i=1,min.p=1e-50)
 .vis(DE1,i=2,min.p=1e-50)

dev.off()

#this calculates pvalues for base-level error rates
depth=.readH5All(transcripts, chroms = attr(transcripts,"chroms"),depth_thresh = 200)

##this is visualisation
  pdf("error_associations.pdf")
  pv_inds = grep("pv", names(depth))
  .qqplot(depth$pv1,min.p=1e-100)
  .qqplot(depth$pv2,min.p=1e-100)
   hist(depth[depth$pv1<1e-5,]$base)
   hist(depth[depth$pv2<1e-5,]$base)
  .vis(depth,i=1,min.p=1e-100)
  .vis(depth,i=2,min.p=1e-100)
  dev.off()

chr_inds=attr(depth,"chr_inds")
chroms = attr(depth,"chroms")
chroms1 = names(chroms)
mi = match(chr_inds, chroms)
mi1 = unlist(lapply(mi, function(x) chroms1[x]))
pos100k = apply(cbind(as.character(mi1), as.character(round(depth$pos/1e5))),1,paste,collapse=".")

#pos100k = apply(cbind(chr_inds,round(depth$pos/1e5)),1,paste,collapse=".")
depth1 = cbind(depth, pos100k)
#chr_inds[which(depth$pv2<1e-5)]
#chr_inds[which(depth$pv1<1e-5)]

sigChr1 = findSigChrom(depth1, fdr_thresh=1e-5, go_thresh=1e-3, nme="pv1", nme2="pos100k")
DE_sig1 = DE2[DE2$pos1M %in% sigChr1$chrs,]
sigChr2 = findSigChrom(depth1, fdr_thresh=1e-5, go_thresh=1e-3, nme="pv2", nme2="pos100k")
DE_sig2 = DE2[DE2$pos1M %in% sigChr2$chrs,]
DE_sig2 = getDescr(DE_sig2, mart,thresh = 1e-5, prefix=prefix)


#CHECK DISTRIBUTION (OPTIONAL)
.qqplot(DE1$FDR)
#.qqplot(DE1$FDR2)




DE1 = getDescr(DE1, mart,thresh = 1e-10, prefix=prefix)

#FOLLOWING COMMAND OFTEN FAILS!  NEED TO RETRY FEW TIMES with different mirrors
geneNames = unique(grep(prefix ,DE1$geneID,v=T))
goObjs = getGoIDs( geneNames,mart)

sigGo1 = lapply(goObjs ,findSigGo_,DE1, fdr_thresh = 1e-5, go_thresh = 1e-4, nme="pvals1")
sigGo2 = lapply(goObjs ,findSigGo_,DE1, fdr_thresh = 1e-5, go_thresh = 1e-4, nme="pvals2")

names(sigGo1) = names(goObjs)
names(sigGo2) = names(goObjs)




DE1_na = DE1[is.na(DE1$type),]
.qqplot(DE1_na$FDR)
sigChr1 = findSigChrom(DE1_na,fdr_thresh = 1e-5, go_thresh = 1e-4, lessThan =T)
sigChr2 = findSigChrom(DE1_na,fdr_thresh = 1e-5, go_thresh = 1e-4, lessThan =F)
findGenesByChrom(DE1_na,"26", fdr_thresh = 1e-5)
findGenesByChrom(DE1_na,"AQIB01159108.1", fdr_thresh = 1e-5)


findGenesByChrom(DE2,"MT", fdr_thresh = 1e-5)

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
