library(biomaRt)



if(!is.null(dataset)){
mart <- useEnsembl(biomart = "ensembl", 
                   dataset =dataset, 
                   mirror = mirror) #asia useast
}else{
  mart=NULL
}


desciGT = .getDescEnrich(DE2,mart,thresh = 1e-5, go_thresh=1e-2,nme="pvals1",nme2="pos1M")
desciLT = .getDescEnrich(DE2,mart,thresh = 1e-5, go_thresh=1e-2,nme="pvals2",nme2="pos1M")

if(!.is.null(mart)){
  DE_sig1 = getDescr(DE_sig1, mart,thresh = 1e-3, prefix=prefix)
  
  DE_sig2 = getDescr(DE_sig2, mart,thresh = 1e-3, prefix=prefix)
}


DE1 = getDescr(DE1, mart,thresh = 1e-10, prefix=prefix)

#FOLLOWING COMMAND OFTEN FAILS!  NEED TO RETRY FEW TIMES with different mirrors
geneNames = unique(grep(prefix ,DE1$geneID,v=T))
goObjs = getGoIDs( geneNames,mart)

sigGo1 = lapply(goObjs ,findSigGo_,DE1, fdr_thresh = 1e-5, go_thresh = 1e-4, nme="pvals1")
sigGo2 = lapply(goObjs ,findSigGo_,DE1, fdr_thresh = 1e-5, go_thresh = 1e-4, nme="pvals2")

names(sigGo1) = names(goObjs)
names(sigGo2) = names(goObjs)




.qqplot(DE1$FDR1)
sigChr1 = findSigChrom(DE1,thresh = 1e-5, go_thresh = 1e-4, nme="FDR1",nme2="chrs")
sigChr2 = findSigChrom(DE1,thresh = 1e-5, go_thresh = 1e-4, nme="FDR2", nme2="chrs")
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
