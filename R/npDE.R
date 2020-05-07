
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
read.gff(file, na.strings = c(".", "?"), GFF3 = TRUE)
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
#source(.findFile(src, "transcript_functions.R"))

source(.findFile(src, "diff_expr_functs.R"))

prefix = "ENSC"; #for vervet monkey
#set up mart
dataset="csabaeus_gene_ensembl"
if(!is.null(dataset)){
mart <- useEnsembl(biomart = "ensembl", 
                   dataset =dataset, 
                   mirror = "uswest") #asia useast
}else{
  mart=NULL
}
files = dir()
chroms = NULL
control_inds = 1
infected_inds = 2
binsize = 1e5
isVirus=  FALSE;
start_text = "start"
mergeByPosAndGene = FALSE
filter = NULL
if(isVirus){
  filter = list("type_nme"="5_3")
  mergeByPosAndGene = TRUE
  binsize = 100  #1e5 for euk
  start_text ="endBreak"

}

#seqfile= "../Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa.gz"
#fasta = read.fasta(seqfile)
#seqlen = length(fasta[[1]])


target= list( chrom="character", 
             leftGene="character", rightGene="character", start = "numeric", 
             end="numeric", ID="character", isoforms="numeric" ,type_nme="character")


##READ TRANSCRIPT DATA

 gfft = read.table("annotation.csv.gz", sep="\t", header=T, fill=T, quote='\"')

          

#"gene:ncRNA_gene:pseudogene"
infilesT = grep("transcripts.txt", dir(), v=T)
transcripts = readTranscriptHostAll(infilesT, start_text = start_text,target = target,   filter = filter, combined_depth_thresh = 1)
                                  
                                   
info = attr(transcripts,'info')

transcripts = cbind(gfft[match(transcripts$geneID, gfft$ID),],transcripts)
ord = order(transcripts$countT, decreasing=T)
head(transcripts[ord,])
##find DE genes
DE1 = DEgenes(transcripts, control_inds, infected_inds,edgeR = F, reorder=F);
pos1M = apply(cbind(as.character(DE1$chrs), as.character(binsize*round(DE1$start/binsize))),1,paste,collapse=".")
#if(mergeByPosAndGene){
#  pos1M = apply(cbind(as.character(DE1$chrs), DE1$geneID, DE1$rightGene,  as.character(binsize*round(DE1$start/binsize))),1,paste,collapse=".")
#}
DE2 = cbind(DE1,pos1M)

desciGT = .getDescEnrich(DE2,mart,thresh = 1e-5, go_thresh=1e-2,nme="pvals1",nme2="pos1M")
desciLT = .getDescEnrich(DE2,mart,thresh = 1e-5, go_thresh=1e-2,nme="pvals2",nme2="pos1M")




pdf("qq.pdf")
.qqplot(DE1$pvals1, min.p= 1e-200,main="both")
.qqplot(DE1$pvals2, min.p= 1e-200,main="infected less",add=T)
.vis(DE1,i=1,min.p=1e-50)
 .vis(DE1,i=2,min.p=1e-50)

dev.off()

#this calculates pvalues for base-level error rates
depth=.readH5All(transcripts, chroms = attr(transcripts,"chroms"),thresh = 100)

depth1 = .appendGeneNamesToDepth(depth, transcripts, sort=NULL)
#depth2 = .appendGeneNamesToDepth(depth, transcripts, sort="pv2")

sigChr1 = findSigChrom(depth1, thresh=1e-3, go_thresh=1e-3, nme="pv1", nme2="gene_names")
sigChr2 = findSigChrom(depth1, thresh=1e-3, go_thresh=1e-3, nme="pv2", nme2="gene_names")

##this is visualisation
  pdf("error_associations.pdf")
  pv_inds = grep("pv", names(depth))
  .qqplot(depth$pv1,min.p=1e-100)
  .qqplot(depth$pv2,min.p=1e-100)
   hist(depth[depth$pv1<1e-3,]$base) ## this is up in controls
   hist(depth[depth$pv2<1e-3,]$base) ## this is up in cases
  .vis(depth,i=1,min.p=1e-20)
  .vis(depth,i=2,min.p=1e-20)
  dev.off()

 # di_2 =   depth$pv2<1e-5 | depth$pv1<1e-5
  
  
chr_inds=attr(depth,"chr_inds")
chroms = attr(depth,"chroms")
chroms1 = names(chroms)
mi = match(chr_inds, chroms)
mi1 = unlist(lapply(mi, function(x) chroms1[x]))
pos100k = apply(cbind(as.character(mi1),as.character(binsize*floor(depth1$pos/binsize))),1,paste,collapse=".")

if(mergeByPosAndGene){
  pos100k = apply(cbind(as.character(mi1),as.character(depth1$gene_names), as.character(binsize*floor(depth1$pos/binsize))),1,paste,collapse=".")
  
}
depth2 = cbind(depth1, pos100k)
sigChr1 = findSigChrom(depth2, thresh=1e-5, go_thresh=1e-3, nme="pv1", nme2="pos100k")
DE_sig1 = DE2[DE2$pos1M %in% sigChr1$chrs,]
sigChr2 = findSigChrom(depth2, thresh=1e-5, go_thresh=1e-3, nme="pv2", nme2="pos100k")
DE_sig2 = DE2[DE2$pos1M %in% sigChr2$chrs,]

#depth2[depth2$pos100k %in% sigChr2$chrs,]

if(!.is.null(mart)){
  DE_sig1 = getDescr(DE_sig1, mart,thresh = 1e-3, prefix=prefix)
  
  DE_sig2 = getDescr(DE_sig2, mart,thresh = 1e-3, prefix=prefix)
}

#CHECK DISTRIBUTION (OPTIONAL)
.qqplot(DE1$FDR1)
#.qqplot(DE1$FDR2)




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
