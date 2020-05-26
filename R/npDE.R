install= FALSE
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
#library(VGAM)
#library(biomaRt)
#library(edgeR)
library(stats)
library(rhdf5)
#read.gff(file, na.strings = c(".", "?"), GFF3 = TRUE)
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
source(.findFile(src, "transcript_functions.R"))
source(.findFile(src, "diff_expr_functs.R"))


prefix = "ENSC"; #for vervet monkey




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



target= list( chrom="character", 
             ORFs ="character",start = "numeric", 
             end="numeric", ID="character", isoforms="numeric" ,type_nme="character", countTotal="numeric")


##READ TRANSCRIPT DATA



infilesT = grep("transcripts.txt", dir(), v=T)
transcripts = readTranscriptHostAll(infilesT, start_text = start_text,target = target,   filter = filter, combined_depth_thresh = 1)
names(transcripts)[grep("count[0-9]",names(transcripts))] = sub("_leftover", "" ,attr(transcripts,"info"))
geneID= as.character(unlist(lapply(strsplit(transcripts$ORFs,";"), function(v) v[1])))
rightGene = as.character(unlist(lapply(strsplit(transcripts$ORFs,";"), function(v) v[length(v)])))
transcripts = cbind(transcripts, geneID, rightGene)


head(getlev(transcripts$chr))      
#head(getlev(transcripts[,names(transcripts) %in% c("chrs","start")]))                
                                   
info = attr(transcripts,'info')

gfft = read.table("annotation.csv.gz", sep="\t", header=F, fill=T, quote='\"')
#gfft[,1] = gsub("transcript:", "", as.character(gfft[,1]))
gfft = gfft[match(transcripts$geneID, gfft[,1]),]
gfft[,1] = transcripts$ID
names(gfft)[1] = "ID"
names(gfft)[2] = "Name"

transcripts = cbind(gfft,transcripts)
ord = order(transcripts$countT, decreasing=T)
head(transcripts[ord,])



ord = order(transcripts$countT, decreasing=T)
head(transcripts[ord,])
##find DE genes
DE1 = DEgenes(transcripts, control_inds, infected_inds,edgeR = F, reorder=F);
pos1M = apply(cbind(as.character(DE1$chrs), as.character(binsize*round(DE1$start/binsize))),1,paste,collapse=".")
#if(mergeByPosAndGene){
#  pos1M = apply(cbind(as.character(DE1$chrs), DE1$geneID, DE1$rightGene,  as.character(binsize*round(DE1$start/binsize))),1,paste,collapse=".")
#}
DE2 = cbind(DE1,pos1M)




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

#CHECK DISTRIBUTION (OPTIONAL)
.qqplot(DE1$FDR1)
#.qqplot(DE1$FDR2)


###BIOMART ANALYSIS
dataset="csabaeus_gene_ensembl"
 mirror = "uswest"
print("####BIOMART ANALYISIS #### ")
#source(.findFile(src, "biomar_analysis.R"))

