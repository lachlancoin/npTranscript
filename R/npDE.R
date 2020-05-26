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

args = commandArgs(trailingOnly=TRUE)
if(length(args)>0){
	control_names = unlist(strsplit(args[1],':'))
	infected_names = unlist(strsplit(args[2],':'))
	prefix = args[3]
}else{
	control_names = "korean_monkey_control"
	infected_names = "korean_monkey_infected"
	prefix = "ENSC"; #for vervet monkey
}




library(VGAM)
#library(biomaRt)
#library(edgeR)
library(stats)
library(rhdf5)
#read.gff(file, na.strings = c(".", "?"), GFF3 = TRUE)
#library(seqinr)

resdir = "results"
dir.create(resdir);

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
attributes = attributes(transcripts)
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


##find DE genes

DE1 = DEgenes(transcripts, control_names, infected_names,edgeR = F, reorder=F);
#pos1M = apply(cbind(as.character(DE1$chrs), as.character(binsize*round(DE1$start/binsize))),1,paste,collapse=".")
#if(mergeByPosAndGene){
#  pos1M = apply(cbind(as.character(DE1$chrs), DE1$geneID, DE1$rightGene,  as.character(binsize*round(DE1$start/binsize))),1,paste,collapse=".")
#}
#DE2 = cbind(DE1,pos1M)
DE1 = .transferAttributes(DE1, attributes)

write.table(DE1[attr(DE1,"order"),],file=paste(resdir,"results.csv",sep="") , quote=F, row.names=F, sep="\t", col.names=T)


pdf(paste(resdir, "/qq.pdf",sep=""))
.qqplot(DE1$pvals1, min.p= 1e-200,main="both")
.qqplot(DE1$pvals2, min.p= 1e-200,main="infected less",add=T)
.vis(DE1,i=1,min.p=1e-50)
 .vis(DE1,i=2,min.p=1e-50)

dev.off()

#this calculates pvalues for base-level error rates

print("####DEPTH ANALYSIS #### ")
#source(.findFile(src, "depth_association.R"))



###BIOMART ANALYSIS
dataset="csabaeus_gene_ensembl"
 mirror = "uswest"
print("####BIOMART ANALYISIS #### ")
#source(.findFile(src, "biomar_analysis.R"))

