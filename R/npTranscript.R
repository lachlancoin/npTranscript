#should run this in  subdirectory with results from java program
INSTALL = FALSE
if(INSTALL){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
	BiocManager::install("rhdf5")
	BiocManager::install("ggplot2")
	BiocManager::install("gridExtra")
	BiocManager::install("RColorBrewer")
	BiocManager::install("gplots")
	BiocManager::install("seqinr")

}
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(gplots)
library(seqinr)
library(rhdf5)
#args = commandArgs(trailingOnly=TRUE)

#print(args)
#print(length(args))
#if(length(args)==0) stop("need to specify the type nmes, e.g. Cell:Virion")

 
#type_nme=	c( "Cell","Virion") #,"Cell2")

#SHOULD BE RUN IN data/ subdirectory
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


#type_nme = strsplit(args[1], ":")[[1]]
infilesBr = grep("breakpoints.", dir(), v=T)
infilesReads = grep("0readToCluster", dir(), v=T)
infilesT = grep("transcripts.txt.gz$", dir(), v=T)
infiles = grep("clusters.h5", dir(), v=T)
print(infilesBr)
print(infilesReads)
print(infilesT)
print(infiles)
type_nme = NULL
if(length(infilesBr)>0){
type_nme = unlist(lapply(infilesBr, function(x) strsplit(x,"\\.")[[1]][2]))
if(length(infilesBr)!=length(type_nme)){
	print(type_nme)
	print(infilesBr)
       	stop("type nme length does not match number of breakpoint files")
}
print(type_nme)
}
src = c("~/github/npTranscript/R" )
data_src =   c(".","..","~/github/npTranscript/data/SARS-Cov2" )
print("#PRELIMINARIES ....")      
source(.findFile(src, "transcript_functions.R"))
resdir = "results"
dir.create(resdir);
t = readCoords(.findFile(data_src, "Coordinates.csv"))
subt_inds = t$gene!="none" & t$gene!="leader"
t1 = t[subt_inds,]
dimnames(t1)[[1]] = t[subt_inds,7]
fimo = read.table(.findFile(data_src,"FIMO.csv"), sep=",", head=T)
fastafile = .findFile(data_src,".fasta.gz$", exact=F)

if(length(fastafile)!=1) stop("should just have one fasta file in parent directory")
fasta = read.fasta(fastafile)
seqlen = length(fasta[[1]])
readlen = seqlen - t1$Minimum
t1 = cbind(t1, readlen)
# df = data.frame(min = t1$Minimum[-dim(t1)[1]],max = t1$Minimum[-1], gene = as.character(t1$gene[-1]))
fastaseq = paste(fasta[[1]], collapse="")
nmes = c("5_3", "5_no3","no5_3", "no5_no3")
mult = 1;
options("heatmap.2" = TRUE)
leader=tolower("ACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAAC")
leader_ind = gregexpr(pattern = leader, fastaseq)[[1]][1]
leader_ind = c(leader_ind, leader_ind + nchar(leader)-1)





transcripts_all = .readTranscripts(infilesT, seqlen, nmes)
names(transcripts_all)
if(is.null(type_nme)){
type_nme = attr(transcripts_all[[1]], "info")
}
count_df = grep('count[0-9]', names(transcripts_all[[1]]))

transcript_counts = lapply(transcripts_all, function(transcripts)  apply(transcripts[,count_df, drop=F], 2,sum))
total_reads  = apply(data.frame(transcript_counts),1,sum)
names(total_reads) = type_nme
max_h = unlist(lapply(transcripts_all, function(transcripts)  max(transcripts[,grep('count[0-9]', names(transcripts)), drop=F])))
if(dim(transcripts_all[[1]])[1]>0){
  maxpos = max(transcripts_all[[1]]$end)
  minpos = min(transcripts_all[[1]]$start)
}else{
  maxpos = length(fasta[[1]])
  minpos = 1
	total_reads = c(1,1)
}
ml1 =.plotGeneExpr(transcripts_all)




print("###READ LEVEL ANALYSIS")

if(length(infilesReads)==1){
	source(.findFile(src, "read_analysis.R"))
}else{
	print("no reads file")
}
#plotHist=T
#sourcePath(src, "read_java_outputs.R")



print('##COVERAGE ANALYSIS')
if(length(infilesT)==1 && length(infiles)==1){


	HEATMAP = TRUE
	COVERAGE = TRUE
	source(.findFile(src, "coverage_analysis.R"))
}else{
	print("no break point files")
}

print('##BREAKPOINT ANALYSIS')
if(length(infilesBr)>=1 && length(infiles)>=1){

	RUN_ALL = TRUE
	todo = 1:length(type_nme)
	source(.findFile(src, "breakpoint_analysis.R"))
}else{
 	print("no break point files")
}



