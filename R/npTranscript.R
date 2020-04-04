#should run this in  subdirectory with results from java program
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(gplots)
library(seqinr)
library(rhdf5)
args = commandArgs(trailingOnly=TRUE)

type_nme = c( "Cell","Virion") #,"Cell2")

#SHOULD BE RUN IN data/ subdirectory
sourcePath<-function(path, files){
  for(i in 1:length(path)){
    if(file.exists(path[i])) {
      for(j in 1:length(files)){
        source(paste(path[i],files[j],sep="/"))
      }
      return(path[i])
      break
    }
  }
}

src = c("~/github/npTranscript/R" )
      
#PRELIMINARIES      
sourcePath(src, "transcript_functions.R")
resdir = "results"
dir.create(resdir);
t = readCoords("../Coordinates.csv")
subt_inds = t$gene!="none" & t$gene!="leader"
t1 = t[subt_inds,]
dimnames(t1)[[1]] = t[subt_inds,7]
fimo = read.table("../FIMO.csv", sep=",", head=T)
fastafile = paste("../" ,grep(".fasta" , dir("../"),v=T),sep="")
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


#infilesT1 = grep("transcripts1.txt", dir(), v=T)
#infilesE = grep("exons.txt", dir(), v=T)


###READ LEVEL ANALYSIS
infilesReads = grep("0readToCluster", dir(), v=T)
if(length(infilesReads)==1){
	sourcePath(src, "read_analysis.R")
}else{
	print("no reads file")
}
plotHist=T
sourcePath(src, "read_java_outputs.R")


##COVERAGE ANALYSIS
infilesT = grep("transcripts.txt", dir(), v=T)
infiles = grep("clusters.h5", dir(), v=T)
if(length(infilesT)==1 && length(infiles)==1){
	HEATMAP = TRUE
	COVERAGE = TRUE
	sourcePath(src, "coverage_analysis.R")
}else{
	print("no transcripts file")
}

infilesBr = grep("breakpoints.txt.gz.[0-9]", dir(), v=T)
if(length(infilesT)==1 && length(infiles)==1){
	RUN_ALL = TRUE
	todo = 1:length(type_nme)
	sourcePath(src, "breakpoint_analysis.R")
}else{
 	print("no break point files")
}



