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
	BiocManager::install("binom")
}
library(binom)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(gplots)
library(seqinr)
library(rhdf5)
library(VGAM)
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0){
data_src = args[1]  ## location of fasta file and Coordinates file
}else{
data_src = c("C:/Users/LCOIN/github/npTranscript/data/SARS-Cov2/VIC01" ,"~/github/npTranscript/data/SARS-Cov2/VIC01","~/github/npTranscript/data/229E_CoV" )
}


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
		return(paste(path[i], files,sep="/"))
   } 
  }
}


#type_nme = strsplit(args[1], ":")[[1]]
infilesBr = grep("breakpoints.", dir(), v=T)
infilesReads = grep("readToCluster", dir(), v=T)
infilesAnnot = grep("annot.txt.gz$", dir(), v=T)
infilesT = grep("transcripts.txt.gz$", dir(), v=T)
infiles = grep("clusters.h5", dir(), v=T)
infilesAltT = grep("isoforms.h5", dir(), v=T)
print(infilesBr)
print(infilesReads)
print(infilesT)
print(infiles)
print(infilesAltT)
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
src = c("~/github/npTranscript/R", "C:/Users/LCOIN/github/npTranscript/R")
#data_src =  # c(".","..","~/github/npTranscript/data/SARS-Cov2" )
print("#PRELIMINARIES ....")  
source(.findFile(src, "diff_expr_functs.R"))    
source(.findFile(src, "transcript_functions.R"))
resdir = "results"
dir.create(resdir);
t = readCoords(.findFile(data_src, "Coordinates.csv"))
subt_inds = t$gene!="none" & t$gene!="leader"
t1 = t[subt_inds,]
dimnames(t1)[[1]] = t[subt_inds,which(names(t)=='gene')]
fimo_file = .findFile(data_src,"fimo.tsv")

if(!is.null(fimo_file)){
	fimo = read.table(fimo_file, sep="\t", head=T)
}else{
	fimo = NULL
}
fastafile = .findFile(data_src,".fasta.gz$", exact=F)
fastafile = grep('leader',fastafile,v=T,inv=T)

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





transcripts = .readTranscripts(infilesT)
attributes = attributes(transcripts)
info = attr(transcripts, "info" )
#count_names =  paste("count",info, sep="_")
#names(transcripts)[grep("count[0-9]", names(transcripts))] =count_names
#names(transcripts)[grep("errors[0-9]", names(transcripts))] = paste("errors",info, sep="_")
#names(transcripts)[grep("error_ratio[0-9]", names(transcripts))] = paste("error_ratio",info, sep="_")

#.processDE1(transcripts,1,2,resdir)
#.processDE1(transcripts,1,3,resdir)
#.processDE1(transcripts,1,4,resdir)





transcripts_all = .splitTranscripts(transcripts, seqlen, nmes, splice=F)
transcripts_all_splice = .splitTranscripts(transcripts, seqlen, nmes, splice=T)

type_nme = attr(transcripts_all[[1]], "info")

#if(length(infilesAltT)>0){
#isoforms = readIsoformH5(infilesAltT[[1]],  transcripts_all[[1]])
#}
count_df = grep('count[0-9]', names(transcripts_all[[1]]))
transcript_counts = lapply(transcripts_all, function(transcripts)  apply(transcripts[,count_df, drop=F], 2,sum))
total_reads  = apply(data.frame(transcript_counts),1,sum)
names(total_reads) = type_nme
#max_h = unlist(lapply(transcripts_all, function(transcripts)  max(transcripts[,grep('count[0-9]', names(transcripts)), drop=F])))
#commented line above due to error (see github issue #3). max_h is created but not referred to.

maxmin_pos = unlist(lapply(transcripts_all, function(t) c(min(t$start),max(t$end))))
minpos = min(maxmin_pos)
maxpos = max(maxmin_pos)


tocompare = .getCompareVec(type_nme)
ml1 = lapply(.getCompareVec(type_nme), .plotGeneExpr, transcripts_all, todo = 1:length(transcripts_all))
ml1_splice = lapply(tocompare, .plotGeneExpr, transcripts_all_splice, todo = 1:length(transcripts_all_splice))
names(ml1) = unlist(lapply(tocompare, function(x)  paste(type_nme[x], collapse=".")))
names(ml1_splice) = unlist(lapply(tocompare, function(x)  paste(type_nme[x], collapse=".")))

for(i in 1:length(ml1)){
	outfile1 = paste(resdir, "/gene_expr_5_3.", names(ml1)[i],".pdf", sep="");
	try(ggsave(outfile1, plot=ml1[[i]], width = 30, height = 30, units = "cm"))
}
for(i in 1:length(ml1_splice)){
	outfile1 = paste(resdir, "/gene_expr_splice.", names(ml1_splice)[i],".pdf", sep="");
	try(ggsave(outfile1, plot=ml1_splice[[i]], width = 30, height = 30, units = "cm"))
}

if(length(infilesAnnot)>0){
  annot0 = .readAnnotFile(.findFile(data_src, "0.annot.tsv"),plot=T, type_nme=c("Cell","Virion"), showEB=T,conf.level=0.95)
 
	annots = .readAnnotFile(infilesAnnot,plot=T, type_nme=type_nme, annot0 = annot0,conf.level=0.95,showEB=F)
	double_inds = unlist(lapply(annots$annot[1,], typeof))=="double"
	annots$annot[,double_inds] = apply(annots$annot[,double_inds,drop=F],c(1,2), function(x)  gsub(' ','',sprintf("%5.3g",x)))
write.table(annots$annot,file=paste(resdir,"cellular_proportions.txt",sep="/"),sep=",",quote=F,col.names=T, row.names=F)
  outfile1 = paste(resdir, "/splice_vs_unspliced.pdf", sep="");
  try(ggsave(outfile1, plot=annots$ggp, width = 30, height = 30, units = "cm"))
}
print("###READ LEVEL ANALYSIS")

#if(length(infilesReads)>0){
#	source(.findFile(src, "read_analysis.R"))
#}else{
#	print("no reads file")
#}
#plotHist=T
#sourcePath(src, "read_java_outputs.R")
print('##COVERAGE ANALYSIS')
if(length(infilesT)==1 && length(infiles)==1){
	HEATMAP = TRUE
	COVERAGE = TRUE
	transcripts_all1 = transcripts_all
	source(.findFile(src, "coverage_analysis.R"))
	transcripts_all1 = transcripts_all_splice
	source(.findFile(src, "coverage_analysis.R"))
}else{
	print("no transcript files")
}

print('##BREAKPOINT ANALYSIS')
if(length(infilesBr)>=1 && length(infiles)>=1){

	RUN_ALL = TRUE
	todo = 1:length(type_nme)
	source(.findFile(src, "breakpoint_analysis.R"))
}else{
 	print("no break point files")
}

