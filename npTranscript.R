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
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0){
data_src = args[1]  ## location of fasta file and Coordinates file
}else{
data_src = c("~/github/npTranscript/data/SARS-Cov2/VIC01","~/github/npTranscript/data/229E_CoV","../data/SARS-Cov2/VIC01" )
}


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
#	if(length(files)==1){
		return(paste(path[i], files,sep="/"))
#	}
   } 
  }
}


#type_nme = strsplit(args[1], ":")[[1]]
infilesBr = grep("breakpoints.", dir(), v=T)
infilesReads = grep("readToCluster", dir(), v=T)
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
src = c("~/github/npTranscript/R", "../../R","./")
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

primerFile = .findFile(data_src[3],"nCoV-2019.tsv")
if(length(primerFiles)>0){
primers = read.table(primerFile, head=T, sep="\t")
indsp = c(grep("_1_LEFT",primers$name),grep("_93_RIGHT",primers$name))
indsp = c(indsp,grep('_95_', primers$name))
primers1 = primers[indsp,]
primers1Left = primers1[grep("LEFT",primers1$name),]
primers1Right = primers1[grep("RIGHT",primers1$name),]

fastastr = toupper(paste(fasta[[1]],collapse=""))
fastastr1 = (toupper(paste(rev(comp(fasta[[1]])),collapse="")))
coordsL = unlist(lapply(primers1Left$seq,regexpr,fastastr))
coordsR = length(fasta[[1]]) - unlist(lapply(primers1Right$seq,regexpr,fastastr1))

primers1Left = primers[grep("LEFT",primers$name),]
primers1Right = primers[grep("RIGHT",primers$name),]
coordsL = unlist(lapply(primers1Left$seq,regexpr,fastastr))
coordsR = length(fasta[[1]]) - unlist(lapply(primers1Right$seq,regexpr,fastastr1))
joint = cbind(coordsL, coordsR)
joint = cbind(joint, apply(joint,1,diff))
 dimnames(joint)[[1]] = primers1Left$name
joint[ which(joint[,2]>t$Minimum[10]),]

}

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
transcripts_all = .splitTranscripts(transcripts, seqlen, nmes, splice=T)

names(transcripts_all)
if(!is.null(attr(transcripts_all[[1]], "info"))){
type_nme = attr(transcripts_all[[1]], "info")
}

#if(length(infilesAltT)>0){
#isoforms = readIsoformH5(infilesAltT[[1]],  transcripts_all[[1]])
#}
count_df = grep('count[0-9]', names(transcripts_all[[1]]))
transcript_counts = lapply(transcripts_all, function(transcripts)  apply(transcripts[,count_df, drop=F], 2,sum))
total_reads  = apply(data.frame(transcript_counts),1,sum)ls
names(total_reads) = type_nme
#max_h = unlist(lapply(transcripts_all, function(transcripts)  max(transcripts[,grep('count[0-9]', names(transcripts)), drop=F])))
#commented line above due to error (see github issue #3). max_h is created but not referred to.
if(dim(transcripts_all[[1]])[1]>0){
  maxpos = max(transcripts_all[[1]]$end)
  minpos = min(transcripts_all[[1]]$start)
}else{
  maxpos = length(fasta[[1]])
  minpos = 1
	total_reads = c(1,1)
}
if(length(type_nme)==1){
	tocompare = list(c(1,1))
}else if(length(type_nme)==4){
	tocompare = list(c(1,3), c(1,2),c(3,4))
}else{
	tocompare = list();
	for(i in 2:length(type_nme)){
		tocompare[[i-1]] = c(1,i)
	}
}

print(tocompare)
print(transcripts_all)
ml1 = lapply(tocompare, .plotGeneExpr, transcripts_all, todo = 1:length(transcripts_all))
names(ml1) = unlist(lapply(tocompare, function(x)  paste(type_nme[x], collapse=".")))

for(i in 1:length(ml1)){
	outfile1 = paste(resdir, "/gene_expr.", names(ml1)[i],".pdf", sep="");
	try(ggsave(outfile1, plot=ml1[[i]], width = 30, height = 30, units = "cm"))
}

print("###READ LEVEL ANALYSIS")

if(length(infilesReads)>0){
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

# a = read.table("summary1.txt", sep=";")
# b = apply(t(data.frame(strsplit(as.character(a[,2]),','))),c(1,2), as.numeric)
# a =  apply(t(data.frame(strsplit(as.character(a[,1]),','))),c(1,2), as.numeric)
 #c = cbind(a[,3:4], b[,3:4])
#overl = apply(c,1,min(v[4] - v[1], v[2] -v[3])) 
.overlap<-function(v){
  
}
