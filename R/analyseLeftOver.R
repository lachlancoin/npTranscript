 #samtools view -F 2048 dRNA_Vero_DavidsonEtAl_guppy351_leftover.HumanRNA_HumanrDNA_SCOV2.bam | grep '>'  | cut -f 3,4 | grep -v 'MT'  | sort -k 3 -k 4,4g | less -S

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
data_src = c("~/github/npTranscript/data/SARS-Cov2","~/github/npTranscript/data/229E_CoV" )
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
src = c("~/github/npTranscript/R", "../../R")
#data_src =  # c(".","..","~/github/npTranscript/data/SARS-Cov2" )
print("#PRELIMINARIES ....")  
source(.findFile(src, "diff_expr_functs.R"))    
source(.findFile(src, "transcript_functions.R"))

a = read.table("tmp.txt")

head(getlev(b),20)
 b[,2] = floor(a[,2]/10)*10
