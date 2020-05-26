
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
#source(.findFile(src, "transcript_functions.R"))
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

#seqfile= "../Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa.gz"
#fasta = read.fasta(seqfile)
#seqlen = length(fasta[[1]])


target= list( chrom="character", 
             ORFs ="character",start = "numeric", 
             end="numeric", ID="character", isoforms="numeric" ,type_nme="character", countTotal="numeric")


##READ TRANSCRIPT DATA



          

#"gene:ncRNA_gene:pseudogene"
infilesT = grep("transcripts.txt", dir(), v=T)
transcripts = readTranscriptHostAll(infilesT, start_text = start_text,target = target,   filter = filter, combined_depth_thresh = 1)
names(transcripts)[grep("count[0-9]",names(transcripts))] = sub("_leftover", "" ,attr(transcripts,"info"))
geneID= as.character(unlist(lapply(strsplit(transcripts$ORFs,";"), function(v) v[1])))
rightGene = as.character(unlist(lapply(strsplit(transcripts$ORFs,";"), function(v) v[length(v)])))
transcripts = cbind(transcripts, geneID, rightGene)


getlev(transcripts$chr)      
head(getlev(transcripts[,names(transcripts) %in% c("chrs","start")]))                
                                   
info = attr(transcripts,'info')

gfft = read.table("annotation.csv.gz", sep="\t", header=F, fill=T, quote='\"')
gfft = gfft[match(transcripts$geneID, gfft[,1]),]
gfft[,1] = transcripts$ID
names(gfft)[1] = "ID"
names(gfft)[2] = "Name"

transcripts = cbind(gfft,transcripts)
ord = order(transcripts$countT, decreasing=T)
head(transcripts[ord,])

