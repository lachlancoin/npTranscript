
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
getlev(transcripts$chr)      
head(getlev(transcripts[,names(transcripts) %in% c("chrs","start")]))                
                                   
info = attr(transcripts,'info')

gfft = gfft[match(transcripts$geneID, gfft$ID),]
gfft[,1] = transcripts$ID

transcripts = cbind(gfft,transcripts)
ord = order(transcripts$countT, decreasing=T)
head(transcripts[ord,])

