install= FALSE
#ggforce for ggsina
if(install){
  install.packages("BiocManager")
  BiocManager::install("VGAM")
 BiocManager::install("ggplot2")
#  BiocManager::install("gplots")
  #BiocManager::install("jsonlite")
  #BiocManager::install("gridExtra")
  #BiocManager::install("GGally")
 # BiocManager::install("binom")
 BiocManager::install("biomaRt")
 BiocManager::install("edgeR")
 BiocManager::install("writexl")
 #BiocManager::install("gridExtra")
 BiocManager::install("ggrepel")
}

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  args = c("control","infected", "betabinom","none")
#  args[1] = "cDNA"
#  args[2] = "RNA"
 # args[4] = "GE"
}
prefix = "ENS"

  control_names = unlist(strsplit(args[1],':'))
  infected_names = unlist(strsplit(args[2],':'))
  type_names = c("control","infected")
  if(length(control_names[1])==1) type_names[1] = control_names[1]
  if(length(infected_names[1])==1) type_names[2] = infected_names[1]
  
  type_names = c(control_names[1], infected_names[1])

	analysis=args[3]
	
	exclude_nme = if(length(args)<4) "do_not_include"  else strsplit(args[4],":")[[1]]
	edgeR = FALSE;
	if(analysis=="edgeR"){
		library(edgeR)
		 edgeR = T
	}else{
	library(VGAM)
	}
library(stats)
	
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
             end="numeric", ID="character" ,type_nme="character", countTotal="numeric")


##READ TRANSCRIPT DATA



infilesT = grep("transcripts.txt", dir(), v=T)
transcriptsl = readTranscriptHostAll(infilesT, start_text = start_text,target = target,   filter = filter, 
                                    combined_depth_thresh = 100)

attributes = attributes(transcriptsl)
filenames = attr(transcriptsl,"info")
control_names = unlist(lapply(control_names, grep, filenames, v=T));#  grep(control_names,filenames,v=T)
infected_names = unlist(lapply(infected_names, grep, filenames, v=T))
for(i in 1:length(exclude_nme)){
control_names = grep(exclude_nme[i], control_names, inv=T,v=T)
infected_names = grep(exclude_nme[i], infected_names, inv=T,v=T)
}

print(paste(type_names[1],paste(control_names,collapse=" ")))
print(paste(type_names[2] , paste(infected_names, collapse=" ")))
if(length(control_names)!=length(infected_names)) error(" lengths different")
transcriptsl = lapply(transcriptsl, .processTranscripts)
if(length(transcriptsl)==1){
transcriptsl[[1]]$ID = paste("ID",1:length(transcriptsl[[1]]$ID),sep="")
}
filtered = .filter(transcriptsl, prefix)
keep = filtered$keep[unlist(lapply(filtered$keep,function(x) dim(x)[[1]]))>0]
transcripts_keep = .process(keep, control_names, infected_names)
remove = filtered$remove[unlist(lapply(filtered$remove,function(x) dim(x)[[1]]))>0]
transcripts_removed = .process(remove, control_names, infected_names)
transcripts_removed$ORFs =  
  apply(cbind(transcripts_removed$chrs, transcripts_removed$start, transcripts_removed$end),1,paste,collapse="_")
pdf(paste(resdir, "/qq.pdf",sep=""))

res_keep = .processDE(transcripts_keep,attributes, resdir, control_names, infected_names, type_names = type_names, outp= "results_keep.csv", type="keep",plot=T)
res_remove = .processDE(transcripts_removed,attributes, resdir, control_names, infected_names, type_names = type_names, outp= "results_removed.csv", type="keep",plot=T)
dev.off()
attr(res_keep$DE1, 'nme') = 'Known genes'
attr(res_remove$DE1, 'nme') = 'Novel genes'

DE_list = list("genes"=res_keep$DE1, 'novel'=res_remove$DE1)


library(writexl)
write_xlsx(DE_list, paste(resdir, "DE.xlsx",sep="/"))
#.volcano(res_keep$DE1, pthresh = 1e-3, prefix="keep")

library(ggplot2)
library(ggrepel)
volcanos = lapply(DE_list, .volcano, pthresh = 1e-5)
pdf(paste(resdir, "/qq.pdf",sep=""))
lapply(volcanos, function(x) print(x))
dev.off()
.head1<-function(transcripts, nme="order1",n=10){
  head(transcripts[attr(transcripts,nme),],n)
}
.head1(res_keep$DE1,"order",10)
.head1(res_remove$DE1,"order",10)

#findSigChrom(res_keep$DE1, thresh = 1e-10, go_thresh = 1e-2,nme="p.adj1", nme2="chrs")

#findSigChrom(res_keep$DE1, thresh = 1e-10, go_thresh = 1e-2,nme="p.adj2", nme2="chrs")

#this calculates pvalues for base-level error rates

print("####DEPTH ANALYSIS #### ")
##
#source(.findFile(src, "depth_association.R"))



####
if(FALSE){
transcripts =.combineTranscripts(transcriptsl, attributes)
ORFs = strsplit(transcripts$ORFs,";")
lens = unlist(lapply(ORFs,length))
mtch = unlist(lapply(ORFs, function(x) x[1]==x[2]))
trans1 = transcripts[mtch==FALSE,]
which(trans1$countTotal==max(trans1[trans1$chrs!="chr11",]$countTotal))

which(transcripts_removed$countTotal==max(transcripts_removed$countTotal))
ord = order(transcripts_removed$countTotal, decreasing=T)
tr = transcripts_removed[ord,]
head(cbind(tr$end - tr$start, tr$countTotal),30)
}
