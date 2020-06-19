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
if(length(args)==0){
  args = c("control=control","infected=infected", "betabinom","none")
}

  control_names = unlist(strsplit(args[1],':'))
  infected_names = unlist(strsplit(args[2],':'))
  type_names = c("control","infected")
  if(length(control_names[1])==1) type_names[1] = control_names[1]
  if(length(infected_names[1])==1) type_names[2] = infected_names[1]
  
  type_names = c(control_names[1], infected_names[1])

	control_names = unlist(strsplit(control_names1[2],':'))
	infected_names =unlist(strsplit(infected_names1[2],':'))
	analysis=args[3]
	
	exclude_nme = if(length(args)<4) "do_not_include"  else args[4]
	edgeR = FALSE;
	if(analysis=="edgeR") edgeR = T

library(VGAM)
library(edgeR)
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
             end="numeric", ID="character", isoforms="numeric" ,type_nme="character", countTotal="numeric")


##READ TRANSCRIPT DATA



infilesT = grep("transcripts.txt", dir(), v=T)
transcriptsl = readTranscriptHostAll(infilesT, start_text = start_text,target = target,   filter = filter, 
                                    combined_depth_thresh = 1)
attributes = attributes(transcriptsl)
filenames = attr(transcriptsl,"info")
control_names = unlist(lapply(control_names, grep, filenames, v=T));#  grep(control_names,filenames,v=T)
infected_names = unlist(lapply(infected_names, grep, filenames, v=T))
control_names = grep(exclude_nme, control_names, inv=T,v=T)
infected_names = grep(exclude_nme, infected_names, inv=T,v=T)


print(paste(type_names[1],paste(control_names,collapse=" ")))
print(paste(type_names[2] , paste(infected_names, collapse=" ")))

transcriptsl = lapply(transcriptsl, .processTranscripts)
filtered = .filter(transcriptsl)

pdf(paste(resdir, "/qq.pdf",sep=""))
.process(filtered$keep,attributes, resdir, control_names, infected_names, type_names= type_names, outp= "results.csv", type="keep")
.process(filtered$remove,attributes, resdir, control_names, infected_names, type_names = type_names, outp= "results_removed.csv", type="keep")
dev.off()

#this calculates pvalues for base-level error rates

print("####DEPTH ANALYSIS #### ")
##library(rhdf5)
#source(.findFile(src, "depth_association.R"))





