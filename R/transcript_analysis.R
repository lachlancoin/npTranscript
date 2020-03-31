#should run this in  subdirectory with results from java program

###PRELIM SECTION
#library(dpylr)
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
      
      
sourcePath(src, "transcript_functions.R")
resdir = "results"
dir.create(resdir);


plotHist=T
sourcePath(src, "read_java_outputs.R")


RUN_ALL = TRUE
todo = 1:length(type_nme)
sourcePath(src, "breakpoint_analysis.R")

#HEATMAP = TRUE
COVERAGE = TRUE
sourcePath(src, "coverage_analysis.R")
