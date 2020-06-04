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
if(length(args)>0){
	control_names = args[1] #unlist(strsplit(args[1],':'))
	infected_names =args[2] # unlist(strsplit(args[2],':'))
	prefix = args[3]
}else{
	control_names = "control"
	infected_names = "infected"
	prefix = "ENSC"; #for vervet monkey
}

library(VGAM)
#library(biomaRt)
#library(edgeR)
library(stats)
#library(rhdf5)
#read.gff(file, na.strings = c(".", "?"), GFF3 = TRUE)
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
#
filenames = attr(transcriptsl,"info")
#names(transcripts)[grep("count[0-9]",names(transcripts))] = filenames
control_names = grep(control_names,filenames,v=T)
infected_names = grep(infected_names,filenames,v=T)
print(paste("control",paste(control_names,collapse=" ")))
print(paste("infected" , paste(infected_names, collapse=" ")))

transcriptsl = lapply(transcriptsl, .processTranscripts)
transcriptsl = lapply(transcriptsl, .mergeRows,sum_names= c(control_names, infected_names), colid="geneID" )
transcripts=.combineTranscripts(transcriptsl)
attributes = attributes(transcripts)
info = attr(transcripts,'info')


transcripts=.addAnnotation("annotation.csv.gz", transcripts, colid="geneID", nmes = c("ID" , "Name" , "Description","biotype"))

ord = order(transcripts$countT, decreasing=T)
head(transcripts[ord,])


##find DE genes

DE1 = DEgenes(transcripts, control_names, infected_names,edgeR = F);

DE1 = .transferAttributes(DE1, attributes)

head(DE1[attr(DE1,"order"),])


.write(DE1, resdir)
  

pdf(paste(resdir, "/qq.pdf",sep=""))
.qqplot(DE1$pvals, min.p= 1e-200,main="both")
.qqplot(DE1$pvals1, min.p= 1e-200,main="infected_more")
.qqplot(DE1$pvals2, min.p= 1e-200,main="infected less",add=T)
#.vis(DE1,i=1,min.p=1e-50)
# .vis(DE1,i=2,min.p=1e-50)

dev.off()

#this calculates pvalues for base-level error rates

print("####DEPTH ANALYSIS #### ")
#source(.findFile(src, "depth_association.R"))



###BIOMART ANALYSIS
dataset="csabaeus_gene_ensembl"
 mirror = "uswest"
print("####BIOMART ANALYISIS #### ")
#source(.findFile(src, "biomar_analysis.R"))

