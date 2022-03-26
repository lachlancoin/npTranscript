library(rhdf5)
#lib="~/R/lib/3.6.0"
#install.packages("BiocManager",lib=lib)
#library(BiocManager, lib=lib)
#BiocManager::install("rhdf5",lib=lib)


args = commandArgs(trailingOnly=TRUE)
if(length(args)>0){
  args = gsub("--","",args)
  argv = (lapply(args, function(x) strsplit(x,"=")[[1]][2]))
  names(argv) = unlist(lapply(args, function(x) strsplit(x,"=")[[1]][1]))
  options(argv)
}


file=getOption("file","chrom_out1")
infile=paste(file, "h5",sep=".")
transcript_file = getOption("transcript_file","../all_transcripts.txt")
print(transcript_file)
h5f = paste(file, "h5",sep=".")         
tf = paste(file, "h5.transcripts.mod.txt.gz",sep=".")


h5ls(h5f)
transcripts = read.table(tf, header=F)[,1]
transcripts1 = read.table(transcript_file)[,1]
mi = match(transcripts, transcripts1)
barcodes = h5read(infile, "0/barcodes")
indices = h5read(infile, "0/indices")
counts = h5read(infile, "0/counts")

matr =data.frame(array(0,dim = c(length(transcripts1), length(barcodes)+1),dimnames = list(NULL, c("GENE",barcodes))))

print(mi)
print(dim(matr))

for(i in 1:length(barcodes)){
  inds = indices[,i]+1
  cnts = counts[,i]
  inds1 = inds>0
  inds_ = inds[inds1]
#  print(i)
#  print(inds_)
  matr[mi[inds_],i+1] = cnts[inds1]
}
matr[,1] = transcripts1
#comma = grep(",", transcripts1)
#if(length(comma)>0){
#  matr = matr[-comma,]
#}

#sars_ind = grep("MN908947", matr[,1])
#sars = matr[sars_ind[1],]
#matr = rbind(sars,matr[-sars_ind[1],])


write.table( matr,  paste0(file, "_tex.txt"), quote=F, row.names=F, sep="\t")

#meta = array(dim = c(length(barcodes)+1, 4), dimnames = list(NULL, c("NAME" ,   "Cohort", "Cell_State", "Patient")))

#meta[-1,1] = barcodes
#meta[-1,4] = "Adult"
#meta[-1,3] = "Alpha"
#meta[which(as.numeric(matr[sars_ind,])>0),2]="Infected"
#meta[which(matr[sars_ind,]==0),2]="Uninfected"
#meta[1,] = rep("group" , ncol(meta))
#meta[1,1]="TYPE"
#write.table( meta,  paste0(file, "_meta.txt"), quote=F, row.names=F,sep="\t")


