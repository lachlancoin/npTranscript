library(rhdf5)
file="chrom_out1"

h5f = paste(file, "h5",sep=".")         
tf = paste(file, "h5.transcripts.mod.txt.gz",sep=".")
h5ls(h5f)
transcripts = read.table(tf, header=F)[,1]
barcodes = h5read("chrom_out1.h5", "0/barcodes")
indices = h5read("chrom_out1.h5", "0/indices")
counts = h5read("chrom_out1.h5", "0/counts")

matr =data.frame(array(0,dim = c(length(transcripts), length(barcodes)+1),dimnames = list(NULL, c("GENE",barcodes))))

for(i in 1:length(barcodes)){
  inds = indices[,i]+1
 cnts = counts[,i]
  inds1 = inds>0
  matr[inds[inds1],i+1] = cnts[inds1]


}
matr[,1] = transcripts
write.csv( matr,  paste0(file, "_tex.txt"), quote=F, row.names=F)
sars_ind = grep("MN908947", transcripts)

meta = array(dim = c(length(barcodes)+1, 4), dimnames = list(NULL, c("NAME" ,   "Cohort", "Cell_State", "Patient")))

meta[-1,1] = barcodes
meta[-1,4] = "Adult"
meta[-1,3] = "Alpha"
meta[which(as.numeric(matr[sars_ind,])>0),2]="Infected"
meta[which(matr[sars_ind,]==0),2]="Uninfected"
meta[1,] = rep("group" , ncol(meta))
meta[1,1]="TYPE"
write.table( meta,  paste0(file, "_meta.txt"), quote=F, row.names=F,sep="\t")


