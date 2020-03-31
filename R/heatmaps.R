library (gplots)
#coref = read.table("data/bin8_sg_reads.fq.bamcoref.txt.gz", sep=",", row.names = 1)
 
coref = read.table( "data/allreads_fq.bamcoref.txt.gz",sep=",",  row.names=1)
dimnames(coref)[[2]] = dimnames(coref)[[1]]
inds = 1:dim(coref)[1] 
 corefM = apply(as.matrix(coref), c(1,2), function(x) if(x==0) NA else log10(x))

 heatmap.2(corefM[1:10,1:10], dendrogram = "none", trace ="none", Rowv = NULL, Colv = NULL, scale="none")
 
 heatmap(corefM[rev(inds), inds], Rowv = NA, Colv = NA, scale="none")

 
 

 dimnames(coref)[[2]] = dimnames(coref)[[1]]
 coref = apply(coref, c(1,2), function(x) if(x==0) NA else log10(x))
 corefM = as.matrix(coref)
 inds1 = inds[1:500]
 heatmap(corefM[rev(inds), inds], Rowv = NA, Colv = NA, scale="none")
 
 
 
 setwd("C:/Users/LCOIN/bitbucket/rapids_rna_seq/corona/data")
 infiles = grep("txt.gz",grep("allreads_fq.bamcoref.txt.", dir(), v=T), v=T, inv=T)
 for(i in 1:length(infiles)){
   infile = infiles[i]
   print(infile)
   outfile = paste(infile,"pdf", sep=".")
   coref = read.table( infile,sep=",",  row.names=1)
   inds = 1:dim(coref)[1]
   dimnames(coref)[[2]] = dimnames(coref)[[1]]
   corefL = apply(coref, c(1,2), function(x) if(x==0) NA else log10(x))
 
   pdf(outfile)
   heatmap.2(as.matrix(corefL), dendrogram = "none", trace ="none", Rowv = NULL, Colv = NULL, scale="none")
   
   #heatmap(corefM[rev(inds), inds], Rowv = NA, Colv = NA, scale="none")
   dev.off()
 }
 
 