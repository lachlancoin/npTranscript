
##READING INPUT FILES REQUIRED FOR LATER ANALYSIS (transcript  information)
t = readCoords("../Coordinates.csv")
subt_inds = t$gene!="none" & t$gene!="leader"
t1 = t[subt_inds,]
dimnames(t1)[[1]] = t[subt_inds,7]
fimo = read.table("../FIMO.csv", sep=",", head=T)

fastafile = paste("../" ,grep(".fasta" , dir("../"),v=T),sep="")
if(length(fastafile)!=1) stop("should just have one fasta file in parent directory")
fasta = read.fasta(fastafile)
seqlen = length(fasta[[1]])
readlen = seqlen - t1$Minimum
t1 = cbind(t1, readlen)

fastaseq = paste(fasta[[1]], collapse="")

infiles = grep("clusters.h5", dir(), v=T)
infilesT = grep("transcripts.txt", dir(), v=T)
infilesE = grep("exons.txt", dir(), v=T)
infilesBr = grep("breakpoints.txt.[0-9]", dir(), v=T)
infilesReads = grep("0readToCluster", dir(), v=T)

nmes = c("5_3", "5_no3","no5_3", "no5_no3")
mult = 1;
options("heatmap.2" = TRUE)

sumT_all = matrix(nrow = length(infiles), ncol = length(type_nme))
dimnames(sumT_all)[[2]] = type_nme
leader=tolower("ACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAAC")
leader_ind = gregexpr(pattern = leader, fastaseq)[[1]][1]
leader_ind = c(leader_ind, leader_ind + nchar(leader)-1)


#if(plotHist) pdf(paste(resdir,'hist.pdf',sep="/"));

 df = data.frame(min = t1$Minimum[-dim(t1)[1]],max = t1$Minimum[-1], gene = as.character(t1$gene[-1]))

reads=  data.frame(read.table(infilesReads, sep="\t",  head=T))
rightGene = rep(NA, dim(reads)[1])
for(i in 1:dim(df)[1]){
indsi = which(reads$breakEnd-5<= df[i,2]  & reads$breakEnd> df[i,1]+5)
rightGene[indsi] = rep(as.character(df$gene[i]), length(indsi))
}
reads = cbind(reads, rightGene)
reads$source = as.factor(reads$source)
write.table(reads, "0readToCluster.mod.txt",sep="\t", row.names=F, col.names=T, quote=F)
reads_leader = reads[ reads$startPos < 100 & reads$endPos >seqlen-100 ,]
reads0  = reads[reads$source==0 & reads$startPos < 100 & reads$endPos >seqlen-100 ,]
reads1  = reads[reads$source==1 & reads$startPos < 100 & reads$endPos >seqlen-100 ,]
pdf("output.pdf")
plot(reads$source, reads$errorRatio, main = "Error vs type")
plot(reads0$rightGene, reads0$errorRatio, main= type_nme[1])
plot(reads1$rightGene, reads1$errorRatio, main = type_nme[2])
dev.off()


reads_ml = plotLengthHist(reads,t1, seqlen, type_nme)
  outfile1 = paste(resdir, "/length_hist.pdf", sep="");
  try(ggsave(outfile1, plot=reads_ml, width = 30, height = 30, units = "cm"))



transcripts = read.table( infilesT,sep="\t", head=T)
o = order(transcripts$countTotal, decreasing=T)
transcripts = transcripts[o,]
err_ratio_inds = grep("error_ratio", names(transcripts))
transcripts[,err_ratio_inds] =apply(transcripts[,err_ratio_inds,drop=F], c(1,2), function(x) if(is.na(x)) -0.01 else x)

transcripts_all = list()
transcripts_all[[1]] = (transcripts[which(transcripts$start<100 & transcripts$end > seqlen -100),, drop=F])
transcripts_all[[2]] = (transcripts[which(transcripts$start<100 & transcripts$end <= seqlen -100),, drop=F])
transcripts_all[[3]] =( transcripts[which(transcripts$start>=100 & transcripts$end > seqlen -100),, drop=F])
transcripts_all[[4]] = (transcripts[which(transcripts$start>=100 & transcripts$end <= seqlen -100),, drop=F])
 

 #ggp2<-ggplot(transcripts_all[[1]],aes_string(x="error_ratio0",y="error_ratio1",fill =  "rightGene" , colour = "rightGene" ))+geom_point(size=.1)# + theme_bw(); #+ggtitle(title2)


  
  #exons = read.table( infilesE,sep="\t", head=F, skip=1)
  #names(exons) = c(read.table( infilesE,sep="\t", head=F, nrows=1,as.is=T), "type")
  #exons_all[[ik]] = exons
transcript_count_matrix = transcripts[,grep('count[0-9]', names(transcripts))]
total_reads  = apply(transcript_count_matrix,2,sum)
max_h = max(transcript_count_matrix)

maxpos = max(transcripts$end)
minpos = min(transcripts$start)
 dev.off()

if(dim(transcripts_all[[1]])[1]>0){
  maxpos = max(transcripts_all[[1]]$end)
  minpos = min(transcripts_all[[1]]$start)
}else{
  maxpos = length(fasta[[1]])
  minpos = 1
total_reads = c(1,1)
}
