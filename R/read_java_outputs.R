
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

 

reads=  data.frame(read.table(infilesReads, sep="\t",  head=T))
reads_ml = plotLengthHist(reads,t1, seqlen)
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
 

 ggp2<-ggplot(transcripts_all[[1]],aes_string(x="error_ratio0",y="error_ratio1",fill =  "rightGene" , colour = "rightGene" ))+geom_point(size=.1)# + theme_bw(); #+ggtitle(title2)

if(plotHist){
  for(i in 1:2){
  hist(rat_m[,i],main = paste(nmes[ik],type_nme[i]),br=20)
  }
}
  sumT = apply(transcripts[,grep("count[0-9]", names(transcripts)), drop=F],2,sum)
  sumT_all[ik,1:length(sumT)] = sumT
 
  transcripts_all[[ik]] = transcripts
  
  exons = read.table( infilesE,sep="\t", head=F, skip=1)
  names(exons) = c(read.table( infilesE,sep="\t", head=F, nrows=1,as.is=T), "type")
  #exons_all[[ik]] = exons
}
total_reads  = apply(sumT_all,2,sum, na.rm=T)
for(i in 1:(dim(sumT_all)[2])){
  sumT_all[,i] = sumT_all[,i]/total_reads[i]
  
}
if(plotHist)  dev.off()

if(dim(transcripts_all[[1]])[1]>0){
  maxpos = max(transcripts_all[[1]]$end)
  minpos = min(transcripts_all[[1]]$start)
}else{
  maxpos = length(fasta[[1]])
  minpos = 1
total_reads = c(1,1)
}
