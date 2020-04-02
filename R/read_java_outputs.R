
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
infilesT1 = grep("transcripts1.txt", dir(), v=T)
infilesE = grep("exons.txt", dir(), v=T)
infilesBr = grep("breakpoints.txt.gz.[0-9]", dir(), v=T)
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
gen_inds = c( grep("upstream", names(reads)),grep("downstream", names(reads)))
comb_l = as.factor(apply(reads[,gen_inds],1,function(x) paste(as.character(x), collapse=".")))
reads = cbind(reads,comb_l)
reads$source =  factor(reads$source, levels = 0:(length(type_nme)-1), labels=type_nme)

reads$type = as.factor(reads$type)
reads$upstream = as.factor(reads$upstream)
reads$downstream = as.factor(reads$downstream)
#write.table(reads, "0readToCluster.mod.txt",sep="\t", row.names=F, col.names=T, quote=F)
reads_leader = reads[ reads$type=="5_3",]
reads0  = reads[reads$source=="Cell" & reads$type=="5_3" ,]
reads1  = reads[reads$source=="Virion" & reads$type=="5_3" ,]


   # geom_bar(stat="identity", position=position_dodge()) +
   # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
    #              position=position_dodge(.9))


#ggplot(reads_leader, aes(x=source, y=errorRatio)) + geom_violin()+ggtitle("Error vs type")
 outfile0 = paste(resdir, "/transcript_eror.pdf", sep="");
ggp<-ggplot(reads[reads$upstream=="LEADER",], aes(x=downstream, y=errorRatio, fill = source)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+ggtitle("Error vs transcript")
  try(ggsave(outfile0, plot=ggp, width = 50, height = 30, units = "cm"))



reads_ml = plotLengthHist(reads,t1, seqlen, type_nme)
  outfile1 = paste(resdir, "/length_hist.pdf", sep="");
  try(ggsave(outfile1, plot=reads_ml, width = 30, height = 30, units = "cm"))



transcripts_all = .readTranscripts(infilesT, seqlen, nmes)
names(transcripts_all)



transcript_counts = lapply(transcripts_all, function(transcripts)  apply(transcripts[,grep('count[0-9]', names(transcripts)), drop=F], 2,sum))
total_reads  = apply(data.frame(transcript_counts),1,sum)
max_h = unlist(lapply(transcripts_all, function(transcripts)  max(transcripts[,grep('count[0-9]', names(transcripts)), drop=F])))




if(dim(transcripts_all[[1]])[1]>0){
  maxpos = max(transcripts_all[[1]]$end)
  minpos = min(transcripts_all[[1]]$start)
}else{
  maxpos = length(fasta[[1]])
  minpos = 1
total_reads = c(1,1)
}
