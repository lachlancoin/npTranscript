
reads=  data.frame(read.table(infilesReads, sep="\t",  head=T))
gen_inds = c( grep("upstream", names(reads)),grep("downstream", names(reads)))
comb_l = as.factor(apply(reads[,gen_inds],1,function(x) paste(as.character(x), collapse=".")))
reads = cbind(reads,comb_l)
reads = .appendGenePosition(reads, t1)
reads$source =  factor(reads$source, levels = 0:(length(type_nme)-1), labels=type_nme)

reads$type = as.factor(reads$type)
reads$upstream = as.factor(reads$upstream)
reads$downstream = as.factor(reads$downstream)



reads_leader = reads[ reads$startPos <= 100,]

reads_no_leader = reads[ reads$startPos > 100 ,]




outfile0 = paste(resdir, "/transcript_error.pdf", sep="");
ggp = plotErrorViolin(reads, reads_no_leader,  inds1 = reads$upstream=="leader" & !is.na(reads$downstream))
try(ggsave(outfile0, plot=ggp, width = 35, height = 15, units = "cm"))

 outfile0 = paste(resdir, "/transcript_expression.pdf", sep="");
if(length(type_nme)>1){
 ggps = list()
 ggps[[1]]<-.plotGeneExpression(reads_leader, nme1 = "comb_l", target = "source", limit = 10, mincount = 0.1)
 ggps[[2]]<-.plotGeneExpression(reads_no_leader, "comb_l", "source", limit = 10, mincount = 0.1)
 try(ggsave(outfile0, plot=marrangeGrob(ggps,nrow = 1, ncol = 2), width = 50, height = 15, units = "cm"))
}

 #ggp1 = plotErrorViolin(reads, inds1 = reads$start < 100,x = "reorder(comb_l, length)",  y = "length")





reads_ml = plotLengthHist(reads,t1, seqlen, type_nme)
  outfile1 = paste(resdir, "/length_hist.pdf", sep="");
  try(ggsave(outfile1, plot=reads_ml, width = 30, height = 30, units = "cm"))

