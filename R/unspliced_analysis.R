if(is.null(reads)){


reads=  data.frame(read.table(infilesReads, sep="\t",  head=T))
#names(reads) =  strsplit("readID\tclusterId\tsubID\tsource\tlength\tstart_read\tend_read\ttype_nme\tchrom\tstartPos\tendPos\tbreakStart\tbreakEnd\tbreakStart2\tbreakEnd2\terrorRatio\tupstream\tdownstream\tupstream2\tdownstream2\tstrand\tbreaks","\t")[[1]]


reads$source =  factor(reads$source, levels = 0:(length(type_nme)-1), labels=type_nme)

reads$type = as.factor(reads$type)

reads_l = list()
for(i in 1:length(transcripts_all)){
reads_l[[i]] = reads[reads$clusterId %in% transcripts_all[[i]]$ID,,drop=F]
}
names(reads_l)  = names(transcripts_all)

}


 getlev(as.character(reads[reads$downstream=="ORF10",]$upstream))
breakInds = which(names(reads) %in% c("breakStart", "breakEnd"))
breakInds2 = which(names(reads) %in% c("breakStart2", "breakEnd2"))

 ORF1ab_ORF10 =.stEndKmer(.expandLev(getlev(reads[reads$downstream=="ORF10" & reads$upstream=="ORF1ab",breakInds])))
 ORF7a_ORF10 =.stEndKmer(.expandLev(getlev(reads[reads$downstream2=="ORF10" & reads$upstream2=="ORF7a",breakInds2])))
 ORF6_ORF10 =.stEndKmer(.expandLev(getlev(reads[reads$downstream2=="ORF10" & reads$upstream2=="ORF6",breakInds2])))

kmers_l = getKmer(fasta[[1]], lev_a[,1],v = 0:8)


 head(getlev(reads[reads$downstream2=="ORF10" & reads$upstream2=="ORF7a",breakInds2]))
 head(getlev(reads[reads$downstream2=="ORF10" & reads$upstream2=="ORF6",breakInds2]))



a = plotHist((reads[reads$downstream=="ORF10" & reads$upstream=="ORF1ab",]$breakEnd), br=seq(25000,seqlen, by=10))





 plotHist(reads[reads$breakSt<0 & reads$start_read > 100,]$startPos,1:seqlen,t, ylog=F, xlim = c(25000,30000))

binsize = 10
breaks = seq(1,seqlen+binsize, binsize)-0.5
ml  = plotHist(reads[reads$downstream=="ORF10" & reads$upstream=="leader",]$breakEnd,breaks, t, ylog=F, xlim = c(28000,seqlen))
ggsave("ORF10_breakEnd.pdf", plot=ml, width = 30, height = 30, units = "cm")

kmer_background8 = getlev(getKmer(fasta[[1]], 1:seqlen,v = 0:8))




 unspliced = reads[reads$breakStart<0 ,]

plotHist(reads[reads$type=="5_3",]$startPos,1:seqlen,t)
plotHist(reads$startPos[reads$start_read > 200], 1:seqlen,t, )
plotHist(unspliced$startPos[unspliced$start_read > 200], 1:seqlen,t, )

kComp0 = .getKComp(fasta[[1]], unspliced, k=8)
compLevs0 = lapply(kComp0, .compareLevs, kmer_background8)
names(compLevs0) = names(kComp0)

compLevs1 = list()
compLevs2 = list()
diff = unspliced$length-unspliced$end_read

#kCompSub_ = .getKComp(fasta[[1]], unspliced,inds =which(unspliced$start_read>29000),  k=8)
#comp1 = .compareLevs(kCompSub_, kmer_background8)

kCompSub = .getKComp(fasta[[1]], unspliced,inds =which(unspliced$start_read>200),  k=8)
kCompSub1 = .getKComp(fasta[[1]], unspliced,inds =which(diff>50),  k=8)
for(i in 1:length(kCompSub)) {
compLevs1[[i]] = .compareLevs( kCompSub[[i]],kComp0[[i]])
compLevs2[[i]] = .compareLevs( kCompSub1[[i]],kComp0[[i]])
}
names(compLevs1) = names(kComp0)
names(compLevs2) = names(kComp0)

#mat = gregexpr('taaacgaac', fastaseq)
#mat = gregexpr('aaacgaaca', fastaseq)

