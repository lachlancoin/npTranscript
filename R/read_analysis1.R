
if(FALSE){
ml = plotJoins(reads, thresh = c(10,200),  t=t, ylog=T, xlim = c(20000,seqlen))
subr = reads[reads$breakStart == -1 & reads$breakEnd==-1 &  reads$start_read>50,]
ggp<- plotJoins1(subr, t=t, ylog=T, xlim=c(0,seqlen)) #, clusterId="ID0.61")
subr1 = reads[reads$breakStart == -1 & reads$breakEnd==-1 &  reads$length - reads$end_read>50,]
ggp<- plotJoins1(subr1, t=t, ylog=T, xlim=c(28000,seqlen), right=T) #, clusterId="ID0.61")

levs=list()
lev_source = levels(subr$source)
for(i in 1:length(lev_source)){
levs[[i]] = getlev(subr[subr$source==lev_source[i],,drop=F]$startPos)
}
names(levs) = lev_source
}
