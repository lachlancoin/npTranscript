###
##NEW WAY OF READING BREAKPONTS

if(FALSE){

.modifyH<-function(breakP_){
	breakP_$heatm = apply(breakP_$heatm,c(1,2), function(x) 1/exp(x))
	breakP_
}


plotRegion<-function(brPs, region, logT=F, pdf=T){
	file = paste(paste(region, collapse="_"),"pdf",sep= ".")
	if(pdf) pdf(file)
	for(i in 1:length(brPs)){
		
		plotBreakPIntrons(brPs[[i]],region = region, title=names(brPs)[[i]],subtitle = max(brPs[[i]]$heatm),  logT=logT)
	}
	if(pdf) dev.off()
}

findMaxSeqsAll<-function(brPs, fasta, region=c(60,80,1,28240,28260,1)){
 	res = list()
	for(i in 1:length(brPs)){
	 res[[i]] = findMaxSeqs(brPs[[i]], region = region, fasta, nme = names(brPs)[[i]])
	}
	names(res) = names(brPs)
	res
}
h5fileB = paste(currdir,"0.breakpoints.h5",sep="/") 
inputs = c("",grep("chr", grep("scores", h5ls(h5fileB)$group, v=T), inv=T, v=T))
names(inputs) = gsub("/","",inputs)
j1=1
brPs =  lapply(inputs, function(prefix) readBreakPointsH5(h5fileB,"chrMT007544", "cellular", j1, prefix=prefix))
brPs0 = brPs[1]
brPs = lapply(brPs[-1],.modifyH)
brPs = c(brPs0, brPs)
names(brPs)[[1]] = "reads"




.calcCor(brPs[[1]]$heatm, brPs[[2]]$heatm, offrow=0, offcol=-1, thresh =1)
.calcCor(brPs[[1]]$heatm, brPs[[5]]$heatm, offrow=12, offcol=12, thresh =1)


plot(as.vector(brPs[[1]]$heatm[indsk]), as.vector(brPs[[2]]$heatm[indsk+1]),log="xy")

summary(lm(as.vector(brPs[[1]]$heatm[indsk])~ as.vector(brPs[[2]]$heatm[indsk+1])))




.calcCor<-function(hm0,hm1, offrow=0, offcol=0, thresh=1){
  if(offrow>0) hm1 = hm1[-(1:offrow),]
  if(offcol>0) hm1 = hm1[,-(1:offcol)]
  if(offrow<0) hm0 = hm0[-(1:-(offrow)),]
  if(offcol<0) hm0 = hm0[,-(1:-(offcol))]
  inds = which(hm0>=thresh)
 print( summary(lm(hm0[inds]~ hm1[inds])))
#  plot(hm0[inds],hm1[inds])
  cor(hm0[inds], hm1[inds], use="pairwise.complete.obs")
}

regions = list(a =  c(60,80,1,28200,28400,1),
               a1=c(60,80,1,26440,26520,1),
		b=c(1000,10000,10,28274,29533,10),
		c=c(1000,7000,100,28274,29533,10),
		d = c(2020,2060,1,28500,28560,1),
		e= c(6930,6960,1,29450,29490,1),
		f= c(1800,2000,1,28274,29533,1)
)
plotRegion(brPs,region=regions$f)
plotRegion(brPs, region=c(1,100,1,25000,30000,100))
plotRegion(brPs, region=c(1,100,1,25000,30000,10), logT=F)

plotRegion(brPs[2], region=regions$a, logT=F, pdf=F)
maxseqs = findMaxSeqsAll(brPs,fasta,  region = regions$a)
lapply(maxseqs, names)

plotRegion(brPs[2], region=regions$b, logT=F, pdf=F)
maxseqs = findMaxSeqsAll(brPs,fasta,  region = regions$b)
lapply(maxseqs, names)

plotRegion(brPs, region=regions$c, pdf=T)
plotRegion(brPs, region=regions$e, pdf=T)
plotRegion(brPs, region=regions$f, pdf=T)

 plotRegion(brPs,region=c(6700,7300,10,28274,29533,10),logT=F,pdf=T)
 plotRegion(brPs[1],region=c(1950,2000,1,29450,29500,1),logT=F,pdf=T)


plotRegion(brPs, region=c(60,80,1,28240,28260,1), logT=F)
plotRegion(brPs[1], region=c(60,80,1,29500,29600,1), logT=F, pdf=F)
maxseqs = findMaxSeqsAll(brPs[1],fasta,  region = c(60,80,1,29533,1))
region = 

plotRegion(brPs, region, logT=F)
maxseqs = findMaxSeqsAll(brPs[1],fasta,  region )

maxseqs = findMaxSeqsAll(brPs[1],fasta,  region = c(1700,2300,1,28274,29533,1))

}
# brP1 = readBreakPointsH5(h5file,"chrMT007544", "virion", 0)
#plots_i =  plotBreakPIntrons(brP1)
#region =  c(1,5000,100,25000,30000,100)
#plots = plotBreakPIntrons(brP1, t,fimo,region =  region,   logT=logT, title= names(special)[i1], subtitle = type_nme[i])

##BREAKPOINT ANALYSIS
#infilesBr = grep("breakpoints.", dir(), v=T)
infilesBr0 = paste("0.breakpoints.",type_nme,".0.txt.gz",sep="")
infilesBr1 = paste("0.breakpoints.",type_nme,".1.txt.gz",sep="")


breakPs = vector("list", length = length(infilesBr0))
incl = rep(T, length(infilesBr0))
for(i in 1:length(infilesBr0)) {
  brP0 = try( readBreakPoints(infilesBr0[i], i, addOne = FALSE))
  brP1 = try( readBreakPoints(infilesBr1[i], i, addOne = FALSE))
  
 
  if(is.null(brP0) ){
    incl[i] = F
   # breakPs[[i]] = NULL
  }else if(is.null(brP1)){
    breakPs[[i]] = brP0;
  }  else{
    breakPs[[i]] = .mergeBreak(list(brP0, brP1))
  }
}
breakPs = breakPs[incl]
names(breakPs) = sub("0.breakpoints.","",sub(".0.txt.gz" ,"" ,infilesBr0))[incl]



 total_reads = rep(1, length(breakPs))

genes = grep("leader",grep("UTR", grep("none",t$gene[-2], inv=T, v=T), inv=T, v=T),v=T, inv=T)

endcs =getGeneBP(t, genes,left = 300, right=10)
endcs =getGeneBP(t, genes,left = 300, right=10)

endcs_n = dim(endcs)[1]
endcs_1 = cbind(endcs, rep(1, endcs_n))
endcs_1 = data.frame(t(cbind(t(matrix(rep(c(1,100,1),endcs_n), ncol=endcs_n)), endcs_1)))
endcs_2 = cbind(endcs, rep(1, endcs_n))
endcs_2 = data.frame(t(cbind(t(matrix(rep(c(1,10000,10),endcs_n), ncol=endcs_n)), endcs_2)))

endcs_b = getGeneBP(t, genes,left = 1000, right=10, left_buffer = 10)
endcs_b_1 = cbind(endcs_b, rep(1, dim(endcs_b)[1]))
endcs_b_1 = data.frame(t(cbind(t(matrix(rep(c(1,10000,1),dim(endcs_b)[1]), ncol=dim(endcs_b)[1])), endcs_b_1)))

#UTR = c(1,10000,1,29675-300, 29893,1)

endcs_2[5,dim(endcs_2)[2]] = length(fasta[[1]])
endcs_1[5,dim(endcs_1)[2]] = length(fasta[[1]])

#endcs_2[[length(endcs_2)+1]] = UTR
#names(endcs_2)[length(endcs_2)] = "5UTR"

special = data.frame(
  "Leader_100b'"  = c( minpos,100,1, maxpos-10000,maxpos,100),
  "Leader_10kb"= c(1,10000,50, maxpos-10000,maxpos,100),
  "full"= c(1,seqlen,1000, 1,seqlen,1000),
  "start"= c(1,10000,100, 1,10000,100),
"end"= c(27000,seqlen,100, 27000,seqlen,100),
"end1"= c(27000,seqlen,50, 27000,seqlen,50),
"end2"= c(27700,27950,10, 29000,29500,10),
"end3"= c(27300,27500,10, 29400,29600,10)
)


if(RUN_ALL){
	outfile3a = paste(resdir, "/expression1.pdf", sep="");
	outfile3b = paste(resdir, "/expression2.pdf", sep="");
	prot_file =  paste(resdir, "/expression.csv", sep="");
	break_file =  paste(resdir, "/breaks.csv", sep="");
	print("special")
	 plotAllHM(special, "full" , resdir,  breakPs, t, fimo, total_reads, type_nme = names(breakPs), log=T)	

	print("endcs_1")
	 plotAllHM(endcs_1, "regional_100b_leader" , resdir,  breakPs, t, fimo, total_reads, type_nme = names(breakPs),  log=T)	
print("endcs_2")
	plotAllHM(endcs_2, "regional_10kb_leader" , resdir,  breakPs, t, fimo, total_reads, type_nme = names(breakPs), plotHM = F,  log=T)	
	
	

	all = list(all = c(1,10000,1,20000,30000,1))
	hmClust_c = getHMClust(breakPs, all,fasta,t, mind = c(1,1))
  merged = hmClust_c$merged
	prot_l = apply(merged[,1:3], c(1,2), function(x) nchar(x)-1)
	names(prot_l) = c("prot0_length", "prot1_length", "prot2_length")
	write.csv(cbind(merged, prot_l),file = prot_file, quote=FALSE, sep=",",row.names=T, col.names=T)
	
	breaks_all = hmClust_c$merged_se
	print_inds = c(grep("^genes$", names(breaks_all)),grep("start", names(breaks_all)),grep("end", names(breaks_all)))
	br_all = breaks_all[	order(breaks_all$genes,breaks_all$start, breaks_all$end),print_inds]
	chrom = (rep(0, dim(br_all)[1]))
 write.table(cbind(chrom,br_all),break_file, quote=FALSE, row.names=F,col.names=T,sep=",")
	if(length(type_nme)>1){
		ml = plotHMClust1(hmClust_c, total_reads,names(breakPs), nudge_y = 0.0, nudge_x =0.25, logT = T, plotDepth = F)
		ml2 = plotHMClust1(hmClust_c, total_reads,names(breakPs), nudge_y = 0.0, nudge_x =0.25, logT = T, plotDepth = T)
		n = length(breakPs)
		try(ggsave(outfile3a, plot=ml, width = 15, height = 15*(n+1)*(n/4), units = "cm",limitsize=F))
		try(ggsave(outfile3b ,plot=ml2, width = 15, height = 15*(n+1)*(n/4), units = "cm", limitsize=F))
	}
	
#	hmClust = getHMClust(breakPs, endcs_2)
}else{
  orf10 = data.frame(endcs_2[9])
  #orf10[[1]][1]=c(5000)
ml1 = plotAllHM(orf10, breakPs, t, fimo, total_reads, type_nme, NULL,NULL, log=T, depth=T)	
hmClust_b = getHMClust(breakPs, endcs_b_1,fasta,t, mind = c(2,2))
hmCplotHMClust1(hmClust_b, total_reads)

regions = list(all=c(6300,6500,1,27760,  27763,1))


hmClust = getHMClust(breakPs, regions,fasta,t, mind = c(5,5))
coords=  hmClust$coords$cell



#prots = getProt(fasta,t[11,3:4] )

orf10[[1]][4:5]=c(29500, 29675)
ml1 = plotAllHM(orf10, breakPs, t, fimo, total_reads, type_nme, NULL,NULL, log=T, depth=T)	

ml1
print(ml1)
}
