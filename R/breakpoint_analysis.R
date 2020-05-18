##BREAKPOINT ANALYSIS
breakPs = list()
for(i in 1:length(infilesBr))  breakPs[[i]] = readBreakPoints(infilesBr[i], i, addOne = FALSE) 

if(length(breakPs)==length(type_nme)){
 names(breakPs) = type_nme[1:length(breakPs)] 
}else{
 names(breakPs) = infilesBr

}
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
		ml = plotHMClust1(hmClust_c, total_reads,type_nme, nudge_y = 0.0, nudge_x =0.25, logT = T, plotDepth = F)
		ml2 = plotHMClust1(hmClust_c, total_reads,type_nme, nudge_y = 0.0, nudge_x =0.25, logT = T, plotDepth = T)
		n = length(type_nme)
		try(ggsave(outfile3a, plot=ml, width = 15, height = 15*(n+1)*(n/4), units = "cm"))
		try(ggsave(outfile3b ,plot=ml2, width = 15, height = 15*(n+1)*(n/4), units = "cm"))
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
