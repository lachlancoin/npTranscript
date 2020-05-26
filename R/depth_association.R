###H5 info
depth=.readH5All(transcripts,chroms= attributes[which(names(attributes)=="chroms")][[1]] ,thresh = 100)

depth1 = .appendGeneNamesToDepth(depth, transcripts, sort=NULL)
#depth2 = .appendGeneNamesToDepth(depth, transcripts, sort="pv2")

sigChr1 = findSigChrom(depth1, thresh=1e-3, go_thresh=1e-3, nme="pv1", nme2="gene_names")
sigChr2 = findSigChrom(depth1, thresh=1e-3, go_thresh=1e-3, nme="pv2", nme2="gene_names")

##this is visualisation
  pdf("error_associations.pdf")
  pv_inds = grep("pv", names(depth))
  .qqplot(depth$pv1,min.p=1e-100)
  .qqplot(depth$pv2,min.p=1e-100)
   hist(depth[depth$pv1<1e-3,]$base) ## this is up in controls
   hist(depth[depth$pv2<1e-3,]$base) ## this is up in cases
  .vis(depth,i=1,min.p=1e-20)
  .vis(depth,i=2,min.p=1e-20)
  dev.off()

 # di_2 =   depth$pv2<1e-5 | depth$pv1<1e-5
  
  
chr_inds=attr(depth,"chr_inds")
chroms = attr(depth,"chroms")
chroms1 = names(chroms)
mi = match(chr_inds, chroms)
mi1 = unlist(lapply(mi, function(x) chroms1[x]))
pos100k = apply(cbind(as.character(mi1),as.character(binsize*floor(depth1$pos/binsize))),1,paste,collapse=".")

if(mergeByPosAndGene){
  pos100k = apply(cbind(as.character(mi1),as.character(depth1$gene_names), as.character(binsize*floor(depth1$pos/binsize))),1,paste,collapse=".")
  
}
depth2 = cbind(depth1, pos100k)
sigChr1 = findSigChrom(depth2, thresh=1e-5, go_thresh=1e-3, nme="pv1", nme2="pos100k")
DE_sig1 = DE2[DE2$pos1M %in% sigChr1$chrs,]
sigChr2 = findSigChrom(depth2, thresh=1e-5, go_thresh=1e-3, nme="pv2", nme2="pos100k")
DE_sig2 = DE2[DE2$pos1M %in% sigChr2$chrs,]

