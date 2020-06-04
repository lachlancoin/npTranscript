###H5 info
chroms =  attributes[which(names(attributes)=="chroms")][[1]]
chroms=chroms[which(names(chroms)=="MT")]
depth=.readH5All(transcripts,chroms=chroms ,thresh = 1000)

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
   
   d1 = data.frame(depth[,which(names(depth) %in% c("pos","pv1"))],type=rep(1,dim(depth)[1]))
   d2 = data.frame(depth[,which(names(depth) %in% c("pos","pv2"))],type=rep(2,dim(depth)[1]))
   names(d1)[2] = "pv"; names(d2)[2] = "pv"
   ggp<-ggplot(rbind(d1,d2),aes(x=pos, y=pv1,fill=type))+geom_point() + theme_bw()+ggtitle("pv")+scale_y_continuous(trans='log10')
   
  .vis(depth,i=1,min.p=1e-20)
  .vis(depth,i=2,min.p=1e-20)
  dev.off()

 # di_2 =   depth$pv2<1e-5 | depth$pv1<1e-5
  
