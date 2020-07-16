library(rhdf5)


##example of diff methylation 
#d1 = rep(100,3)
#d2 = rep(100,3)
#m1 = c(10,10,10)
#m2 = c(10,20,30)
#dfi= cbind(d1,d2,m1,m2)
#.fishert<-function(v){
#  ft = fisher.test(matrix(v,nrow=2,ncol=2))
#  ft$p.value
#}
#pv = apply(dfi,1,.fishert)
#cbind(dfi,pv)


depth_inds_i = grep("d",names(dfi))
meth_inds_i = grep("m",names(dfi))
res = t(apply(dfi, 1, betaBinomialP2, depth_inds_i, meth_inds_i,1,2,binom=F, log=F))
res = data.frame(t(apply(dfi, 1, betaBinomialP2, depth_inds_i, meth_inds_i,1,2,binom=F, log=F)))
names(res) =  c("pv_less","pv_more","ratio_control","ratio_case")?fie
##


###H5 info

#chroms =  attributes[which(names(attributes)=="chroms")][[1]]
#chroms=chroms[which(names(chroms)=="MT")]

depths=.readH5All(transcripts,attributes,filenames, thresh = 1000, chrs=NULL)
depth = .combineTranscripts(depths, attributes)
depth = .transferAttributes(depth, attributes)
attr1 = attributes(depth)

DE = DEdepth(depth, control_names, infected_names, tojoin=1:3)

#DEs1 = DEdepth(depths[[1]], control_names, infected_names )
#DEs =   lapply(depths, DEdepth, control_names, infected_names,tojoin=1:3)
#DE = .combineTranscripts(DEs, attributes)
DE = DE[,grep(exclude_nme,names(DE),inv=T)]
DE = .transferAttributes(DE, attr1)



outp="depth_association.csv"
.write(DE ,resdir,outp, numeric_inds = which(!(names(DE)%in% c("IDS", "pos","base"))))
par(mfrow = c(2,2))
nme_ = grep("pv_", names(DE),v=T)
for(k in (nme_)){
.qqplot1(DE,k, min.p= 1e-200,main=k)
}

.plot1<-function(DE, nme,log="", extra=0){
  topl = DE[,which(names(DE)%in%nme)]
  topl=apply(topl,c(1,2),function(x) x+extra)
  print(cor(topl))
  plot(topl,log=log)
}
.hist1<-function(DE, nme,breaks){
  topl = DE[,which(names(DE)%in%nme)]
 hist(topl[,1],br=breaks,main=nme[1])
 hist(topl[,2],br=breaks,main=nme[2])
}
nme_1 = grep("pv_less", names(DE),v=T)
nme_2 = grep("pv_more", names(DE),v=T)
nme_3 = grep("ratio_control", names(DE),v=T)
nme_4 = grep("ratio_case", names(DE),v=T)



par(mfrow=c(2,2))
.plot1(DE, nme_1,log="xy",1e-200)
.plot1(DE, nme_2,log="xy",1e-200)

.plot1(DE, nme_3)
.plot1(DE, nme_4)

dev.off()
.hist1(DE,nme_3,br=seq(0,100))
.hist1(DE,nme_4,br=seq(0,100))


par(mfrow=c(2,1))
chrom = NULL
min.p = 1e-200
.vis(DE, nme_2[1],min.p = min.p,log=F, chroms=chrom,usePos=F)
.vis(DE, nme_2[2],min.p =min.p,log=F, chroms=chrom, usePos = F)
.vis(DE, nme_1[1],min.p = min.p,log=F, chroms=chrom)
.vis(DE, nme_1[2],min.p = min.p,log=F, chroms=chrom)

par(mfrow=c(2,1))
chrom="chrM"
.vis(DE, nme_4[1],min.p = min.p,log=F, chroms=chrom, funct = NULL)
.vis(DE, nme_4[2],min.p = min.p,log=F, chroms=chrom, funct = NULL)
.vis(DE, nme_3[1],min.p = min.p,log=F, chroms=chrom, funct = NULL)
.vis(DE, nme_3[2],min.p = min.p,log=F, chroms=chrom, funct = NULL)



DE1 = DE[,names(DE)%in% nme_2]
more_lt = DE[apply(DE1,1,max)<1e-200,1:3]
getlev(more_lt$chrs)



depth1 = .appendGeneNamesToDepth(depth, transcripts, sort=NULL)
#depth2 = .appendGeneNamesToDepth(depth, transcripts, sort="pv2")

sigChr1 = findSigChrom(depth, thresh=1e-3, go_thresh=1e-3, nme="pv_less", nme2="gene_names")
sigChr2 = findSigChrom(depth1, thresh=1e-3, go_thresh=1e-3, nme="pv_more", nme2="gene_names")

##this is visualisation
  pdf("error_associations.pdf")
  pv_inds = grep("pv", names(depth))
  .qqplot(depth$pv1,min.p=1e-100)
  .qqplot(depth$pv2,min.p=1e-100)
   hist(depth[depth$pv1<1e-3,]$base) ## this is up in controls
   hist(depth[depth$pv2<1e-3,]$base) ## this is up in cases
   
   d1 = data.frame(depth[,which(names(depth) %in% c("pos","pv1"))],type=rep("1",dim(depth)[1]))
   d2 = data.frame(depth[,which(names(depth) %in% c("pos","pv2"))],type=rep("2",dim(depth)[1]))
   names(d1)[2] = "pv"; names(d2)[2] = "pv"
   ggp<-ggplot(rbind(d1,d2),aes(x=pos, y=pv,fill=type, color=type))+geom_point() + theme_bw()+ggtitle("pv")+scale_y_continuous(trans='log10')
   ggp
  .vis(depth,i=1,min.p=1e-20)
  .vis(depth,i=2,min.p=1e-20)
  dev.off()

 # di_2 =   depth$pv2<1e-5 | depth$pv1<1e-5
  

  
  
  
