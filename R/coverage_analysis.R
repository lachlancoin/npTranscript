
##TRANSCRIPT COVERAGE ANALYSIS
min_t_count = 20
maxt_ =10
max_h = unlist(lapply(transcripts_all, function(transcripts)  max(transcripts[,grep('count[0-9]', names(transcripts)), drop=F])))


sumT_all = matrix(nrow = length(infiles), ncol = length(type_nme))
dimnames(sumT_all) = list(infiles, type_nme)

h5file = grep("clusters.h5",dir(),v=T)[1]

#h5f = H5Fopen(cluster_file)
header = h5read(h5file,"header")

dinds = grep("depth", header)
einds = grep("errors", header)

for(ik in 1:length(transcripts_all)){
  print(paste("ik",ik))

  inds_ik = transcripts_all[[ik]]$countTotal>min_t_count
  if(length(which(inds_ik))<1) next;
  transcripts_ = transcripts_all[[ik]][inds_ik,,drop=F]
  if(dim(transcripts_)[1]> maxt_) transcripts_ = transcripts_[1:maxt_,]
  pos1 =  c(0:10000, 20000:seqlen)
   clusters_ = readH5(h5file, header, transcripts_, pos =NULL, span = 0.01, cumul=F)
nme1 = names(transcripts_all)[[ik]] #as.character(transcripts_$type[1])

  #.getRatios(clusters, dinds, einds)
  
  nme=rev(strsplit(infiles,"\\.")[[1]])[2]
  ml1 = makeCovPlots(clusters_, header,  total_reads,  t, paste(type_nme,nme1), logy=T, rawdepth =T,  xlim = list(c(0,10000), c(20000,seqlen)), fill=F)
  ml = makeCovPlots(clusters_, header, total_reads, t, paste(type_nme,nme1), logy=T, rawdepth =T) 

  outfile1 = paste(resdir, "/coverage",".", ik,".",  nme,".",  nme1,  ".1.pdf", sep="");
  outfile2 = paste(resdir,"/coverage", ".", ik, ".", nme,".",  nme1, ".2.pdf", sep="");

  try(ggsave(outfile1, plot=ml, width = 30, height = 30, units = "cm"))
  try(ggsave(outfile2,plot=ml1, width = 30, height = 30, units = "cm"))

if(HEATMAP){
  outfile2 = paste("results/coverage_heatmap_", nme,nme1, ".pdf", sep="");
  
  pdf(outfile2)
  
  numt = length(which(transcripts_$countTotal>min_t_count))
  xlim = list(c(1,200), c(1,10000), c(20000,30000))
  for(jk in 1:length(dinds)){
    hm = plotHeatmap(h5file,header,  transcripts_, jk, xlim = xlim, max_h=max_h[ik],logHeatmap = T, title = paste(type_nme[jk], nme1), featureName = "depth")
    if(length(einds)>0){
    hm = plotHeatmap(h5file,header, transcripts_, jk, xlim = xlim, max_h=max_h[ik],logHeatmap = F, title = paste(type_nme[jk],nme1), featureName = "ratios")
    }
    }
  dev.off();
}
}

