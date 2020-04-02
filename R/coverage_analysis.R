
##TRANSCRIPT COVERAGE ANALYSIS
min_t_count = 20
max_h = 0

h5file = "0clusters.h5"
#h5f = H5Fopen(cluster_file)
header = h5read(h5file,"header")

dinds = grep("depth", header)
  einds = grep("errors", header)

for(ik in 1:length(transcripts_all)){
  print(paste("ik",ik))


  transcripts_ = transcripts_all[[ik]][transcripts_all[[ik]]$countTotal>min_t_count]
   clusters_ = readH5(h5file, header, transcripts_, span = 0.005)


  #.getRatios(clusters, dinds, einds)
  
  nme=rev(strsplit(infiles,"\\.")[[1]])[2]
 
  ml = makeCovPlots(clusters_, header, total_reads, t, type_nme, logy=T, rawdepth =T) 
  ml1 = makeCovPlots(clusters_, header,  total_reads,  t, type_nme, logy=T, rawdepth =T,  xlim = list(c(0,1000), c(20000,seqlen)))
   outfile1 = paste(resdir, "/cov",".", ik, nme, "1.pdf", sep="");
  outfile2 = paste(resdir,"/cov", ".", ik, nme, "2.pdf", sep="");

  try(ggsave(outfile1, plot=ml, width = 30, height = 30, units = "cm"))
  try(ggsave(outfile2, plot=ml1, width = 30, height = 30, units = "cm"))

if(HEATMAP){
  outfile2 = paste("results/hm_", nme,".pdf", sep="");
  
  pdf(outfile2)
  
  numt = length(which(transcripts$countTotal>min_t_count))
  xlim = list(c(1,200), c(1,10000), c(20000,30000))
  for(jk in 1:length(dinds)){
    hm = plotHeatmap(h5file,header,  transcripts_, jk, xlim = xlim, max_h=max_h[ik],logHeatmap = T, title = type_nme[jk], featureName = "depth")
    if(length(einds)>0){
    hm = plotHeatmap(h5file,header, transcripts_, jk, xlim = xlim, max_h=max_h[ik],logHeatmap = F, title = type_nme[jk], featureName = "ratios")
    }
    }
  dev.off();
}
}

