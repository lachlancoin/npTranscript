
##TRANSCRIPT COVERAGE ANALYSIS
min_t_count = 2
max_h = 0

h5file = "0clusters.h5"
#h5f = H5Fopen(cluster_file)
header = h5read(h5file,"header")

dinds = grep("depth", header)
  einds = grep("errors", names(clusters))

for(ik in 1:length(transcripts_all)){
  print(paste("ik",ik))
  transcripts_ik = transcripts_all[[ik]]
  


  .getRatios(clusters, dinds, einds)
  
  nme=rev(strsplit(infiles,"\\.")[[1]])[2]
  outfile1 = paste(resdir,".", ik, "/cov", nme, "1.pdf", sep="");
  outfile2 = paste(resdir,".", ik, "/cov", nme, "2.pdf", sep="");

  ml = makeCovPlots(h5f, header, transcripts_ik[1:10,], total_reads, t, type_nme, logy=F, rawdepth =T, span = 0.02) 
  ml1 = makeCovPlots(clusters, transcripts_ik[1:10,], total_reads,  t, type_nme, logy=F, rawdepth =T, span = 0.02, xlim = list(c(0,100), c(25000,30000)))
  
  try(ggsave(outfile1, plot=ml, width = 30, height = 30, units = "cm"))
  try(ggsave(outfile2, plot=ml1, width = 30, height = 30, units = "cm"))

if(HEATMAP){
  outfile2 = paste("results/hm_", nme,".pdf", sep="");
  dinds = grep("depth", names(clusters))
  pdf(outfile2)
  
  numt = length(which(transcripts$countTotal>min_t_count))
  xlim = list(c(1,10000), c(20000,30000))
  for(jk in 1:length(dinds)){
    hm = plotHeatmap(clusters, transcripts[1:numt,,drop=F], jk, xlim = xlim, max_h=max_h,logHeatmap = T, title = type_nme[jk], featureName = "depth")
    if(length(einds)>0){
    hm = plotHeatmap(clusters, transcripts[1:numt,,drop=F], jk, xlim = xlim, max_h=max_h,logHeatmap = F, title = type_nme[jk], featureName = "ratios")
    }
    }
  dev.off();
}
}

