
##TRANSCRIPT COVERAGE ANALYSIS
min_t_count = 2
max_h = 0
if(HEATMAP){
for(ik in 1:length(infiles)){
  transcripts = transcripts_all[[ik]]
 
 # exons = exons_all[[ik]]
  clusters = data.frame(as.matrix(read.table( infiles[ik],sep=",",head=T)))  ##
  dinds = grep("depth", names(clusters))
  max_h = max(max_h,max(clusters[,dinds],na.rm=T))
}
}

for(ik in 1:length(infiles)){
  print(paste("ik",ik))
  transcripts = transcripts_all[[ik]]
  
  # exons = exons_all[[ik]]

  clusters = data.frame(as.matrix(read.table( infiles[ik],sep=",",head=T)))  ##
  dinds = grep("depth", names(clusters))
  if(max_h==0)  max_h = max(clusters[,dinds],na.rm=T)

  if(addOne){
    pos_ind = which(names(clusters)=="pos")
    clusters[,pos_ind] = clusters[,pos_ind]+1
  }
  clusters[,2] = as.factor(clusters[,2])
  
  
  
  if(mult>1)  clusters = apply(clusters,c(1,2), function(a) a/mult)
  
  dinds = grep("depth", names(clusters))
  einds = grep("errors", names(clusters))
  if(length(einds)>0){
     ratios=  data.frame(matrix(nrow = dim(clusters)[1], ncol = length(dinds)))
    for(i in 1:length(dinds)){
      ratios[,i] = clusters[,einds[i]]/clusters[,dinds[i]]
      ratios[ clusters[,dinds[i]]==0,i] = NA
    }
     names(ratios) = paste("ratios",1:length(dinds), sep="")
    clusters = cbind(clusters,ratios)
  }
  
  nme=rev(strsplit(infiles[ik],"\\.")[[1]])[2]
  outfile1 = paste(resdir, "/cov", nme, "1.pdf", sep="");
  outfile2 = paste(resdir, "/cov", nme, "2.pdf", sep="");
if(COVERAGE){  
  ml = makeCovPlots(clusters, transcripts[1:10,], total_reads, t, type_nme, logy=F, rawdepth =T, span = 0.0) 
  ml1 = makeCovPlots(clusters, transcripts[1:10,], total_reads,  t, type_nme, logy=F, rawdepth =T, span = 0.02, xlim = list(c(0,100), c(25000,30000)))
  
  try(ggsave(outfile1, plot=ml, width = 30, height = 30, units = "cm"))
  try(ggsave(outfile2, plot=ml1, width = 30, height = 30, units = "cm"))
}
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

