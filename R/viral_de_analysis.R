library(abind)
library(ggrepel)
library(writexl)
library(grDevices)


attributes = attributes(transcripts)
info = attr(transcripts, "info")
filenames = info
count_names =  paste("count",info, sep="_")
names(transcripts)[grep("count[0-9]", names(transcripts))] =count_names
names(transcripts)[grep("errors[0-9]", names(transcripts))] = paste("errors",info, sep="_")
names(transcripts)[grep("error_ratio[0-9]", names(transcripts))] = paste("error_ratio",info, sep="_")
todo=.getAllPairwiseComparisons(info)

## first DE of expression
pdf(paste(resdir,"DE.pdf",sep="/"))
DE1 = lapply(todo, function(x) .processDE1(transcripts, count_names,x[1],x[2], resdir, top=5, pthresh=1e-2,plot=F))
if(length(DE1)==1) DE1 = DE1[[1]]
volcanos = lapply(DE1, .volcano, top=10)
lapply(volcanos, function(x) print(x))
dev.off()
names(DE1) = unlist(lapply(DE1, function(x) gsub("_merged","",gsub("_pass","",gsub("_infected","",attr(x,"nme"))))))
names(volcanos) = names(DE1)  
write_xlsx(DE1[[1]], paste(resdir, "DE.xlsx",sep="/"))

##DIFF METH  ON COMBINED
if(length(todo)>0){
depth_combined=readH5_c("0.clusters.h5",transcripts,filenames, thresh=as.numeric(getOption("np.depth_thresh","100")))
  #.readH5All(transcripts,attributes,filenames, thresh = 100,  readH5_ = readH5_c)[[1]]
DE2 = lapply(todo, function(x) .processDM(depth_combined, filenames, filenames[x[1]],filenames[x[2]], method=".exact", thresh =100))
DE2 = DE2[unlist(lapply(DE2, function(x) dim(x)[[1]]))>0]
for(i in 1:length(DE2)) attr(DE2[[i]],"nme") = names(DE2)[i]
#names(DE2[[1]]) = unlist(lapply(DE2[[1]], function(x) .shorten(attr(x,"nme"),31)))
#volcanos = lapply(DE2, .volcano, pthresh = 1e-5)
volcanos = lapply(DE2, .volcano, top=5, logFCthresh = 0.5)
pdf(paste(resdir,"DM_combined.pdf",sep="/"))
lapply(volcanos, function(x) print(x))
dev.off()
write_xlsx(.joinSS(lapply(DE2,.xlim,1e-5), sort=T), paste(resdir, "DM_combined.xlsx",sep="/"))

ggp  = .plotError(depth_combined,t,extend=F,range =28250:30000, method="bayes", ci=0.995, thresh = 1000, adj=T, diff_thresh =0.2, pval_thresh = 1e-6,lower_thresh =  0.2 )
outfile1 = paste(resdir, "/error_profiles.pdf", sep="");
try(ggsave(outfile1, plot=ggp, width = 60, height = 20, units = "cm"))
}

##ALL TRANSCRIPTS AND ALL POSITIONS COMPARED BETWEEN SAMPLES
filenames=info
if(length(todo)>0){
depth=readH5_h("0.clusters.h5",transcripts,filenames, thresh = as.numeric(getOption("np.depth_thresh","100")))
#  .readH5All(transcripts,attributes,filenames, thresh = 1000, chrs=NULL, readH5_ = readH5_h)[[1]]
#depth = .transferAttributes(depth, attributes)
DE2 = lapply(todo, function(x) .processDM(depth,filenames,  filenames[x[1]],filenames[x[2]], method=".exact", thresh =500))
DE2 = DE2[unlist(lapply(DE2, function(x) dim(x)[[1]]))>0]
for(i in 1:length(DE2)) attr(DE2[[i]],"nme")=names(DE2)[[i]]
#names(DE2) = unlist(lapply(DE2, function(x) .shorten(attr(x,"nme"),31)))
volcanos = lapply(DE2, .volcano, top=5, logFCthresh =0.5)
pdf(paste(resdir,"DM_1.pdf",sep="/"))
lapply(volcanos, function(x) print(x))
dev.off()
write_xlsx(.joinSS(lapply(DE2,.xlim,1e-5), sort=F), paste(resdir, "DM_1.xlsx",sep="/"))

}



##LOOK at effect of splicing
#depths_combined_spliced= readH5_c( "0.clusters.h5", transcripts_all_splice[[1]], attributes, filenames, thresh = 100)

depths_combined_spliced=
  lapply(transcripts_all_splice,function(x)  readH5_c( "0.clusters.h5", x,  filenames, thresh = 100))
#depth_combined=readH5_c("0.clusters.h5",transcripts,filenames, thresh=as.numeric(getOption("np.depth_thresh","100")))

    #.readH5(x, "0.clusters.h5" , attributes, filenames, thresh=100,  readH5_ = readH5_c)[[1]])
depth_spliced = lapply(1:length(info), .mergeDepth, depths_combined_spliced)
names(depth_spliced)= info

ggp  = .plotError(depth_spliced[[1]][,,1:2],t,extend=T,log=F, range = 28100:28350, method="bayes", ci=0.95, 
                  thresh = 500, diff_thresh =0.1, lower_thresh =  0.2, pval_thresh = 1e-6, adj=T )

if(TRUE){
nme_spl = names(transcripts_all_splice)
todo2 = lapply((1:length(nme_spl))[-2], function(x)  c(2,x))
names(todo2) = paste(nme_spl[2],nme_spl[-2],sep=" v ")
todo2 = todo2[1]
#depth_spliced[[length(depth_spliced)+1]]= depth_spliced[[2]]+ depth_spliced[[3]]
#names(depth_spliced)[[length(depth_spliced)]] = paste(names(depth_spliced)[2:3],collapse="_")
pdf(paste(resdir,paste("DM_spliced","pdf",sep="."),sep="/"))
DE2_split_all = list()
for(k in 1:length(depth_spliced)){
  print(k)
  depthk = depth_spliced[[k]]
  fn = dimnames(depth_spliced[[k]])[[3]]

  DE2_split = lapply(todo2, function(x) .processDM(depth_spliced[[k]], fn, fn[x[1]],fn[x[2]], method="fisher.test", thresh_min =10000,plot=F))
  names(DE2_split) = lapply(DE2_split, function(x) .shorten(paste(names(depth_spliced)[k],attr(x,"nme"),sep=":")))
  DE2_split = DE2_split[unlist(lapply(DE2_split, function(x) dim(x)[[1]]))>0]
  volcanos_split = lapply(DE2_split, .volcano, top=10, prefix =names(depth_spliced)[k] )
  if(length(DE2_split)>0){
    lapply(volcanos_split, function(x) print(x))
    DE2_split_all = c(DE2_split_all, DE2_split)
  }
  #.extractFromDepth(depth_spliced[[k]], 1:5,ORF=NULL, pos=28265)
  
}
dev.off()
write_xlsx(lapply(DE2_split_all,.xlim,1e-5), paste(resdir, "DE_splicing.xlsx",sep="/"))
}

if(TRUE){
  ##src_i refers to which entry in info
depth_split=readH5_s("0.clusters.h5",transcripts[,1:2],filenames, thresh = 10,  src_i = 1)
todo1 = lapply(2:length(dimnames(depth_split)[[3]]), function(x) c(1,x))
names(todo1) = dimnames(depth_split)[[3]][-1]
fn = dimnames(depth_split)[[3]]
DE2_split = lapply(todo1, function(x) .processDM(depth_split, fn, fn[x[1]],fn[x[2]], method=".exact", thresh_min =0)[[1]])
#if(length(DE2_split)==1) DE2_split = DE2_split[[1]]
volcanos_split = lapply(DE2_split, .volcano, top=10, logFCthresh = 1)
pdf(paste(resdir,"DM_split.pdf",sep="/"))
lapply(volcanos_split, function(x) print(x))
dev.off()
ggp  = .plotError(depth_split,t,extend=F, range = 28250:30000, method="bayes", ci=0.995, thresh = 200, diff_thresh =0.2, lower_thresh =  0.2 )
outfile1 = paste(resdir, "/error_profiles.pdf", sep="");
try(ggsave(outfile1, plot=ggp, width = 40, height = 30, units = "cm"))
}
#.extractFromDepth(depth_split1, c(1,grep("^M;end",dimnames(depth_split1)[[3]])),ORF=NULL, pos=28254)


#depth = depths[[1]]
#DE2 = lapply(todo1, function(x) .processDM(depth,  x[1],x[2], method=.chisq, thresh =100))

#dev.off()
#pdf(paste(resdir,"DM_volcano.pdf",sep="/"))
#dev.off()

#DE2_=lapply(DE2,.xlim)
#names(DE2) = unlist(lapply(DE2, function(x) gsub("_merged","",gsub("_pass","",gsub("_infected","",attr(x,"nme"))))))
#names(DE2_) = names(DE2)
#write_xlsx(DE2_, paste(resdir, "DM1.xlsx",sep="/"))

#attr1 = attributes(depth)

#.extractFromDepth(depth, info[todo[[7]]],ORF="st;leader,M;end", pos=26536)

#.extractFromDepth(depth, info[todo[[13]]],ORF="", pos=29759)