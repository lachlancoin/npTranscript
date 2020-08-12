attributes = attributes(transcripts)
info = attr(transcripts, "info")
count_names =  paste("count",info, sep="_")
names(transcripts)[grep("count[0-9]", names(transcripts))] =count_names
names(transcripts)[grep("errors[0-9]", names(transcripts))] = paste("errors",info, sep="_")
names(transcripts)[grep("error_ratio[0-9]", names(transcripts))] = paste("error_ratio",info, sep="_")

todo = list()
for(i in 2:length(info)){
  for(j in 1:(i-1)){
    todo[[length(todo)+1]] = c(i,j)
  }
}
names(todo) = unlist(lapply(todo, function(x) paste(info[x],collapse=" v ")))


## first DE of expression
pdf(paste(resdir,"DE.pdf",sep="/"))
DE1 = lapply(todo, function(x) .processDE1(transcripts, count_names,x[1],x[2], resdir, top=5, pthresh=1e-2,plot=F))
volcanos = lapply(DE1, .volcano, pthresh = 1e-5)
lapply(volcanos, function(x) print(x))
dev.off()
names(DE1) = unlist(lapply(DE1, function(x) gsub("_merged","",gsub("_pass","",gsub("_infected","",attr(x,"nme"))))))
names(volcanos) = names(DE1)  
write_xlsx(DE1, paste(resdir, "DE.xlsx",sep="/"))

.shorten<-function(str, len=31){
  str = gsub("_pass","",str)
  str = gsub("_","",str)
  str = gsub("vs","v",str)
  if(nchar(str)>len)str = substr(str,1,len)
  str
}
##DIFF METH  ON COMBINED
depth_combined=.readH5All(transcripts,attributes,filenames, thresh = 100, chrs=NULL, readH5_ = readH5_c)[[1]]
DE2 = lapply(todo, function(x) .processDM(depth_combined,  x[1],x[2], method=.exact, thresh =100))
DE2 = DE2[unlist(lapply(DE2, function(x) dim(x)[[1]]))>0]
names(DE2) = unlist(lapply(DE2, function(x) .shorten(attr(x,"nme"),31)))
volcanos = lapply(DE2, .volcano, pthresh = 1e-5)
pdf(paste(resdir,"DM_combined.pdf",sep="/"))
lapply(volcanos, function(x) print(x))
dev.off()
write_xlsx(.joinSS(lapply(DE2,.xlim,1e-5), sort=T), paste(resdir, "DM_combined.xlsx",sep="/"))


##ALL TRANSCRIPTS AND ALL POSITIONS COMPARED BETWEEN SAMPLES
depth=.readH5All(transcripts,attributes,filenames, thresh = 1000, chrs=NULL, readH5_ = readH5_h)[[1]]
#depth = .transferAttributes(depth, attributes)
DE2 = lapply(todo, function(x) .processDM(depth,  x[1],x[2], method=.exact, thresh =1000))
DE2 = DE2[unlist(lapply(DE2, function(x) dim(x)[[1]]))>0]
names(DE2) = unlist(lapply(DE2, function(x) .shorten(attr(x,"nme"),31)))
volcanos = lapply(DE2, .volcano, pthresh = 1e-5)
pdf(paste(resdir,"DM_1.pdf",sep="/"))
lapply(volcanos, function(x) print(x))
dev.off()
write_xlsx(.joinSS(lapply(DE2,.xlim,1e-5), sort=F), paste(resdir, "DM_1.xlsx",sep="/"))





##LOOK at effect of splicing
depths_combined_spliced=
  lapply(transcripts_all_splice,function(x) .readH5All(x, attributes, filenames, thresh=100, chrs=NULL, readH5_ = readH5_c)[[1]])
depth_spliced = lapply(1:length(info), .mergeDepth, depths_combined_spliced)
names(depth_spliced)= info
nme_spl = names(transcripts_all_splice)
todo2 = lapply((1:length(nme_spl))[-2], function(x)  c(2,x))
names(todo2) = paste(nme_spl[2],nme_spl[-2],sep=" v ")
#depth_spliced[[length(depth_spliced)+1]]= depth_spliced[[2]]+ depth_spliced[[3]]
#names(depth_spliced)[[length(depth_spliced)]] = paste(names(depth_spliced)[2:3],collapse="_")
pdf(paste(resdir,paste("DM_spliced","pdf",sep="."),sep="/"))
DE2_split_all = list()
for(k in 1:length(depth_spliced)){
  print(k)
  depthk = depth_spliced[[k]]
  DE2_split = lapply(todo2, function(x) .processDM(depth_spliced[[k]],  x[1],x[2], method=.exact, thresh =100,plot=F))
  names(DE2_split) = lapply(DE2_split, function(x) .shorten(paste(names(depth_spliced)[k],attr(x,"nme"),sep=":")))
  DE2_split = DE2_split[unlist(lapply(DE2_split, function(x) dim(x)[[1]]))>0]
  volcanos_split = lapply(DE2_split, .volcano, pthresh = 1e-5, prefix =names(depth_spliced)[k] )
  if(length(DE2_split)>0){
    lapply(volcanos_split, function(x) print(x))
    DE2_split_all = c(DE2_split_all, DE2_split)
  }
  #.extractFromDepth(depth_spliced[[k]], 1:5,ORF=NULL, pos=28265)
  
}
dev.off()
write_xlsx(lapply(DE2_split_all,.xlim,1e-5), paste(resdir, "DE_splicing.xlsx",sep="/"))





depths_split=.readH5All(transcripts[1:20,],attributes,filenames, thresh = 100, chrs=NULL, tokeepi = 5, readH5_ = readH5_s)
depth_split = depths_split[[1]]
depth_split1 = abind(depth_combined[,,5,drop=F],depth_split,  along=3)
depth_split1 = .transferAttributes(depth_split1, attributes)

todo1 = lapply(2:length(dimnames(depth_split1)[[3]]), function(x) c(1,x))
names(todo1) = dimnames(depth_split1)[[3]][-1]

DE2_split = lapply(todo1, function(x) .processDM(depth_split1,  x[1],x[2], method=.chisq, thresh =100))
volcanos_split = lapply(DE2_split, .volcano, pthresh = 1e-5)
pdf(paste(resdir,"DM_split.pdf",sep="/"))
lapply(volcanos_split, function(x) print(x))
dev.off()
.extractFromDepth(depth_split1, c(1,grep("^M;end",dimnames(depth_split1)[[3]])),ORF=NULL, pos=28254)


#depth = depths[[1]]
#DE2 = lapply(todo1, function(x) .processDM(depth,  x[1],x[2], method=.chisq, thresh =100))

dev.off()
pdf(paste(resdir,"DM_volcano.pdf",sep="/"))
dev.off()

DE2_=lapply(DE2,.xlim)
names(DE2) = unlist(lapply(DE2, function(x) gsub("_merged","",gsub("_pass","",gsub("_infected","",attr(x,"nme"))))))
names(DE2_) = names(DE2)
write_xlsx(DE2_, paste(resdir, "DM1.xlsx",sep="/"))

attr1 = attributes(depth)

.extractFromDepth(depth, info[todo[[7]]],ORF="st;leader,M;end", pos=26536)

.extractFromDepth(depth, info[todo[[13]]],ORF="", pos=29759)