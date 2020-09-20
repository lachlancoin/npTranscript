options("np.install"="FALSE")
options("np.libs_to_install"="rhdf5,VGAM,ggplot2,writexl,ggrepel,abind");
options("np.results"="results")
options("np.control"="control")
options("np.case"='infected')
#options("np.comparisons"=NULL)

options("np.casecontrol" = "casecontrol.txt")
options("np.exclude"="none")
options("np.virus"="FALSE")
options("np.analysis"="betabinom")
options("np.combined_depth_thresh"="100")
options("np.depth_thresh"="1000")
options("np.prefix_keep"=NULL)
options("np.prefix_remove"="^[0-9]{1,}\\.")
options("np.prefix_sequins"="^R[0-9]_")
options('np.isoformDepth'=100)
options("np.dm.test"="chisq.test") #fisher.test
options("np.maxIsoformGroups"=5)
options("np.adjustMethod"="BH");
##specific options for running on windows machine
#options("np.source"="../../R")
#options("np.libdir"="C:/Users/LCOIN/R-4.0.2/library")

args = commandArgs(trailingOnly=TRUE)
if(length(args)>0){
  args = gsub("--","",args)
  argv = (lapply(args, function(x) strsplit(x,"=")[[1]][2]))
  names(argv) = unlist(lapply(args, function(x) strsplit(x,"=")[[1]][1]))
  options(argv)
}


print(.libPaths())
vers = R.Version()
version=paste(vers$major,vers$minor,sep=".")
np.libdir = getOption("np.libdir", "~/R/lib")
libdir=paste(np.libdir,version,sep="/")
dir.create(libdir , recursive=T)
.libPaths(libdir)
.libPaths()
#options("np.source"="../../R")

libs_to_install = unlist(strsplit(getOption("np.libs_to_install"),","))
if(getOption("np.install","FALSE")=="TRUE"){
  install.packages("BiocManager", lib=libdir, repos="https://cran.ms.unimelb.edu.au/")
  library("BiocManager")
  for(i in libs_to_install){
    BiocManager::install(i, update=F, ask=F, lib=libdir)
  }
} 

annotFile = "annotation.csv.gz"
control_names = unlist(strsplit(getOption("np.control"),':'))
infected_names = unlist(strsplit(getOption("np.case"),':'))
type_names = c(control_names[1], infected_names[1])

if(getOption("np.analysis","betabinom")=="edgeR"){
		library(edgeR)
	}else{
	  library(VGAM)
}
library(stats)
library(writexl)
library(rhdf5)
library(ggplot2)
library(ggrepel)
library(abind)

resdir = getOption("np.results","results")
dir.create(resdir);


.findFile<-function(path, file, exact = T){
  for(i in 1:length(path)){
    if(exact){
      res = paste(path[i],file,sep="/") 
      if(file.exists(res)) { 
        return(res)
      }
    }else{
      files = grep(file, dir(path[i]) , v=T)
      if(length(files)==1){
        return(paste(path[i], files,sep="/"))
      }
    } 
  }
}


source(.findFile(getOption("np.source","~/github/npTranscript/R"), "transcript_functions.R"))
source(.findFile(getOption("np.source","~/github/npTranscript/R"), "diff_expr_functs.R"))
optionFile = paste(resdir,"options.txt",sep="/")
.printOptions(optionFile)


files = dir()
chroms = NULL
isVirus=  getOption("np.virus","FALSE")=="TRUE"
start_text = "start"
filter = NULL
#if(isVirus){
#  start_text ="endBreak"

#}

.rename<-function(x) {
  strand = rep("",length(x))
  strand[regexpr("\\+$",x$ORFs)>=0]="+"
  strand[regexpr("\\-$",x$ORFs)>=0]="-"
  apply(cbind(x$chrs, x$start, x$end, strand),1,paste,collapse="_")
}


target= list( chrom="character",
             span ="character",ORFs="character",start = "numeric", 
             end="numeric", ID="character" ,type_nme="character", countTotal="numeric")



##READ TRANSCRIPT DATA

infilesT = grep("transcripts.txt", dir(), v=T)
combined_depth_thresh = as.numeric(getOption("np.combined_depth_thresh","100"))
transcriptsl_unmerged = try(.readTranscriptsHost(infilesT,target=target,filter = filter, 
                                                 combined_depth_thresh = combined_depth_thresh))

transcriptsl = transcriptsl_unmerged
attributes = attributes(transcriptsl)
filenames = attr(transcriptsl,"info")
if(!isVirus){
jointind = grep(';', transcriptsl$span)
transcriptsl$span[jointind] = unlist(lapply(transcriptsl$span[jointind],function(x) strsplit(x,";")[[1]][1]))
nullind = grep('^-$',transcriptsl$span)
transcriptsl$span[nullind] = .rename(transcriptsl[nullind,])
transcriptsl =.mergeRows(transcriptsl,sum_names = filenames,append_names = "ID", colid="span")
transcriptsl = transcriptsl[,-which(names(transcriptsl)=="ORFs")]
names(transcriptsl) = sub("span","ORFs",names(transcriptsl))
}
control_names = unlist(lapply(control_names, grep, filenames, v=T));#  grep(control_names,filenames,v=T)
infected_names = unlist(lapply(infected_names, grep, filenames, v=T))
exclude_nme = getOption("np.exclude",'none')
for(i in 1:length(exclude_nme)){
  control_names = grep(exclude_nme[i], control_names, inv=T,v=T)
  infected_names = grep(exclude_nme[i], infected_names, inv=T,v=T)
}

if(!is.null(getOption("np.casecontrol",NULL))){
  casecontrol = .readCaseControl(getOption("np.casecontrol",NULL))  
  if(!is.null(casecontrol)){
  control_names = as.character(casecontrol$control)
  infected_names =as.character( casecontrol$case)
  }
}



if(length(control_names)!=length(infected_names)) error(" lengths different")
transcriptsl =  .processTranscripts(transcriptsl)

  remove_inds =  regexpr("^[0-9]{1,}\\.", as.character(transcriptsl$ORFs))>=0
  keep_inds = !remove_inds
  sequins_inds = regexpr(getOption("np.prefix_sequins"), transcriptsl$ORFs)>=0
  grp = rep(NA, dim(transcriptsl)[[1]])
  grp[!remove_inds]="genes"
  grp[remove_inds] = "unknown"
  grp[sequins_inds]="sequins"
  grp = as.factor(grp)
#transcripts_all1 = lapply(transcripts_all, .mergeRows,sum_names= c(control_names, infected_names), colid="geneID" )
  dimnames(transcriptsl)[[1]] = transcriptsl$ID
# 
#info = attr(transcripts,'info')
  if(!isVirus){
  transcriptsl1 =.addAnnotation(annotFile, transcriptsl,grp,  colid="geneID", nmes = c("chr" , "ID" , "Name" , "Description","biotype"))
}else{
  transcriptsl1 = cbind(transcriptsl,grp)
}

  DE_list =  .processDE(transcriptsl1, attributes, resdir, control_names, infected_names, type_names, type="")
  
#  DE2 = .transferAttributes(DE2, attributes)
  
 
##OUTPUT FILE
h5DE = paste(resdir,"DE.h5",sep="/")
file.remove(h5DE)
h5createFile(h5DE)
h5write(transcriptsl1, h5DE, "transcripts")
h5createGroup(h5DE,"DE")
for(i in 1:length(DE_list)){
  grpnme = paste("DE",names(DE_list)[i], sep="/")
  h5write(DE_list[[i]], h5DE, grpnme)
}
prefix = "" #paste(control_names[1], infected_names[1], sep=" v " )
towrite = lapply(DE_list,.xlim, pthresh = 1.0, col = "p.adj")
#towrite[[length(towrite)+1]] = transcriptsl1
#names(towrite)[length(towrite)] = "transcripts"
.write_xlsx1(towrite,paste(resdir, "DE.xlsx",sep="/") )
DE2 = data.frame(unlist(DE_list,recursive=FALSE))
#write_xlsx(lapply(list(transcripts=transcriptsl1,DE=DE2),function(x) x[attr(x,"order"),,drop=F]), paste(resdir, "DE.xlsx",sep="/"))
volcanos = lapply(DE_list, .volcano, logFCthresh = 0.5, prefix=prefix, top=10, exclude=c())
todo=.getAllPairwiseComparisons(names(DE_list), start=2)
comparisonPlots = lapply(todo, function(x) .comparisonPlot(DE2,transcriptsl1 , inds  = x, excl=c()))

pdf(paste(resdir, "/DE.pdf",sep=""))
if(TRUE){
  for(i in 1:length(DE_list)).qqplot(DE_list[[i]]$pvals, min.p= 1e-200,main="" ,add=i>1, col=i+1)
}
lapply(volcanos, function(x) print(x))
lapply(comparisonPlots, function(x) print(x[[1]]))
lapply(comparisonPlots, function(x) print(x[[2]]))


dev.off()


print("####DEPTH ANALYSIS #### ")
##isoform analysis
isofile = grep("isoforms.h5" , dir(),v=T)
isoforms_i=  try(readIsoformH5(transcriptsl1, isofile[1],depth=getOption("np.depth_thresh",1000)))
if(inherits(isoforms_i,"try-error")) {
  print(paste(" could not read h5 file", isofile[1]))
}else{
  inds_i = attr(isoforms_i,"inds")
  pvs_all =  .testIsoformsAll(transcriptsl1, isoforms_i,filenames, control_names, infected_names, n=getOption("np.maxIsoformGroups",5), test_func =chisq.test)
  pvs_all = pvs_all[unlist(lapply(pvs_all, length))>0]
  pvs_all2 = data.frame(unlist(pvs_all,recursive=FALSE))
  
  for(i in 1:length(pvs_all)) attr(pvs_all[[i]],"nme") = names(pvs_all)[i]
  
  volcanos =  lapply(pvs_all, .volcano, top=10, useadj = F, prefix =prefix,  logFCthresh = 0.1)
  todo=.getAllPairwiseComparisons(names(pvs_all), start=2)
  comparisonPlots = lapply(todo, function(x) .comparisonPlot(pvs_all2,transcriptsl1[inds_i,] , inds  = x, excl=c()))
  
  pdf(paste(resdir, "/isoDE.pdf",sep=""))
  for(i in 1:length(pvs_all)) .qqplot1(pvs_all[[i]],"pvals", main = names(pvs_all)[[i]], add = i>1)
  lapply(volcanos, function(x) print(x))
  lapply(comparisonPlots, function(x) print(x[[1]]))
  lapply(comparisonPlots, function(x) print(x[[2]]))
  dev.off()

  h5createGroup(h5DE,"isoDE")
  lapply(pvs_all, function(x) h5write(x, h5DE,paste("isoDE",attr(x, "nme"),sep="/")))
  h5ls(h5DE)
  write_xlsx(pvs_all, paste(resdir, "isoDE.xlsx",sep="/"))
}
#h5closeAll()

######
if(file.exists("0.clusters.h5")){
transcripts_ = transcriptsl_unmerged[transcriptsl_unmerged$countTotal>1000,,drop=F]
if(isVirus){
depth_combined =  try(readH5_c( "0.clusters.h5", transcripts_,  filenames, thresh = 1000))
}else{
  depth_combined=  try(readH5_h("0.clusters.h5",transcripts_,filenames, thresh = 1000))
}
if(inherits(depth_combined,"try-error")) {
 print(" skipping error analysis ") 
}else{
  if(!isVirus){
    depth_combined = .mergeDepthByPos(depth_combined)
  }
  DE2 =.processDM(depth_combined, filenames, control_names ,infected_names, method=getOption("np.dm.test", "chisq.test"), thresh_min =100,plot=T,adjust="none")
  
  volcanos = lapply(DE2, .volcano, top=10,prefix=prefix, logFCthresh = 0.1)
  todo=.getAllPairwiseComparisons(names(DE2), start=2)
  DE3 = data.frame(unlist(DE2,recursive=FALSE))
  
  comparisonPlots = lapply(todo, function(x) .comparisonPlot(DE3,transcripts_ , inds  = x, excl=c()))
  
  
  pdf(paste(resdir,"DM_combined.pdf",sep="/"))
  
  lapply(volcanos, function(x) print(x))
  lapply(comparisonPlots, function(x) print(x[[1]]))
  lapply(comparisonPlots, function(x) print(x[[2]]))
  dev.off()
  
  h5createGroup(h5DE,"DM")
  lapply(DE2, function(x) h5write(x, h5DE,paste("DM",attr(x, "nme"),sep="/")))
  h5ls(h5DE)
  H5close()
  #h5closeAll()
  ord=order(unlist(lapply(names(DE3),function(x) strsplit(x,"\\.")[[1]][2])))
  write_xlsx(lapply(DE2,.xlim,1e-2), paste(resdir, "DM_combined.xlsx",sep="/"))
#write_xlsx(lapply(list(DE3=DE3[,ord]),.xlim,pthresh=1e-2,col="meta.p.adj" ), paste(resdir, "DM_combined1.xlsx",sep="/"))
}
}

if(FALSE){
  .extractFromDepth(depth_combined,1:8,"28782")[3,,]
  .extractFromDepth(depth_combined, 5:6,"MT.1047.0")
  .extractFromDepth(depth_combined1, 5:6,"MT.2684")
  
  
  .extractFromDepth(depth_combined, 5:6,"19.48966785")
fisher.test(.extractFromDepth(depth_combined, 1:2,"12.56159663")[1:2,,])
}

