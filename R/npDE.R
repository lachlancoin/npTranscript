
options("np.install"="FALSE")
options("np.libs_to_install"="")
options("np.results"="results")
options("np.control"="control")
options("np.case"='infected')
options("np.exclude"="none")
options("np.analysis"="betabinom")
options("np.depth_thresh"="100")
options("np.prefix_keep"=NULL)
options("np.prefix_remove"="^[0-9]{1,}\\.")
options("np.prefix_sequins"="^R[0-9]_")
options('np.isoformDepth'=100)
options("np.dm.test"="chisq.test") #fisher.test
options("np.maxIsoformGroups"=5)
options("np.adjustMethod"="BH");
print(.libPaths())
vers = R.Version()
version=paste(vers$major,vers$minor,sep=".")
libdir=paste("~/R/lib",version,sep="/")
dir.create(libdir , recursive=T)
.libPaths(libdir)
.libPaths()
#options("np.source"="../../R")
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0){
  argv = (lapply(args, function(x) strsplit(x,"=")[[1]][2]))
  names(argv) = unlist(lapply(args, function(x) strsplit(x,"=")[[1]][1]))
  options(argv)
}
libs_to_install = unlist(strsplit(getOption("np.libs_to_install"),","))
if(getOption("np.install","FALSE")=="TRUE"){
  if(length(libs_to_install)==0){
    libs_to_install = c("rhdf5","VGAM","ggplot2","biomaRt","edgeR","writexl","ggrepel","abind")
  }
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
binsize = 1e5
isVirus=  FALSE;
start_text = "start"
mergeByPosAndGene = FALSE
filter = NULL
if(isVirus){
  filter = list("type_nme"="5_3")
  mergeByPosAndGene = TRUE
  binsize = 100  #1e5 for euk
  start_text ="endBreak"

}

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
transcriptsl_unmerged = readTranscriptHostAll(infilesT, start_text = start_text,target = target,   filter = filter, 
                                    combined_depth_thresh = getOption("np.depth_thresh",100))[[1]]
transcriptsl = transcriptsl_unmerged
jointind = grep(';', transcriptsl$span)
transcriptsl$span[jointind] = unlist(lapply(transcriptsl$span[jointind],function(x) strsplit(x,";")[[1]][1]))
nullind = grep('^-$',transcriptsl$span)
transcriptsl$span[nullind] = .rename(transcriptsl[nullind,])
attributes = attributes(transcriptsl)
filenames = attr(transcriptsl,"info")
transcriptsl_unm = transcriptsl
transcriptsl =.mergeRows(transcriptsl,sum_names = filenames, colid="span")
transcriptsl = transcriptsl[,-which(names(transcriptsl)=="ORFs")]
names(transcriptsl) = sub("span","ORFs",names(transcriptsl))
control_names = unlist(lapply(control_names, grep, filenames, v=T));#  grep(control_names,filenames,v=T)
infected_names = unlist(lapply(infected_names, grep, filenames, v=T))
exclude_nme = getOption("np.exclude",'none')
for(i in 1:length(exclude_nme)){
  control_names = grep(exclude_nme[i], control_names, inv=T,v=T)
  infected_names = grep(exclude_nme[i], infected_names, inv=T,v=T)
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
 transcriptsl1 =.addAnnotation(annotFile, transcriptsl,grp,  colid="geneID", nmes = c("chr" , "ID" , "Name" , "Description","biotype"))


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
DE2 = data.frame(unlist(DE_list,recursive=FALSE))
write_xlsx(lapply(list(transcripts=transcriptsl1,DE=DE2),function(x) x[attr(x,"order"),,drop=F]), paste(resdir, "DE.xlsx",sep="/"))
volcanos = lapply(DE_list, .volcano, logFCthresh = 0.5, top=20, exclude=c())
todo=.getAllPairwiseComparisons(names(DE_list), start=2)
comparisonPlots = lapply(todo, function(x) .comparisonPlot(DE2,transcriptsl1 , inds  = x, excl=c()))

pdf(paste(resdir, "/DE1.pdf",sep=""))
if(FALSE){
  for(i in 1:length(DE1)).qqplot(DE1[[i]]$p.adj, min.p= 1e-200,main=type,add=T, col=i+1)
}
lapply(volcanos, function(x) print(x))
lapply(comparisonPlots, function(x) print(x[[1]]))
lapply(comparisonPlots, function(x) print(x[[2]]))


dev.off()


print("####DEPTH ANALYSIS #### ")
##isoform analysis
pdf(paste(resdir, "/isoDE.pdf",sep=""))
isofile = grep("isoforms.h5" , dir(),v=T)
isoforms_i=  readIsoformH5(transcriptsl_unmerged, isofile[1],depth=getOption("np.depth_thresh",100))
pvs_all =  .testIsoformsAll(isoforms_i,n=getOption("np.maxIsoformGroups",5), test_func =chisq.test)
pvs_all = pvs_all[unlist(lapply(pvs_all, length))>0]
for(i in 1:length(pvs_all)) attr(pvs_all[[i]],"nme") = names(pvs_all)[i]
for(i in 1:length(pvs_all)) .qqplot1(pvs_all[[i]],"p", main = names(pvs_all)[[i]], add = i>1)

volcanos = lapply(pvs_all, .volcano, pthresh = 1e-5, useadj = T)
lapply(volcanos, function(x) print(x))

dev.off()
h5createGroup(h5DE,"isoDE")
lapply(pvs_all, function(x) h5write(x, h5DE,paste("isoDE",attr(x, "nme"),sep="/")))
h5ls(h5DE)
write_xlsx(lapply(pvs_all,function(x) x[order(x$p),,drop=F]), paste(resdir, "isoDE.xlsx",sep="/"))

#h5closeAll()

######
pdf(paste(resdir,"DM_combined.pdf",sep="/"))

depth_combined=lapply(transcripts_all,
                      function(transcripts_)   .readH5All(transcripts_,attributes,filenames, thresh = 1000,  readH5_ = readH5_h))

DE2 = lapply(depth_combined,.processDM,  1,2, method=getOption("np.dm.test", "chisq.test"), thresh =1000,plot=T,adjust="none")
DE2 = DE2[unlist(lapply(DE2, function(x) !is.null(x) && dim(x)[[1]]))>0]
for(i in 1:length(DE2)) attr(DE2[[i]],"nme") = names(DE2)[i]

volcanos = lapply(DE2, .volcano, pthresh = 1e-4)

lapply(volcanos, function(x) print(x))
dev.off()

h5createGroup(h5DE,"DM")
lapply(DE2, function(x) h5write(x, h5DE,paste("DM",attr(x, "nme"),sep="/")))
h5ls(h5DE)
H5close()
#h5closeAll()

write_xlsx(lapply(DE2,.xlim,1e-2), paste(resdir, "DM_combined.xlsx",sep="/"))

if(FALSE){
fisher.test(.extractFromDepth(depth_combined[[1]], 1:2,"12.56159663")[1:2,,])
}

