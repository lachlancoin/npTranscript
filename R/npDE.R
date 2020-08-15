install = Sys.getenv(c("INSTALL_R_PCKGS", "INSTALL_UPDATE"))
libs_to_install = unlist(strsplit(Sys.getenv("R_PCKGS_TO_INSTALL"),","))
if(install[1]=="TRUE"){
  if(length(libs_to_install)==0){
      libs_to_install = c("VGAM","ggplot2","biomaRt","edgeR","writexl","ggrepel")
  }
   update = install[2]=="TRUE"
    install.packages("BiocManager", update=F,ask=F)
    for(i in libs_to_insall){
      BiocManager::install(i, update=update, ask=F)
    }
} 

options("np.qThresh"=2)

options("np.control"="control")
options("np.case"='infected')
options("np.exclude"="none")
options("np.analysis"="betabinom")
options("np.depth_thresh"="100")
options("np.prefix_keep"=NULL)
options("np.prefix_remove"="^[0-9]{1,}\\.")
options("np.prefix_sequins"="^R[0-9]_")
options('np.isoformDepth'=100)
options("np.isoform.test"="chisq.test") #fisher.test
options("np.dm.test"="chisq.test") #fisher.test
options("np.maxIsoformGroups"=5)
options("np.adjustMethod"="BH");

#options("np.source"="../../R")
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0){
  argv = (lapply(args, function(x) strsplit(x,"=")[[1]][2]))
  names(argv) = unlist(lapply(args, function(x) strsplit(x,"=")[[1]][1]))
  options(argv)
}

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



target= list( chrom="character", 
             ORFs ="character",start = "numeric", 
             end="numeric", ID="character" ,type_nme="character", countTotal="numeric")


##READ TRANSCRIPT DATA
infilesT = grep("transcripts.txt", dir(), v=T)
transcriptsl = readTranscriptHostAll(infilesT, start_text = start_text,target = target,   filter = filter, 
                                    combined_depth_thresh = getOption("np.depth_thresh",100))

attributes = attributes(transcriptsl)
filenames = attr(transcriptsl,"info")
control_names = unlist(lapply(control_names, grep, filenames, v=T));#  grep(control_names,filenames,v=T)
infected_names = unlist(lapply(infected_names, grep, filenames, v=T))
exclude_nme = getOption("np.exclude",'none')
for(i in 1:length(exclude_nme)){
  control_names = grep(exclude_nme[i], control_names, inv=T,v=T)
  infected_names = grep(exclude_nme[i], infected_names, inv=T,v=T)
}

#print(paste(type_names[1],paste(control_names,collapse=" ")))
#print(paste(type_names[2] , paste(infected_names, collapse=" ")))
if(length(control_names)!=length(infected_names)) error(" lengths different")
transcriptsl = lapply(transcriptsl, .processTranscripts)

if(!is.null(getOption("np.prefix_remove",NULL))){
  remove = lapply(transcriptsl,.filter, getOption("np.prefix_remove"), inv=F)
  keep = lapply(transcriptsl,.filter, getOption("np.prefix_remove"), inv=T)
}else {
  remove = lapply(transcriptsl,.filter, getOption("np.prefix_keep","\\."), inv=T)
  keep = lapply(transcriptsl,.filter, getOption("np.prefix_keep","\\."), inv=F)
}

.rename<-function(x) {
  strand = rep("",length(x))
  strand[regexpr("\\+$",x$ORFs)>=0]="+"
  strand[regexpr("\\-$",x$ORFs)>=0]="-"
  apply(cbind(x$chrs, x$start, x$end, strand),1,paste,collapse="_")
}
##CHANGE NAME OF NOVEL ORFS TO BE MORE READABLE
for(i in 1:length(remove)) remove[[i]]$ORFs = .rename(remove[[i]])

if(!is.null(getOption("np.prefix_sequins",NULL))){
  sequins = lapply(remove,.filter, getOption("np.prefix_sequins"), inv=F)
  remove = lapply(remove,.filter, getOption("np.prefix_sequins"), inv=T)
}else{
  sequins = list();
}

transcripts_all = list(genes=keep, novel=remove, sequins = sequins)
transcripts_all =  lapply(transcripts_all, .process, control_names, infected_names)

pdf(paste(resdir, "/DE.pdf",sep=""))
DE_list = lapply(transcripts_all, .processDE, attributes, resdir, control_names, infected_names, type_names, type="known",plot=T)
for(i in 1:length(res_)) attr(DE_list[[i]],"nme") = names(transcripts_all)[i]

##OUTPUT FILE
h5DE = "DE.h5"
file.remove(h5DE)
h5createFile(h5DE)
h5createGroup(h5DE,"DE")
lapply(DE_list, function(x) h5write(x, h5DE,paste("DE",attr(x, "nme"),sep="/")))
#h5ls(h5DE)
#


write_xlsx(lapply(DE_list,function(x) x[attr(x,"order"),,drop=F]), paste(resdir, "DE.xlsx",sep="/"))
#.volcano(res_keep$DE1, pthresh = 1e-3, prefix="keep")


volcanos = lapply(DE_list, .volcano, pthresh = 1e-5)
lapply(volcanos, function(x) print(x))
dev.off()


print("####DEPTH ANALYSIS #### ")
##isoform analysis
pdf(paste(resdir, "/isoDE.pdf",sep=""))
infilesAltT = grep("isoforms.h5" , dir(),v=T)
pvs_all = lapply(transcripts_all, .testIsoformsAll, infilesAltT,n=getOption("np.maxIsoformGroups",5), test_func = getOption("np.isoform.test","chisq.test"))
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
h5closeAll()

write_xlsx(lapply(DE2,.xlim,1e-2), paste(resdir, "DM_combined.xlsx",sep="/"))

if(FALSE){
fisher.test(.extractFromDepth(depth_combined[[1]], 1:2,"12.56159663")[1:2,,])
}

