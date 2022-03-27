
library("rhdf5" ,lib="/usr/local/easybuild-2019/easybuild/software/mpi/gcc/8.3.0/openmpi/3.1.4/r-bundle-bioconductor/3.9-r-3.6.0")
setwd("/data/scratch/projects/punim1068/analysis_genome/resdir4")
library(jsonlite)

## FOLLOWING MERGES SCP data and can be run independently of rest of code
.mergeSCP<-function( scp_dir="../seurat_meta_data/",  out_scp = "scp_meta.txt"){
  if(!file.exists(out_scp)){
    scp_files = grep("scp_meta.txt",dir(scp_dir),v=T)
    names(scp_files) = unlist(lapply(scp_files, function(x) gsub("lib","",gsub("_scp_meta.txt","",x))))
    tab = NULL
    for(i in 1:length(scp_files)){
      x = scp_files[i]
      nme = gsub("lib","",gsub("_scp_meta.txt","",x))
      t = read.table(paste0("../seurat_meta_data/",x), header=T, sep="\t")
      t[,1] = paste(t[,1],nme, sep="-")
      tab = if(is.null(tab)) t else rbind(tab,t)
    }
    write.table(tab,out_scp,sep="\t",quote=F, row.names=F)
  }else{
    print(" scp_meta.txt already exists")
  }
}
.mergeSCP()
##################



args = commandArgs(trailingOnly=TRUE)
if(length(args)>0){
  args = gsub("--","",args)
  argv = (lapply(args, function(x) strsplit(x,"=")[[1]][2]))
  names(argv) = unlist(lapply(args, function(x) strsplit(x,"=")[[1]][1]))
  options(argv)
}


file_=getOption("file","spliced_out")
infile=paste(file_, "h5",sep=".")
h5f = paste(file_, "h5",sep=".")         

#transcript_file = getOption("transcript_file",paste0(file_, ".all_transcripts.txt"))
#all_transcripts = read.table(transcript_file)[,1]
dirs = dir()
dirs = dirs[unlist(lapply(dirs, function(dir){
  file.exists(paste(dir,h5f, sep="/"))
}))]


all_transcripts = sort(unique(unlist(lapply(dirs, function(dir){
  tf = paste0(dir,"/",file_, ".h5.transcripts.mod.txt.gz")
  read.table(tf, header=F)[,1]
}))))

file_out = vector("list", length(dirs))
for(k in 1:length(dirs)){
  dir = dirs[[k]]
  print(dir)
  nme1 = rev(strsplit(strsplit(dir,"\\.")[[1]][1],"_")[[1]])[1]
  nme= gsub("lib","",nme1) 
  h5f_1 = paste(dir, h5f,sep="/")
  tf = paste0(dir,"/",file_, ".h5.transcripts.mod.txt.gz")
  print(h5f_1)
  groups=sort(grep("DATA_TYPES",grep("^/$",unique(h5ls(h5f_1)$group),inv=T,v=T),inv=T,v=T))
  transcripts = read.table(tf, header=F)[,1]
  mi = match(transcripts, all_transcripts)
 all_barcodes = unique(unlist(lapply(groups, function(x){
   gsub("@","",h5read(h5f_1, paste(x,"barcodes",sep="/")))
  }), recursive =T))
 all_barcodes1 = paste(all_barcodes,nme, sep="-")
 print(head(all_barcodes1))
  matr =data.frame(array(0,dim = c(length(all_transcripts), length(all_barcodes)+1)))
                 
  matr[,1] = all_transcripts
 ## fill in the matrix
 for(grp in groups){
   barcodes = gsub("@","",h5read(h5f_1, paste(grp,"barcodes",sep="/")))
   mi_barcodes = match(barcodes, all_barcodes)
   indices = h5read(h5f_1, paste(grp,"indices",sep="/"))
   counts = h5read(h5f_1, paste(grp,"counts",sep="/"))
   for(i in 1:length(barcodes)){
      inds = indices[,i]+1
      cnts = counts[,i]
      inds1 = inds>0
      inds_ = inds[inds1]
      i1 = mi_barcodes[i]
      matr[mi[inds_],i1+1] = cnts[inds1]
   }
 }
  file_out_i =  paste0(dir,"/",file_, "_tex.txt.gz")
  outf = gzfile(file_out_i, open="w")
  #print(file_out_i)
  #outf = gzfile(file_out_i, open="w")
  dimnames(matr)[[2]] = c("GENE",all_barcodes1)
  write.table( matr, outf, quote=F, row.names=F, sep="\t")
  close.connection(outf)
  file_out[[k]] = file_out_i
}



## NOW MERGE THE OUTFILES
#.updateLine<-function(line, k, collapse="\t"){
#  header=unlist(lapply(strsplit(line,"\t")[[1]], function(x) paste(x,k,sep="-")))
#  header[1] = "GENE"
#  paste(header, collapse=collapse)
#}
.subsetLine<-function(line, collapse="\t"){
  header=strsplit(line,"\t")[[1]][-1]
  paste(header, collapse=collapse)
}

.mergeFiles<-function(outfiles, fi_out){
  closeAllConnections()
  outf = gzfile(fi_out, open="w")
  connections = lapply(outfiles, function(fi){
    gzfile( fi,"rb")
  })
  min_length=1
  max_lines = 30000
  i=1
  len = length(connections)
  while(min_length>0 && i<max_lines){
    print(i)
   brk=FALSE
   for(k in 1:length(connections)){
     line=readLines(connections[[k]],1)
     if(length(line)==0){
       brk=TRUE;
       break;
     }
    # if(i==1){
    #   line=.subsetLine(line)
    # }
     if(k>1){
      line=.subsetLine(line) 
    }
     coll = if(k==len) "\n" else "\t"
     writeLines(line, outf,sep=coll)
   }
   if(brk) break;
   i = i+1
 }
closeAllConnections()
}



fi_out=paste0(file_, "_tex.txt.gz")

.mergeFiles(file_out, fi_out)

params = read_json("params.json")
params$input_file = fi_out
params$prefix=file_
params$outdir = paste(params$outdir,file_,sep="_")
params_file = paste0(file_,".json")

write_json(params, params_file, pretty=T, auto_unbox=T)



