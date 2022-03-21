#library(jsonlite)
#setwd( "/data/scratch/projects/punim1068/analysis_DTU/resdir/resdir1")
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0){
  args = gsub("--","",args)
  argv = (lapply(args, function(x) strsplit(x,"=")[[1]][2]))
  names(argv) = unlist(lapply(args, function(x) strsplit(x,"=")[[1]][1]))
  options(argv)
}


split=getOption("split","\t")
prefix=getOption("prefix","merged_")
collapse=getOption("collapse","\t")
input_file <-getOption("input_matrix","chrom_out1_tex.txt")
max_lines=30000
dirs= grep(prefix, dir(),v=T)
names(dirs)=dirs

closeAllConnections()

connections = lapply(dirs, function(dir_){
   gzfile( paste(dir_,input_file,sep="/"),"rb")
})


.updateLine<-function(line, k){
  header=unlist(lapply(strsplit(line,"\t")[[1]], function(x) paste(x,k,sep="-")))
  header[1] = "GENE"
  paste(header, collapse=collapse)
}

.subsetLine<-function(line){
  header=strsplit(line,"\t")[[1]][-1]
  paste(header, collapse=collapse)
}



fi="chrom_out_tex.txt"
outf = file(fi, open="w")
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
    if(i==1){
      line=.updateLine(line,k)
    }
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

tab_all = NULL
for(i in 1:length(dirs)){
  tab = read.table(paste(dirs[i],"scp_meta.txt", sep="/"), header=T)
  tab[,1] = paste(tab[,1],"-",i,sep="")
  #tab[1,1] = "NAME"
  tab_all = if(is.null(tab_all)) tab else rbind(tab_all, tab)
}
write.table(tab_all,"scp_meta.txt",sep="\t", row.names=F, quote=F)
