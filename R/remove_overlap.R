
library(rhdf5)
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0) args = c("out4.txt", "out4a.txt", "out4b.txt", "depth4")

src = c("~/github/npTranscript/R", "../../R")

#SHOULD BE RUN IN data/ subdirectory
.findFile<-function(path, file, exact = T){
  for(i in 1:length(path)){
    if(exact){
      res = paste(path[i],file,sep="/") 
      if(file.exists(res)) { 
        return(res)
      }
    }else{
      files = grep(file, dir(path[i]) , v=T)
      #	if(length(files)==1){
      return(paste(path[i], files,sep="/"))
      #	}
    } 
  }
}

if(FALSE) src = gsub("~","C:/Users/LCOIN", src)
source(.findFile(src, "transcript_functions.R"))

####FUNCTIONS #####
.overl<-function(v1,v2)  min(v2[2]-v1[1], v1[2] - v2[1])
.overlaps<-function(row, b, thresh = 0, chrom_strand_ind = length(row)){
  inds1 = which(b[,chrom_strand_ind] == as.character(row[ chrom_strand_ind]))
  if(length(inds1)==0) return (c())
  b1 = apply(b[inds1,3:4],c(1,2), as.numeric)
  overl= apply(b1,1,.overl,as.numeric(row[3:4]))
  inds2 = which(overl>thresh)
  if(length(inds2)==0) return(c())
  inds1[inds2]
}
.readT<-function(infile){
  table1 = read.table(infile, skip=2, sep="\t", head=F, as.is=T,fill=T)
  header = names(read.table(infile, nrows=3, sep="\t", head=T, as.is=T, fill=T))
  names(table1) = header[1:dim(table1)[2]]
  chrom_index = unlist(lapply(table1$ORFs,function(x) strsplit(x,"\\.")[[1]][1]))
  strand = unlist(lapply(table1$ORFs,function(x) rev(strsplit(x,"")[[1]])[1]))
  chrom_strand = as.character(apply(cbind(chrom_index, strand),1,paste, collapse="_"))
  res = data.frame(cbind(table1, chrom_strand))
  res$chrom_strand = as.character(res$chrom_strand)
  res
}
#returns indices of transcripts with overlapping exons and filtercount >0
.overlappingTranscripts<-function(isoentry, v, filter="depth4"){
  if(is.array(isoentry$breaks)) breaks = data.frame(isoentry$breaks) else breaks = isoentry$breaks
  cnts = isoentry$counts
  inds_iso = if(is.null(filter)) 1:(dim(cnts)[1]) else which(cnts[,grep(filter,names(cnts))]>0)
  if(length(inds_iso)==0){
    return (c())
  }else{
    br=breaks[inds_iso] 
    overlaps = lapply(br,.overlapExon,v)
  }
  inds_iso[which(overlaps>0)]
}
.overlapExon<-function(br1, v){
  indsv = seq(1, length(br1),2)
  o = unlist(lapply(indsv, function(j)  .overl(br1[j:(j+1)],v)))
  max(o)
}
##################


table1 = .readT(args[1])
table2 = .readT(args[2])
filter=args[4] 
#filter="depth4"

chrom_strand_ind = grep("chrom_strand", names(table1))

if(FALSE){
  ###dEBUG
  test = .overlaps(table1[1,],table1,thresh=0, chrom_strand_ind = chrom_strand_ind)
  print(table1[test,1:5])
}

res = apply(table1,1, .overlaps,table2,thresh =0, chrom_strand_ind = chrom_strand_ind)
reslen = lapply(res, length)
inds1 = which(reslen>0)
table1a = table1[which(reslen==0),]
table1b = table1[which(reslen>0),]
res = res[reslen>0]
chrom_index = as.factor(table1b$chrom_strand)
chromlev = levels(chrom_index)
tokeep = rep(F, dim(table1b)[1])
for(i in 1:length(chromlev)){
  inds_i = which(chrom_index==chromlev[i])
  isofile =paste(strsplit(chromlev[i],"_")[[1]][1],"isoforms.h5",sep=".")
  if(file.exists(isofile)){
  for(j in 1:length(inds_i) ){
    row = table1b[inds_i[j],]
    # table1b[inds_i[j],3:4]
   print(isofile)
    iso2 = try(readIsoformH5(isofile,  table2[res[[inds_i[j]]],]))
   if(!inherits(iso2,"try-error")) {
    H5close()
    if(!is.null(iso2)){
      inds_iso2 = lapply(iso2,.overlappingTranscripts, v=row[3:4],filter=filter)
      tokeep[inds_i[j]] = length(which(unlist(lapply(inds_iso2, length))>0))==0
    }else{
	#print(isofile)
    print( table2[res[[inds_i[j]]],1])
      print("did not find isoform model, need to rerun npTranscript with lower threshold")
      tokeep[inds_i[j]]=FALSE
    }
  }else{
	  print(paste(isofile,"error"))
  }
  }
  }else{
    print(paste(isofile," does not exist"))
  }
}

print(paste("rescued ",table1b[which(tokeep),]$ID))



table1a = rbind(table1a, table1b[tokeep,,drop=F])

write.table(table1a, file = args[3] , quote=F, sep="\t", row.names=F, col.names=F)
