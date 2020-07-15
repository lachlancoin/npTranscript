install= FALSE
if(install){
  install.packages("BiocManager")
}

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  args[1] = "out2.txt"
  args[2] = "out3a.txt"
  args[3] = "out3b.txt"
}

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


src = c( "../../R" , "~/github/npTranscript/R" )
source(.findFile(src, "transcript_functions.R"))
source(.findFile(src, "diff_expr_functs.R"))

t =read.table(args[1], head=T,as.is=T)
sum_names = c(grep("count",names(t),v=T), grep("depth" , names(t),v=T), grep("errors" , names(t),v=T))
for(i in 1:length(sum_names)){
  t[,sum_names[i]] = as.numeric(t[,sum_names[i]])
  t[is.na(t[,sum_names[i]]),sum_names[i]]=0
}
t2 = .mergeRows(t,sum_names = sum_names,  colid='span')
write.table(t2,args[2],col.names=TRUE, sep="\t", quote=F, row.names=F)

#t3 = t2[t2$count4==0 & t2$count5<5 & t2$countTotal >=5 ,]
#t3_1 = t2[t2$count4==0 & t2$countTotal >=5,]
#t4 = t3[order(t3$countTotal, decreasing=T),]
#t4_1 = t3_1[order(t3_1$countTotal, decreasing=T),]
#write.table(t4,args[2],col.names=TRUE, sep="\t", quote=F, row.names=F)
#write.table(t4_1,paste(args[2],"1.txt",sep="."),col.names=TRUE, sep="\t", quote=F, row.names=F)

