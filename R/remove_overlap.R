.overl<-function(v1,v2)  min(v2[2]-v1[1], v1[2] - v2[1])

overlaps<-function(row, b, thresh = 0){
	inds1 = which(b[,2] ==as.character(row[2]))
	if(length(inds1)==0) return (0)
	b1 = apply(b[inds1,3:4],c(1,2), as.numeric)
	overl= apply(b1,1,.overl,as.numeric(row[3:4]))
	return(length(which(overl>thresh)))
}

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0) args = c("out4.txt", "out4a.txt", "out4b.txt")
header = names(read.table(args[1], nrows=3, sep="\t", head=T, as.is=T, fill=T))
table1 = read.table(args[1], skip=2, sep="\t", head=F, as.is=T)
table2 = read.table(args[2], skip=2, sep="\t", head=F, as.is=T)
names(table1) = header[1:dim(table1)[2]]
names(table2) =header[1:dim(table1)[2]]

res = apply(table1,1, overlaps,table2)

table1 = table1[res==0,]
write.table(table1, file = args[3] , quote=F, sep="\t", row.names=F, col.names=F)
