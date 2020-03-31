library(ggplot2)
library(gridExtra)

source("test_functs.R")



dir_ = "data"
files1 = grep(".txt",dir(dir_),v=T)
files = paste(dir_,files1,sep="/")
resdir = "res"
dir.create(resdir)
#depth_thresh = 10
#frac_thresh = 0.4
#depth_thresh1 = 10



tables = readTables(files, files1)
results = list()
kmers = getKmerOffsets(3)
for(j in 1:length(kmers)){
	results[[j]] = processMods(tables,  depth_thresh = 0, frac_thresh =0.4, depth_thresh1 = 10, skip=grep("bin", files1, inv=T), kmer =kmers[[j]])
}
names(results) = names(kmers)



##TOMBO ANALYSIS
ref = results[[1]]$kmer
tombo = read.table("data/Subgenome_5mC.csv", sep=",", head =T)
ref_pos = grep("Position_ref", names(tombo))
indsToFix = tombo[,ref_pos]>71
tombo[indsToFix,ref_pos] = tombo[indsToFix,ref_pos] + 60
tom1 = tombo[which(tombo$Fraction > 0.9),]
tom2 = tombo[which(tombo$Fraction > 0.5),]
pos_ref = unique(sort(tom1$Position_reference))
pos_ref2 = unique(sort(tom2$Position_reference))


kmer10 = getKmer(ref, pos_ref, v= -10:10)
writeFasta(kmer10, "kmer10.fasta")




getCountBin<-function(tombo_, pos_ref_, bins = levels(as.factor(tombo_$Subgenome))){
 tab = array(NA, dim = c(length(pos_ref_), ncol = length(bins)), dimnames = list(pos_ref_, bins))
 for(i in 1:length(pos_ref_)){
 	bi = tombo_[which(tombo_$Position_reference==pos_ref_[i]),2:3]
	 tab[i,match(bi$Subgenome, bins) ] = bi[,1]
 }
 cbind(pos_ref,tab)
}


matr =  getCountBin(tombo, pos_ref)

write.table(matr, "mods_table.csv", sep=",", quote=F, col.names=T, row.names=F)

