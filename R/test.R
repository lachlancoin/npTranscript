library(ggplot2)
library(gridExtra)
library(RColorBrewer)



path = "C:/Users/LCOIN/bitbucket/rapids_rna_seq/corona" 
source(paste(path,"test_functs.R",sep="/"))

dir_ = "bins"
files1 = grep(".txt",dir(dir_),v=T)
files = paste(dir_,files1,sep="/")
resdir = "res"
dir.create(resdir)
#depth_thresh = 10
#frac_thresh = 0.4
#depth_thresh1 = 10



tables = readTables(files, files1)
results = list()
kmers = getKmerOffsets(4)
for(j in 1:length(kmers)){
	results[[j]] = processMods(tables,  depth_thresh = 0, frac_thresh =0.4, depth_thresh1 = 10, skip=grep("bin", files1, inv=T), kmer =kmers[[j]])

}
names(results) = names(kmers)


avgError = getAvg(results[[1]]$big_t)




ks = 5
offset = 4
base = "C"

motifAnalysis(results, ks, offset, base)


result = processMods(files, files1, depth_thresh = 0, frac_thresh =0.4, depth_thresh1 = 10, skip=grep("bin", files1, inv=T), kmer = seq(-2,2,1))
or = order(result$kmer_ratios[,1], decreasing = T)
or1 = order(result$kmer_abs[,1], decreasing = T)
result$kmer_ratios[head(or,20),]
result$kmer_abs[head(or1,20),]

head(result$big_t[which(result$kmer=="CGCGG"),col_inds])

result1 = processMods(files, files1, depth_thresh = 0, frac_thresh =0.4, depth_thresh1 = 10, skip=grep("bin", files1, inv=F), kmer = kmer)

small_t = result$big_t[result$all_pos[,1],,drop=F]

#small_t1 = cbind(small_t[,c(1,2)], apply(small_t[,-c(1,2)], c(1,2), function(s) sprintf("%5.3g", s)))

#names(all_mod) = nmes
write.table(small_t, paste1(resdir,"mods.csv"), sep=",", quote=F, col.names=T, row.names=F)


small_t1 = result1$big_t[result1$all_pos[,1],,drop=F]
write.table(small_t1, paste1(resdir,"mods1.csv"), sep=",", quote=F, col.names=T, row.names=F)


ml0 = plotFAll(result$long_t, c(0,100), exclude="allreads") 
ml1 = plotFAll(result$long_t, c(0,30000), exclude="bin") 
ml2= plotFAll(result$long_t, c(20000,30000), exclude="allreads") 

ml3= plotFAll(result$long_t, c(29040,29050), exclude="allreads") 
ml3



try(ggsave( paste1(resdir,"modsStart.pdf"), plot=ml0, width = 30, height = 30, units = "cm"))
try(ggsave( paste1(resdir,"modsAll.pdf"), plot=ml1, width = 30, height = 30, units = "cm"))
try(ggsave( paste1(resdir,"modsEnd.pdf"), plot=ml2, width = 30, height = 30, units = "cm"))
	




####TOMBO ANALYSIS
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


tom_kmer = list()
tom_kmer_lev = list()
for(i in 1:length(results)){	
	tom_kmer[[i]] = results[[i]]$kmer[pos_ref]
	tom_kmer_lev[[i]]= getlev(tom_kmer[[i]])
}

names(tom_kmer) = names(results)
names(tom_kmer_lev) = names(results)
lapply(tom_kmer_lev, head, 1)
kmerPerc = lapply(results, getKmerPerc)
tomb_ratios = list()
baseTombo = "C"
out_res = list()
for(i in 1:length(kmerPerc)){
	offset = which(as.numeric(strsplit(names(results)[i],"_")[[1]])==0)-1
	rownmes = data.frame(strsplit(names(kmerPerc[[i]]), ""))
	ba = as.character(unlist(rownmes[offset+1,]))
	kmerPerc[[i]][which(ba != baseTombo)] =0
	kmerPerc[[i]] = kmerPerc[[i]]/sum(kmerPerc[[i]])
 	tr = matchKmerPerc(kmerPerc[[i]], tom_kmer_lev[[i]])
 	tomb_ratios[[i]] = tr[order(tr[,4], decreasing=T),]
	 out_res[[i]] = head(tomb_ratios[[i]],2)
 
}
names(out_res) = names(results)
names(tomb_ratios) = names(results)






