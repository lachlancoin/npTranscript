.metapv<-function(pvi){
  nonNA = which(!is.na(pvi))
  if(length(nonNA)==0) return(NA)
  pchisq(sum(unlist(lapply(pvi[nonNA],qchisq,lower.tail=F, df=1 ))),df=length(nonNA), lower.tail=F)
  
}
readH5_h<-function(h5file, df, thresh =100,log=F){
  IDS = df$ID
  header = h5read(h5file,"header")
  
  if(header[2]!="base") header = c(header[1], "base", header[-1])
  depth_inds = grep("depth[0-9]", header)
  error_inds = grep("errors[0-9]", header)
  dinds  = grep("depth", header)
  einds  = grep("error", header)
  
  countT = apply(df[,grep("count[0-9]", names(df)),drop=F],1,min)
  pos_ind = which(header=="pos")
  names = h5ls(h5file)$name
  
  inds = which(countT>thresh & df$ID %in% names)
  if(length(inds)==0) return(data.frame(matrix(nrow = 0, ncol = length(header)+3)))
  
  mat_all = NULL; 
  #
  for(i in 1:length(inds)){
    print(paste(i, length(inds), sep=" of "))
    ID = as.character(IDS[inds[i]])
    mat = t(h5read(h5file,as.character(ID)))
    countT1 = apply(mat[,depth_inds,drop=F],1,sum)
    indsi = which(countT1>thresh)
    if(length(indsi)>0){
    mat1 = mat[indsi,,drop=F]
  #  v = unlist(mat1[163,])
  #  betaBinomialP2(v,depth_inds, error_inds, control=1, case=2,binom=F, log=F)
    pv1 = t(apply(mat1,1,betaBinomialP2, depth_inds,error_inds , control=1, case=2, binom=F,log=log))
    #dimnames(pv1) = list(mat$pos[indsi],c("pv1","pv2"))
    IDv = rep(ID, length(indsi))
  #  chrv=rep(df$chrs[i], length(indsi))
    mat2 = data.frame(IDv ,mat1,pv1)
    #print(head(mat2))
    if(dim(mat2)[1]>0){
      if(is.null(mat_all)) mat_all = mat2 else mat_all = rbind(mat_all,mat2)
    }
    }
  } 
  if(is.null(mat_all)) mat_all = data.frame(matrix(nrow = 0, ncol = length(header)+3))
  names(mat_all) = c("IDS",header, "pv1","pv2")
  #print(head(mat_all[mat_all$pv2<1e-100,]))
  mat_all
}

DE_err<-function(DE, inds_control=1, inds_case=2, sum_thresh = 100){
  
  inds2 = grep("error_ratio[0-9]", names(DE))
  inds1 = grep("count[0-9]", names(DE))
  if(length(inds2)>0){
    m = max(DE[,inds2,drop=F],na.rm=T)
    if(m<=1.01){
    for(j in 1:length(inds2)){
      DE[,inds2[j]] =  DE[,inds2[j]] * DE[,inds1[j]] 
    } 
    }
  }
  
  countT = apply(DE[,grep("sum_", names(DE)),drop=F],1,sum)
#  print(head(sort(countT, decreasing=T),20))
  inds_ = which(countT>sum_thresh)
 # print(head(inds_),20)
  if(length(inds_)==0) return (NULL)
  inds2 = grep("error_ratio[0-9]", names(DE))
  inds1 = grep("count[0-9]", names(DE))
  DE_a = DE[,c(inds1,inds2),drop=F]
  inds1 = 1:length(inds1)
  inds2 = (1:length(inds2))+length(inds1)
  pv_res = t(apply(DE_a[inds_,,drop=F],1,betaBinomialP2, inds1,inds2, 1, 2))
  dimnames(pv_res) =  list(inds_, c("pv1","pv2"))
  pv_res
}

#betaBinomialP2(v,depth_inds, error_inds, control=1, case=2,binom=F, log=F)
betaBinomialP2<-function(v, indsDepth, indsError, control, case,binom=F,log=F){
  pv1 = rep(NA, length(control))
  pv2 = rep(NA, length(control))
  
  for(j in 1:length(control)){
    indsj = c(indsDepth[control[j]], indsError[control[j]], indsDepth[case[j]], indsError[case[j]])
    indsj2 = c(indsDepth[case[j]], indsError[case[j]], indsDepth[control[j]], indsError[control[j]])
    pv1[j] = betaBinomialP1(v,indsj , binom=binom, lower.tail=T, log=log)
    pv2[j] = betaBinomialP1(v, indsj, binom=binom, lower.tail=F, log=log)
  }
  pv1m =chisqCombine(pv1,log=log)
  pv2m = chisqCombine(pv2,log=log)
  c(pv1m,pv2)
#  2*min(pv1m, pv2m, na.rm=T)
}



betaBinomialP1<-function(v, ord, binom=F, lower.tail=T,log=F){
  if(length(which(is.na(v)))!=0) return (NA)
  size = v[ord[1]];
  shape1 = v[ord[2]]
  y = shape1
  shape2 = size-shape1
  sizex = v[ord[3]]
  shape1x = v[ord[4]]
  x = shape1x
  
  shape2x = sizex -shape1x
  if(sizex==0 || size==0) return (1)
  x1 = x;
  if(lower.tail==F)x1 = x-1 #so its greater than or equal to 
  if(binom){
    pv_res = pbinom(x1,size = sizex,prob = y/size,log=log)
  }else{
    pv_res = pbetabinom.ab(x1,size = sizex,shape1 = shape1,shape2 =shape2,log=log)
    
  }
  if(log){
    if(lower.tail==F) pv_res =   log(max(0,1-exp(pv_res)))
  }else{
    if(lower.tail==F) pv_res =   max(0,1-pv_res)
  }
  pv_res
}
  

betaBinomialP<-function(x,y, binom=F, lower.tail=T,log=F){
  size = ceiling(sum(y))
  # geneID = df$geneID
  
  shape1 = y 
  shape2 = size -y 
  proby = y/size;
  zeros = y==0 & x== 0
  sizex = ceiling(sum(x))
  shape1x = x 
  shape2x = sizex -x 
  probx = x/sizex
  
  pbb =rep(NA, length(x))
  if(binom){
    if(lower.tail==FALSE){
      pbb[!zeros] = pbinom(x[!zeros]-1,size = sizex,prob = proby[!zeros],log.p=log)
    }else{
      pbb[!zeros] = pbinom(x[!zeros],size = sizex,prob = proby[!zeros],log.p=log)
      
    }
  }else{
    if(lower.tail==FALSE){
      pbb[!zeros] = pbetabinom.ab(x[!zeros]-1,size = sizex,shape1 = shape1[!zeros],shape2 =shape2[!zeros],log=log)
    }else{
      pbb[!zeros] = pbetabinom.ab(x[!zeros],size = sizex,shape1 = shape1[!zeros],shape2 =shape2[!zeros],log=log)
    }
  }
  if(log){
    if(lower.tail==F) pbb[!zeros] = log(1-exp(pbb[!zeros]))
    pbb[!zeros] =  pbb[!zeros]/log(10)
  }
else{
  if(lower.tail==F) pbb[!zeros] = 1-pbb[!zeros]
  pbb[pbb>1] = 1
  pbb[pbb<0] = 0
} 
  pbb
}

chisqCombine<-function(pv,log=log){
  #if(TRUE) return(pv[1])
  nonNA = which(!is.na(pv))
  
  if(length(nonNA)==0) return(NaN) else if(length(nonNA)==1) return (pv[nonNA]);
  if(log){
    resp = pchisq(sum(unlist(lapply(exp(pv[nonNA]),qchisq,lower.tail=F, df=1 ))),df=length(nonNA), lower.tail=F,log=T)
    
  }else{
  resp = pchisq(sum(unlist(lapply(pv[nonNA],qchisq,lower.tail=F, df=1 ))),df=length(nonNA), lower.tail=F)
  }
}

.write<-function(DE1, resdir, filename="results.csv"){
  DE1[,1:9] = apply(DE1[,1:9], c(1,2), function(x) sub(' ' , '', sprintf("%5.3g",x)))
  write.table(DE1[attr(DE1,"order"),],file=paste(resdir,filename,sep="/") , quote=F, row.names=F, sep="\t", col.names=T)
}

#which x is significiantly more or less than expected given y
#if(lower.tail=T returns p(x<=y) else p(x>=y)
##ASSUMES MATCHED DATA BETWEEN CONTROL  AND INFECTED
DEgenes<-function(df,control_names,infected_names, edgeR = F,  type="lt", binom=F, log=F,
                  remove=c(control_names, infected_names, "countT", "ID","type_nme", "ORFs")
                  ){
  lower.tail = T
  control_inds = rep(NA, length(control_names))
 infected_inds = rep(NA, length(infected_names))
  for(i  in 1:length(control_inds)){
	control_inds[i] = which(names(df)==control_names[i])
infected_inds[i] = which(names(df)==infected_names[i])
 }
  if(!edgeR){
    pvalsM1 = matrix(NA,nrow = dim(df)[1], ncol = length(control_inds))
    pvalsM2 = matrix(NA,nrow = dim(df)[1], ncol = length(control_inds))
    for(i in 1:length(control_inds)){
      print(i)
        x = df[,control_inds[i]]
        y = df[,infected_inds[i]]
        pvalsM1[,i] = betaBinomialP(x,y, binom=binom, lower.tail=lower.tail,log=log)
        pvalsM2[,i] = betaBinomialP(y,x, binom=binom, lower.tail=lower.tail,log=log)
        
    }
    pvals1 = apply(pvalsM1, 1, chisqCombine,log=log)
    pvals2 = apply(pvalsM2, 1, chisqCombine,log=log)
   # pvals = 2*apply( cbind(pvals1,pvals2),1,min)
  #  lessThan = pvals2<pvals1
  }else{
    qlf = DE_egdeR(df, control_inds, infected_inds)
    pvals1 = qlf$table$P
    pval2 = pvals1
#    lessThan = qlf$coefficients[,2]<0
  }
 pvals = apply(cbind(pvals1,pvals2),1,min, na.rm=T)
 p.adj = p.adjust(pvals,method="BH")
 p.adj1 = p.adjust(pvals1, method="BH");
 p.adj2 = p.adjust(pvals2, method="BH");
 
 x = apply(df[,control_inds,drop=F],1,sum)
 y = apply(df[,infected_inds,drop=F],1,sum)
 
  tpm_control = (x/sum(x))*1e6
  tpm_infected = (y/sum(y))*1e6
  probX1 = (x+0.5)/sum(x+.5)
  probY1 = (y+0.5)/sum(y+.5)
  ratio1 = probY1/probX1
  logFC = log(ratio1)/log(2)
  output =  data.frame(pvals, p.adj, pvals1,pvals2,p.adj1, p.adj2, tpm_control, tpm_infected, ratio1,logFC, sum_control=x,sum_infected=y)
  output= cbind(output, df[, -which(names(df) %in% remove)])
  
  
  
 
    attr(output, "order") = order(pvals)
    attr(output, "order1") = order(pvals1)
    attr(output, "order2") = order(pvals2)
    
  
  #print('h')
  att = grep('class' ,grep('names', names(attributes(df)), inv=T, v = T),inv=T,v=T)
  if(length(att)>0){
  for(i in 1:length(att)){
    
    attr(output,att[i]) = attr(df,att[i])
    
  }
}
  output
#  output[orders[,1],,drop=F]
}

.transferAttributes<-function(output, attributes){
 att = grep('class' ,grep('names', names(attributes), inv=T, v = T),inv=T,v=T)
 attributes1 = attributes[which(names(attributes) %in% att)]
 att = names(attributes1)
 print(att)
if(length(att)>0){
  for(i in 1:length(att)){
    
   attr(output,att[i]) = attributes1[[i]]
    
  }
}
 output
}


DE_egdeR<-function(df, inds_control, inds_infected){
  groups = c(rep(1,length(inds_control)),rep(2,length(inds_infected)))
  y <- DGEList(counts=df[,c(inds_control,inds_infected)],group=groups)
  y <- calcNormFactors(y, fdr_thresh = 0.1)
  design <- model.matrix(~groups)
  y <- estimateDisp(y,design)
  if(is.na(y$common.dispersion)) return (rep(NA, dim(df)[1]))
  #To perform quasi-likelihood F-tests:
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
 qlf
}

getDescr<-function(DE,mart, thresh = 1e-10, prefix="ENSCS"){
  inds = which(DE$FDR1<thresh | DE$FDR2<thresh)
  print(length(inds))
  subDE = DE[inds,,drop=F]
  genenames1 = as.character(subDE$geneID)
  
  subinds = grep(prefix, genenames1)
  desc1 = rep("", dim(DE)[1])
  
  if(length(subinds)>0){
    genenames = genenames1[subinds];
    attr =  c('ensembl_gene_id','description')# 'go_id') #, "name_1006", "namespace_1003") #, "definition_1006")
    filt = c('ensembl_gene_id')
    #	goids = getBM(attributes =attr,   filters = filt,   values = list(ensg) ,    mart = mart) 
   print(genenames)
     desc =  biomaRt::getBM(attributes=attr, filters = filt, mart = mart, values = list(genenames)) 
    #FDR = DE$FDR
    desc1[inds][subinds] = desc[match(genenames, desc[,1]),2]
    
  }
  data.frame(cbind(DE,desc1))
}



.getDescEnrich<-function(DE2,mart, nme="pvals1",nme2="pos1M",thresh = 1e-5, go_thresh = 1e-2){
  sigChrGT = findSigChrom(DE2,thresh = thresh, go_thresh = go_thresh,  nme=nme,nme2=nme2)
  if(is.null(mart)) return(sigChrGT)
  desci = list()
  for(j in 1:(dim(sigChrGT)[1])){
    genesInChrom=findGenesByChrom(DE2,chrom=as.character(sigChrGT$chrs[j]), nme2=nme2, nme=nme, type="FDR", thresh = 1e-5)
    desci[[j]] = getDescr(genesInChrom, mart,thresh = 1e-5, prefix=prefix)
    #print(desci)
  }
  names(desci) = sigChrGT$chrs
  desci
}



getChromIDs<-function(ensg, mart){
  goids = getBM(attributes = c('ensembl_gene_id', 'chromosome_name'),   filters = c('ensembl_gene_id'),   values = list(ensg) ,    mart = mart) 
  goids = goids[goids[,2]!="",]
  lev_all = getlev(goids$chromosome)
  #dimnames(lev_all)[[1]] = goids[match(lev_all[,1], goids$chromosome),3]
  
  chromObj = list(goids=goids, lev_all = lev_all)
  chromObj
}

getGoIDs<-function(genenames, mart){
  ensg = genenames
  #genenames = exons_$GENENAME.1
  attr =  c('ensembl_gene_id', 'go_id', "name_1006", "namespace_1003") #, "definition_1006")
  filt = c('ensembl_gene_id')
  #	goids = getBM(attributes =attr,   filters = filt,   values = list(ensg) ,    mart = mart) 
  
  goids =  biomaRt::getBM(attributes=attr, filters = filt, mart = mart, values = list(genenames)) 
  goids2 = goids[goids[,2]!="",]
  
  gn = genenames[match(goids2[,1], ensg)]
  goids = cbind(goids2, gn)
  lev_all = getlev(goids$go_id)
  lev_all = cbind(lev_all,goids[match(lev_all[,1], goids$go_id),3])
  names(lev_all)[3] = "description"
  goObjs = list()
  goObjs[[1]] = list(goids=goids, lev_all = lev_all )
  ns = as.factor(goids$namespace)
  lev = levels(ns)
  
  for(i in 1:length(lev)){
    goids1 = goids[goids$namespace==lev[i],]
    lev_all1 = getlev(goids1$go_id)
     lev_all1 = cbind(lev_all1,goids1[match(lev_all1[,1], goids1$go_id),3])
     names(lev_all1)[3] = "description"
    goObjs[[i+1]] = list(goids = goids1, lev_all = lev_all1)
  }
  lev[lev==""]  = "blank"
  names(goObjs) = c("combined", lev)
  goObjs = goObjs[names(goObjs)!="blank"]
  goObjs
}

findGenesByChrom<-function(DE,chrom="MT", thresh = 1e-10,nme2="chrs", nme="FDR1", type="FDR"){
  col_ind = grep(nme, names(DE))[i]
  chr_ind = grep(nme2,names(DE))[1]
  inds = which(DE[,chr_ind]== chrom & DE[,col_ind]<thresh)
  print(inds)
  if(length(inds)==0) return (NULL)
  DE[inds,,drop=F]
  
}

.readFeatureCounts<-function(files){
  a = read.table(files[1], head=T)
  b = read.table(files[2], head=T)
  geneID = as.character(a$Geneid)
  chrs = unlist(lapply(a$Chr,function(x) strsplit(as.character(x),";")[[1]][1]))
  df = data.frame(control = a[,7],infected = b[,7],chrs = chrs, geneID = geneID)
  df
}
.processTranscripts<-function(transcript, prefix="ENS"){
  geneID= as.character(unlist(lapply(strsplit(transcript$ORFs,";"), function(v) v[1])))
  rightGene = as.character(unlist(lapply(strsplit(transcript$ORFs,";"), function(v) v[length(v)])))
indsL =  grep(prefix,geneID,inv=T)
indsR =  grep(prefix,rightGene)
comb = which(indsL %in% indsR)
if(length(comb)>0) geneID[comb] = rightGene[comb]
  cbind(transcript, geneID)
}

.addAnnotation<-function(annotfile, transcripts, colid="geneID",nmes = c("ID" , "Name" , "Description","biotype")){
  match_ind = which(names(transcripts)==colid)[1]
  gfft = read.table("annotation.csv.gz", sep="\t", header=F, fill=T, quote='\"')
  names(gfft) = nmes
  
  #gfft[,1] = gsub("transcript:", "", as.character(gfft[,1]))
  gfft = gfft[match(as.character(transcripts[,match_ind]), as.character(gfft[,1])),]
  gfft[,1] = transcripts$ID
  transcripts = cbind(gfft,transcripts)
  return(transcripts)
}

.mergeRows<-function(transcripts1,sum_names = c(),  colid='geneID'){
  sum_names = unique(c(sum_names,"countTotal"));
  ind = which(names(transcripts1) %in% colid)
  ind_s = which(names(transcripts1) %in% sum_names)
  lev = getlev(transcripts1[,ind])
  
  todo = as.character(lev[lev$cnts>1,]$lev)
  if(length(todo)==0) return(transcripts1)
  extract_inds = lapply(todo, function(x) which(transcripts1[,ind]==x))
  torem = sort(unlist(extract_inds))
  subt = transcripts1[unlist(lapply(extract_inds,function(x) x[1])),]
  for(i in 1:length(extract_inds)){
    subind = extract_inds[[i]]
    subt[i,ind_s] =apply(transcripts1[subind,ind_s,drop=F],2,sum)
  }
  rbind(transcripts1[-torem,,drop=F], subt)
  
}

readTranscriptHostAll<-function(infilesT, 
                                combined_depth_thresh = 100,
                                start_text = "start", 
                                filter  = NULL,
                                target= list(chrom="character", 
                                             leftGene="character", rightGene="character", start = "numeric", 
                                             end="numeric", ID="character", isoforms="numeric" ,error_ratio0 = "numeric",error_ratio1="numeric") ){
  chroms = unlist(lapply(infilesT, function(x) strsplit(x,"\\.")[[1]][[1]]))
  chrom_names = rep("", length(chroms))
  dfs = list()
  ncol = -1
  for(i in 1:length(chroms)){
    infilesT1 = paste(chroms[i],"transcripts.txt.gz", sep=".")
   # print(infilesT1)

    dfi = try(.readTranscriptsHost(infilesT1,target=target,filter = filter, combined_depth_thresh = combined_depth_thresh, start = start_text))
if(inherits(dfi,"try-error")) {
	print(infilesT1)
	dfs[[i]] = NULL
}else{
    dfi = dfi[order(dfi$start),]
    if(dim(dfi)[1]>0){
      chrom_names[i] = dfi$chrs[1]
     # print(chrom_n)
      dfs[[i]] = dfi
	#print(dim(dfi))
      ncol = dim(dfi)[[2]]
      #print(ncol)
      #print(' not null')
    }else{
      print("is null")
      dfs[[i]] = NULL
    }
  }
  }
  lengs = unlist(lapply(dfs,function(x) if(is.null(x)) NA else dim(x)[1]))
  inds = which(!is.na(lengs))
  chroms = chroms[inds]
  chrom_names = chrom_names[inds]
  dfs = dfs[inds]
  #print(chrom_names)
  numeric_names = as.numeric(chrom_names)
 # print(numeric_names)
  ord = c()
  i1 = which(!is.na(numeric_names))
  i2 = which(is.na(numeric_names))
  if(length(i1)>0){
    ord1 = order(numeric_names[i1])
    ord = c(ord,i1[ord1])
  }
  if(length(i2)>0){
    ord2 = order(chrom_names[i2],decreasing=T)
    ord = c(ord,i2[ord2] )
  }
  
  dfs = dfs[ord]
  chroms = chroms[ord]
  chrom_names = chrom_names[ord]
  
  attr(dfs,"info")=attr(dfs[[1]],"info")
  attr(dfs,"chroms")=chroms
  attr(dfs,"chrom_names")=chrom_names
  
  dfs
}
.combineTranscripts<-function(dfs, chrom=attr(df, "chroms"), chrom_names =attr(dfs, "chrom_names") ){
  lengs = unlist(lapply(dfs,function(x) if(is.null(x)) NA else dim(x)[1]))
  ncol = dim(dfs[[1]])[2]
  #lengs = lengs[inds][ord]
  res = data.frame(matrix(nrow  =sum(lengs), ncol = ncol))
  names(res) = names(dfs[[1]])
  start=1;
  ranges = matrix(nrow = length(dfs), ncol=3)
  for(i in 1:length(dfs)){
   # print(i)
    lengi = dim(dfs[[i]])[1]
    for(k in 1:ncol){
      if(is.factor(dfs[[i]][,k])) dfs[[i]][,k] = as.character(dfs[[i]][,k])
    }
    ranges[i,] = c(start, start+lengi-1,dfs[[i]]$chrs[1]);
    res[ranges[i,1]:ranges[i,2],] = dfs[[i]]
    start = start+lengi
  }
  names(chroms) = chrom_names
  m1 =rep(NA, dim(res)[1])# matrix(nrow = dim(depth)[1], ncol=2)
  for(i in 1:(dim(ranges)[1])){
    ri = as.numeric(ranges[i,]) 
    m1[ri[1]:ri[2]]=ri[3]
  }
  attr(res,"chr_inds")=m1
  attr(res,"ranges") = ranges
  attr(res,"chroms") = chroms
  attr(res,"chrom_names") = chrom_names
  attr(res,"info") = attr(dfs,"info")
 # attr(res,"chroms")=chrom_
  dimnames(res)[[1]] = res$ID
  
  res
}


readIsoformH5<-function(h5file,  transcripts_){
  .ext<-function(x) x[unique(c(0,which(x>0)))]
  header = h5read(h5file,"header")
  IDS = transcripts_$ID
  trans = list()
  for(i in 1:length(IDS)){
    ID = IDS[i]
    mat = t(h5read(h5file,as.character(ID)))
    cnts = data.frame(mat[,1:length(header)])
    names(cnts) = header
    dimnames(cnts) = list(cnts$subID, header)
    cnts = cnts[,-1]
    cntT = apply(cnts,1,sum)
    ord = order(cntT, decreasing=T)
    cnts = cnts[ord,,drop=F]
    mat = mat[ord,-(1:length(header))]
    
    transi = apply(mat,1,.ext)
    names(transi) =dimnames(cnts)[[1]] #paste(ID,cnts$subID,sep='.')
    trans[[i]] = list(breaks = transi, counts = cnts)
  } 
  names(trans) = IDS
  trans
  
}
.readH5All<-function(transcripts, attr = attributes(transcripts), thresh = 1000,chroms= attr[which(names(attr)=="chroms")][[1]]){
  depth = NULL
  ranges = matrix(nrow = length(ord), ncol=2)
  start = 1
  depths = list()
  for(i in 1:length(chroms)){
    print(chroms[i])
    infile = paste(chroms[i],"clusters.h5", sep=".")
  
    dfi= try(readH5_h(infile, transcripts, thresh =thresh,log=F))
    if(inherits(dfi,"try-error")) {
depths[[i]] = NULL
	}else{
depths[[i]]  = dfi
	}

  }
H5close();
 # if(TRUE) return(depths[[1]])
  
  lengs = unlist(lapply(depths,function(x) if(is.null(x)) NA else dim(x)[1]))
  inds = which(!is.na(lengs))
  chroms = chroms[inds]
  depths = depths[inds]
  lengs = lengs[inds]
  res = data.frame(matrix(nrow  =sum(lengs), ncol = dim(depths[[1]])[2]))
  names(res) = names(depths[[1]])
  start=1;
  
  ranges = data.frame(matrix(nrow = length(depths), ncol=4))
  names(ranges) = c("start","end","chrom","chrom_name");
  if(length(depths)>0){
  for(i in 1:length(depths)){
    #print(i)
    depths[[i]][,1] = as.character(depths[[i]][,1])
    
    lengi = dim(depths[[i]])[1]
    ranges[i,] = c(start, start+lengi-1, chroms[i], names(chroms)[i]);
    
    if(lengi>0){
     res[ranges[i,1]:ranges[i,2],] = depths[[i]]
     start = start+lengi
    }
  }
}
  attr(res,"ranges") = ranges
  attr(res,"chroms")=chroms
  m1 =rep(NA, dim(res)[1])# matrix(nrow = dim(depth)[1], ncol=2)
  for(i in 1:(dim(ranges)[1])){
   ri = as.numeric(ranges[i,]) 
   m1[ri[1]:ri[2]]=ri[3]
  }
  attr(res,"chr_inds")=m1
  res
  }


.readTranscriptsHost<-function(infilesT1, 
                               filter = NULL,
                  target= list(chrom="character", leftGene="character", rightGene="character", start = "numeric", end="numeric", ID="character")
              ,prefix="ENSC" ,combined_depth_thresh =100  , start_text = "start"                                
  ){

  names(target)[names(target)=="start"] = start_text
  
  header = names(read.table( infilesT1,sep="\t", head=T, nrows = 3, comment.char='#'))
   extra = grep("count[0-9]", header,v=T)
extran = as.list(rep("numeric", length(extra)))
names(extran) = extra
target = c(target,extran)
  #print(header)
  inf = scan(infilesT1, nlines=1, what=character())
  #if(length(grep("#", inf))>0) attr(transcripts,"info") = sub("#", "",inf)
  
  inf = sub('#','',inf)
  types = unlist(lapply(inf, function(x) rev(strsplit(x,"_")[[1]])[1]))
  header_inds = match(names(target),header)
  #print(target)
  colClasses = rep(NA, length(header));
  colClasses[header_inds] = target
  #colClass = cbind(rep("numeric", length(extra)), colClasses)
  transcripts = read.table( infilesT1,sep="\t", head=T, comment.char='#', colClasses= colClasses)
  if(!is.null(filter)){
    for(k in 1:length(filter)){
      nme_ind = grep(names(filter)[k], names(transcripts))
      if(length(nme_ind)==0) stop(paste("not found ",names(filter)[k]))
      indsk=which(transcripts[,nme_ind]==filter[k])
      transcripts = transcripts[indsk,,drop=F]
    }
  }
  names(transcripts)[names(transcripts)==start_text] = "start"
  names(target)[names(target)==start_text] = "start"
  #print(head(transcripts))
  
  
  header_inds1 = match(names(target),names(transcripts))
 # head_inds1 = grep("count[0-9]", names(transcripts));
countT = as.numeric(transcripts$countTotal)
 # countT = apply(transcripts[,head_inds1,drop=F],1,function(x) sum(as.numeric(x)))
 #print(countT)
 
indsk = countT>=combined_depth_thresh
  transcripts = transcripts [indsk,header_inds1, drop=F] 
  names(transcripts)  = sub("chrom", "chrs" ,names(transcripts))
  if(length(grep("#", inf))>0) inf = sub("#", "",inf)
  attr(transcripts,"types")=types
  attr(transcripts,"info")=inf
  names(transcripts)[grep("count[0-9]",names(transcripts))] = inf
  transcripts
}

.appendGeneNamesToDepth<-function(depth, transcripts, sort= "pv1"){
  gene_names=apply(transcripts[match(depth$IDS,
                                     transcripts$ID),names(transcripts) %in% c("chrs","geneID","rightGene")],1,paste,collapse=".")
  
  #d_s = depth[,c(grep("depth", names(depth)), grep("errors", names(depth)))]
  #apply(d_s,1 ,function(v) c(v[3]/v[1], v[4]/v[2]))
  depth1 = cbind(gene_names, depth)
  if(is.null(sort)) return (depth1)
  ord = order(depth[,grep(sort,names(depth))] )
  
  depth1[ord,]
}


findGenes<-function(goid, goObj,DE, fdr_thresh = 1e-10, lessThan = FALSE){
  
  inds =  which(goObj$goids$go_id==goid) 
  genes = goObj$goids[inds,,drop=F]
  inds1 = which(DE$geneID %in% genes$ensembl_gene_id)
  if(!is.null(lessThan)){
    inds2 = which(DE[inds1,]$FDR<fdr_thresh & DE[inds1,]$lessThan==lessThan)
  }else{
    inds2 = which(DE[inds1,]$FDR<fdr_thresh)
    
  }
 # ge = DE$geneID[inds1[inds2]]
#  print(genes[which(genes$ensembl_gene_id %in% ge),])
  DE[inds1[inds2],]
}


getGoGenes<-function(go_categories,goObjs, lessThan = T, fdr_thresh = 1e-5){
  names(go_categories) = lapply(go_categories,function(x, goObj) as.character(goObj$lev_all[which(goObj$lev_all[,1]==x),3]), goObjs[[1]])
  go_genes1 = lapply(go_categories, findGenes, goObjs[[1]],DE1, fdr_thresh = fdr_thresh, lessThan=lessThan)
  names(go_genes1) = names(go_categories)
  go_genes1 = go_genes1[which(unlist(lapply(go_genes1,function(x) dim(x)[1]))>0)]
  go_genes1
}

.compareLevs<-function(lev_, lev_all, k =   sum(as.numeric(as.character(lev_[,2]))), 
mn = sum(as.numeric(as.character(lev_all[,2])))){
print(paste(k,mn))
	 inds_m = match(as.character(lev_[,1]), as.character(lev_all[,1]))
  	lev_ = cbind(lev_,lev_all[inds_m,1:2,drop=F])
  	lev_1 = t(apply(lev_,1,.phyper2,  k = k, mn = mn))
	res = data.frame(cbind(lev_, lev_1))
	res$pv = as.numeric(as.character(res$pv))
	res$enrich = as.numeric(as.character(res$enrich))
		res$enrich99 = as.numeric(as.character(res$enrich99))
	res[order(res$pv),]
}



.getKComp<-function(fa, fh, inds = 1:(dim(fh)[1]), k =8){
	l1 = floor(k/2)
	l2 = k -l1
	kComp = list("stLeft"= getKmer(fasta[[1]], unique(fh$startPos[inds]), v =(-k:0)),
	   "stRight"= getKmer(fasta[[1]], unique(fh$startPos[inds]), v =(0:k)),
		"endLeft"= getKmer(fasta[[1]], unique(fh$endPos[inds]), v =(-k:0)),
		"endRight" = getKmer(fasta[[1]], unique(fh$endPos[inds]), v =(0:k)),
		"st" = getKmer(fasta[[1]], unique(fh$startPos[inds]), v =-l1:l2),
		"end" = getKmer(fasta[[1]], unique(fh$endPos[inds]), v =-l1:l2)
	)
res = lapply(kComp, getlev)
names(res) = names(kComp)
res
}

findSigGo_<-function(goObj, DE1_1, fdr_thresh = 1e-10, go_thresh = 1e-5, prefix="ENSC", nme="FDR1"){
  #DE = DE1[DE1$lessThan==lessThan,,drop=F]
  goids = goObj$goids 
  ensg_inds = grep(prefix, DE1_1$geneID)
  DE1 =DE1_1[ensg_inds,]
  ensg =DE1$geneID
  goidx = rep(FALSE,length(goids$ensembl_gene_id ))
  
  
  pvs = DE1[,grep(nme,names(DE1))]
  sig =  which(pvs<fdr_thresh )
  goidx1 = goids$ensembl_gene_id %in% ensg[sig]
  a = data.frame(goidx1)
  lev_all = goObj$lev_all
  
  suma = apply(a,1,sum)
  subs = which(suma>0)
  go1 = goids[subs,]
  lev1 = getlev(go1[,2])  #go_ids or chromosome_name
  go_todo = lev1[,1]
  go_ = goids[goidx1,]
  
  lev_ = getlev(go_[,2], todo=go_todo)
  inds_m = match(as.character(lev_[,1]), as.character(lev_all[,1]))
  lev_ = cbind(lev_,lev_all[inds_m,1:2,drop=F])
  
  lev_1 = t(apply(lev_,1,.phyper2,  k = length(sig), mn = length(ensg)))
  #	
  if(dim(goids)[[2]]>2){
    descr = goids[match(lev_[,1], goids[,2]),3]
#      apply( cbind(lev_[,1],,1,paste,collapse=":")
  goids = as.character(lev_[,1])
   lev_1 = cbind(goids, lev_1, descr)
  #
  # lev_1[,pv_ind]=sprintf("%5.3g", lev_1[,pv_ind])
   #lev_1$pv = sprintf("%5.3g", as.numeric(as.character(lev_1$pv)))
  }else{
    
    #dimnames(lev_1)[[1]] = lev_[,1]
  }
  pv_ind = which(dimnames(lev_1)[[2]]=="pv")
  pvs = as.numeric(lev_1[,pv_ind])
  len = length(which(pvs<go_thresh))
  #lev_1[,pv_ind] = sprintf("%5.3g",pvs)
  
 # lev_1[order(pvs),,drop=F]
  outp = data.frame(lev_1)
  #outp
  outp[order(pvs)[1:len],]
 
}


findSigChrom<-function( DE, thresh = 1e-10, go_thresh = 1e-5,nme="FDR1", nme2="chrs"){
  mn = dim(DE)[1]
  pvs = DE[,grep(nme, names(DE))[1]]
  chr_ind = grep(nme2,names(DE))
  sig =  which(pvs<thresh)
  if(length(sig)==0) return(NULL)
  lev_all =getlev(DE[,chr_ind, drop=F])
  lev_ = getlev(DE[sig,chr_ind, drop=F])
  go_todo = lev_[,1]
  inds_m = match(as.character(lev_[,1]), as.character(lev_all[,1]))
  lev_ = cbind(lev_,lev_all[inds_m,1:2,drop=F])
  chrs = as.character(lev_[,1])
  lev_1 = cbind(chrs,t(apply(lev_,1,.phyper2,  k = length(sig), mn = mn)))
  pv_ind = which(dimnames(lev_1)[[2]]=="pv")
  pvs = as.numeric(lev_1[,pv_ind])
  len = length(which(pvs<go_thresh))
  if(len==0) return (NULL)
  outp = data.frame(lev_1)
  
  outp1 = outp[order(pvs)[1:len],]
 # attr(outp1,"inds") = order(pvs)[1:len]
  outp1
}

#.qqplot<-function(DE1, nme="p_lt"){
#  i = which(dimnames(DE1)[[2]]==nme)[1]
#  expected = -log10(seq(1:dim(DE1)[1])/dim(DE1)[1])
#  observed = -log10(DE1[,i])
#  plot(expected,observed, main = nme)
#}



.qqplot<-function(pvals1,log=F,min.p = 1e-20,main="", add=F, col=2){
  pvals = pvals1[!is.na(pvals1)]
  expected = -log10(seq(1:length(pvals))/length(pvals))
  observed = if(log) sort(pvals)/log(10) else log10(sort(pvals))
  observed[observed<log10(min.p)]  = log10(min.p)
  if(add){
    lines(expected, -observed,main=main,type="p", col=col)
    
  }else{
  plot(expected, -observed,main=main)
  }
}
.log10p<-function(pv, log,min.p) {
  pv1 =  if(log) pv/log(10) else  log10(pv)
  pv1[pv1<log10(min.p)] = log10(min.p)
  pv1
}
.vis<-function(dpth, i,min.p = 1e-20,log=F, chroms=NULL, xlim = NULL){
  pv_inds = grep("pv", names(dpth))
  pvs = dpth[,pv_inds,drop=F]
#  for(i in 1:(dim(pvs)[2])){
    pvs[,i] = .log10p(pvs[,i], log=log, min.p = min.p)
 # }
    ranges = attr(dpth, "ranges")
    pos=dpth$pos
    
    if(is.null(pos)) pos = dpth$start
    offset = 0
    for(j in 1:(dim(ranges)[1])){
      r2 = as.numeric(ranges[j,1:2])
      pos[r2[1]:r2[2]] = pos[r2[1]:r2[2]]+offset
      offset =  pos[r2[2]]+10e6
      
    }
    chr_inds = attr(dpth,'chr_inds')
    #print(chr_inds)
    toplot = cbind(pos, -pvs[,i], chr_inds)  #, col=0,ylim = c(0,-min(pvs[,i])))
    
    if(!is.null(chroms)){
   
    inds_ = which(chr_inds %in% chroms)
    print(length(inds_))
    toplot1=cbind(pos[inds_], pvs[inds_,i],chr_inds[inds_])  #, col=0,ylim = c(0,-min(pvs[,i])))
    plot(toplot1[,1:2],col=0,ylim = c(0,-min(pvs[,i], na.rm=T) ), xlim = xlim)
    
    }else{
      
    plot(toplot[,1:2],col=0,ylim = c(0,-min(pvs[,i],na.rm=T) ), xlim = xlim)
    
    }
    for(j in 1:(dim(ranges)[1])) {
      r2 = as.numeric(ranges[j,1:2])
      
      if(is.null(chroms) || ranges[j,3] %in% chroms){   
        print(j)
        print(r2)
        lines(toplot[r2[1]:r2[2],1:2],type="p", col=j)  #col=0,ylim = c(0,-min(pvs[,i]) ),col=0)
      }
    }
    dimnames(toplot)[[2]] = c("pos","pv", "chr_ind")
    invisible(data.frame(toplot))
#    invisible(minp)
}

.phyper2<-function(vec1,k,mn){
  vec = as.numeric(vec1[c(2,length(vec1))])
  #print(vec)
  m = vec[2]
  n = mn-vec[2]
  #print(c(vec[1], m,n,k))
  pv =  phyper(vec[1], m = m, n = n, k = k, lower.tail=F) + dhyper(vec[1], m = m, n = n, k=k )
  vals = qhyper(0.99, m = m, n = n, k = k)
  enrich = (vec[1]/k)/(vec[2]/mn)
  enrich1 = vec[1]/vals
  
  res = c(vec, pv, enrich, enrich1)
  names(res) = c("c1", "c2", "pv", "enrich", "enrich99") 
  res
}
