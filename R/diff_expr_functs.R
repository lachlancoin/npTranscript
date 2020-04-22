
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
      pbb_lower =1- pbinom(x[!zeros],size = sizex,prob = proby[!zeros],log.p=F, lower.tail=lower.tail)
      pbb_lower[pbb_lower<=0] = 0
      point =  dbinom(x[!zeros],size = sizex,prob = proby[!zeros],log.p=F)
      pbb[!zeros] =   pbb_lower + point; 
      if(log) pbb[!zeros] = log(pbb[!zeros,1] )
    }else{
      pbb[!zeros] = pbinom(x[!zeros],size = sizex,prob = proby[!zeros],log.p=log, lower.tail=lower.tail)
      
    }
  }else{
    if(lower.tail==FALSE){
      pbb_lower =1- pbetabinom.ab(x[!zeros],size = sizex,shape1 = shape1[!zeros],shape2 =shape2[!zeros],log=F)
      pbb_lower[pbb_lower<=0] = 0
      point = dbetabinom.ab(x[!zeros],size = sizex,shape1 = shape1[!zeros],shape2 =shape2[!zeros],log=F)
      pbb[!zeros] = pbb_lower + point;
      if(log) pbb[!zeros] = log(pbb[!zeros] )
    }else{
      pbb[!zeros] = pbetabinom.ab(x[!zeros],size = sizex,shape1 = shape1[!zeros],shape2 =shape2[!zeros],log=log)
    }
  }
  if(log){
    pbb[!zeros] =  pbb[!zeros]/log(10)
  }
else{
  pbb[pbb>1] = 1
} 
  pbb
}

chisqCombine<-function(pv){
  #if(TRUE) return(pv[1])
  nonNA = which(!is.na(pv))
  
  if(length(nonNA)==0) return(NaN)
  pchisq(sum(unlist(lapply(pv[nonNA],qchisq,lower.tail=F, df=1 ))),df=length(nonNA), lower.tail=F)
}

#which x is significiantly more or less than expected given y
#if(lower.tail=T returns p(x<=y) else p(x>=y)
##ASSUMES MATCHED DATA BETWEEN CONTROL  AND INFECTED
DEgenes<-function(df,inds_control, inds_infected, edgeR = F, log=F,binom=F, lower.tail = T,reorder=T){
  if(!edgeR){
    pvalsM = matrix(NA,nrow = dim(df)[1], ncol = length(inds_control))
    for(i in 1:length(inds_control)){
        x = df[,inds_control[i]]
        y = df[,inds_infected[i]]
        pvalsM[,i] = betaBinomialP(x,y, binom=binom, lower.tail=lower.tail,log=log)
        
    }
    
    pvals = apply(pvalsM, 1, chisqCombine)
  
  }else{
    pvals = DE_egdeR(df, inds_control, inds_infected)
  }
 FDR = p.adjust(pvals, method="BH");

  probX = (x/sum(x))*1e6
  probY = (y/sum(y))*1e6
  probX1 = (x+0.5)/sum(x+.5)
  probY1 = (y+0.5)/sum(y+.5)
  ratio1 = probX1/probY1
  
  
 
  output =  data.frame(pvals,FDR,probX, probY, ratio1,x,y, df)
 # print(names(output))

  names(output)[names(output) %in% c("x","y") ] = names(df)[inds]
  names(output)[names(output) %in% c("probX","probY") ] = paste( names(df)[inds]," TPM", sep="")
  if(reorder){
    orders =order(pvals)
    
    output = output[orders,]
  }
  output
#  output[orders[,1],,drop=F]
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
  qlf$table$P
}

getDescr<-function(DE,mart, thresh = 1e-10, prefix="ENSCS"){
  inds = which(DE$FDR<thresh)
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


getlev<-function(x, todo = NULL){
  lev = levels(as.factor(as.character(x)))
  cnts = rep(0, length(lev))
  for(i in 1:length(lev)){
    cnts[i] = length(which(x==lev[i]))
  }
  res = data.frame(lev,cnts)[order(cnts, decreasing = T),, drop=F]
  if(is.null(todo)) return(res)
  
  matr = data.frame(lev=todo, cnts=rep(0,length(todo)))
  # print(dim(matr))
  # print(dim(res))
  if(dim(res)[1]>1){
    matr[match(res[,1], matr[,1]),2] = res[res[,1] %in% matr[,1],2]
    dimnames(matr)[[2]] = dimnames(res)[[2]]
  }else{
    #print(todo)
    matr[match(res[,1], matr[,1]),2] =res
  }
  
  matr
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

findGenesByChrom<-function(DE,chrom="MT", fdr_thresh = 1e-10){
  inds = which(DE$chrom== chrom & DE$FDR<fdr_thresh)
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

.readTranscriptsHost<-function(infilesT, 
                         target= list(count0="numeric", count1 = "numeric",chrom="character", leftGene="character", rightGene="character")
              ,prefix="ENSC"                                   
  ){
  header = names(read.table( infilesT,sep="\t", head=T, nrows = 3, comment.char='#'))
  inf = scan(infilesT, nlines=1, what=character())
  inf = sub('#','',inf)
  types = unlist(lapply(inf, function(x) rev(strsplit(x,"_")[[1]])[1]))
  header_inds = match(names(target),header)
  colClasses = rep(NULL, length(header));
  colClasses[header_inds] = target
  
  transcripts = read.table( infilesT,sep="\t", head=T, comment.char='#', colClasses= colClasses)
  header_inds1 = match(names(target),names(transcripts))
  
  transcripts = transcripts [,header_inds1] 
  names(transcripts)[1:2] = types
  names(transcripts)  = sub("leftGene", "geneID" ,names(transcripts))
  names(transcripts)  = sub("chrom", "chrs" ,names(transcripts))
  geneID = transcripts$geneID
  type = rep(NA, length(geneID))
  type[ grep(prefix, transcripts$geneID)] = "left"
  missing = grep(prefix, transcripts$geneID,inv=T)
  if(length(missing)>0){
    have = grep(prefix,transcripts$rightGene[missing])
    transcripts$geneID[missing[have]] = transcripts$rightGene[missing[have]]
    type[missing[have]] = "right"
  }
  
  #o = order(countAll, decreasing=T)
  #transcripts = transcripts[o,]
  #err_ratio_inds = grep("error_ratio", names(transcripts))
  #transcripts[,err_ratio_inds] =apply(transcripts[,err_ratio_inds,drop=F], c(1,2), function(x) if(is.na(x)) -0.01 else x)
  
  if(length(grep("#", inf))>0) attr(transcripts,"info") = sub("#", "",inf)
  #print(inf)
  type = as.factor(type)
res = cbind(transcripts, type)
res
}


findGenes<-function(goObj,DE,goid, fdr_thresh = 1e-10){
  
  inds =  which(goObj$goids$go_id==goid) 
  genes = goObj$goids[inds,,drop=F]
  inds1 = which(DE$geneID %in% genes$ensembl_gene_id)
  inds2 = which(DE[inds1,]$FDR<fdr_thresh)
 # ge = DE$geneID[inds1[inds2]]
#  print(genes[which(genes$ensembl_gene_id %in% ge),])
  DE[inds1[inds2],]
}

findSigGo_<-function(goObj, DE, fdr_thresh = 1e-10, go_thresh = 1e-5, prefix="ENSC"){
  goids = goObj$goids 
  ensg =grep(prefix, DE$geneID, v=T)
  goidx = rep(FALSE,length(goids$ensembl_gene_id ))
  
  
  pvs = DE$FDR
  sig =  which(pvs<fdr_thresh)
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


findSigChrom<-function( DE, fdr_thresh = 1e-10, go_thresh = 1e-5){
  ensg =DE$geneID
  pvs = DE$FDR
  sig =  which(pvs<fdr_thresh)
  lev_all =getlev(DE$chrom)
  lev_ = getlev(DE$chrom[sig])
  go_todo = lev1[,1]
  inds_m = match(as.character(lev_[,1]), as.character(lev_all[,1]))
  lev_ = cbind(lev_,lev_all[inds_m,1:2,drop=F])
  chrs = as.character(lev_[,1])
  lev_1 = cbind(chrs,t(apply(lev_,1,.phyper2,  k = length(sig), mn = length(ensg))))
  pv_ind = which(dimnames(lev_1)[[2]]=="pv")
  pvs = as.numeric(lev_1[,pv_ind])
  len = length(which(pvs<go_thresh))
  outp = data.frame(lev_1)
  outp[order(pvs)[1:len],]
  
}

#.qqplot<-function(DE1, nme="p_lt"){
#  i = which(dimnames(DE1)[[2]]==nme)[1]
#  expected = -log10(seq(1:dim(DE1)[1])/dim(DE1)[1])
#  observed = -log10(DE1[,i])
#  plot(expected,observed, main = nme)
#}


.qqplot<-function(pvals1){
  pvals = pvals1[!is.na(pvals1)]
  expected = -log10(seq(1:length(pvals))/length(pvals))
  observed = -log10(sort(pvals))
  plot(expected, observed)
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

