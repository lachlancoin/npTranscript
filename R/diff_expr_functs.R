.metapv<-function(pvi){
  nonNA = which(!is.na(pvi))
  if(length(nonNA)==0) return(NA)
  pchisq(sum(unlist(lapply(pvi[nonNA],qchisq,lower.tail=F, df=1 ))),df=length(nonNA), lower.tail=F)
  
}

.xlim<-function(x, pthresh = 1e-5, col="p.adj"){
  i = which(names(x)==col)[1]
  x1 = x[order(x[,i]),]
  x2 = x1[x1[,i]<=pthresh,]
  x2
}

.joinSS<-function(l, sort=F){
  res = NULL
  for(i in 1:length(l)){
    type=rep(names(l)[[i]],dim(l[[i]])[[1]])
    toadd = l[[i]]
    toadd[[length(toadd)+1]] = type
    names(toadd)[length(toadd)] = "type"
    if(i==1){
      res = toadd
    }else{
      res = rbind(res, toadd)
    }
  }
  if(sort){
    res$ORFs = as.numeric(as.character(res$ORFs))
    res = res[order(res$ORFs),]
  }
res
}

.mergeDepthByPos<-function(depth){
  ddim = dim(depth)
  dimnames = dimnames(depth)
  dn =dimnames(depth)[[2]]
  duplicate = duplicated(dn)
  if(length(which(duplicate))==0) return(depth)
  duplicate_lev = levels(as.factor(dn[duplicate]))
  duplicate_ = dn %in% duplicate_lev
  duplicate_ind = which(duplicate_)
  nonduplicate_ind = which(!duplicate_)
  depth_dupl = depth[,duplicate_ind,]
  dn_dupl = dimnames(depth_dupl)[[2]]
  ddim[2] = length(duplicate_lev)
  dimnames[[2]] = duplicate_lev
  depth2 = array(NA, dim = ddim, dimnames = dimnames)
  dn_dupl = dimnames(depth_dupl)[[2]]
  for(i in 1:length(duplicate_lev)){
    depth2[,i,] = apply(depth_dupl[,which(dn_dupl==duplicate_lev[i]),,drop=F],c(1,3),sum)
  }
  abind(depth[, nonduplicate_ind,,drop=F],depth2,along=2)
}

readH5_h<-function(h5file, df, filenames, thresh =100,tokeepi = NULL, log=F){
  IDS = df$ID
  chrs = df$chrs
  header = h5read(h5file,"header")
  
  if(header[2]!="base") header = c(header[1], "base", header[-1])
  depth_inds = grep("depth[0-9]", header)
  error_inds = grep("errors[0-9]", header)
  dinds  = grep("depth", header)
  einds  = grep("error", header)
  header[dinds] = paste("depth",filenames,sep="_")
  header[einds] =  paste("error",filenames,sep="_")
  countT = apply(df[,grep("count[0-9]", names(df)),drop=F],1,min)
  pos_ind = which(header=="pos")
  names = h5ls(h5file)$name
  inds = which(countT>thresh & df$ID %in% names)
  array_nmes = list(c("depth","error"), c(), c(sub("depth_","",header[dinds]) ))
  mat_all  = array(dim = c(2,0, length(dinds)),dimnames=array_nmes)
  #  mat_all_e = data.frame(matrix(nrow = 0, ncol = length(einds)+1))
  
  if(length(inds)>0){
  for(i in 1:length(inds)){
    print(paste(i, length(inds), sep=" of "))
    ID = as.character(IDS[inds[i]])
    chr =  as.character(chrs[inds[i]])
    if(is.null(chr) || length(chr)==0) chr = as.character(df$ORFs[inds[i]])
    
    mat = t(h5read(h5file,as.character(ID)))
    countT1 = apply(mat[,depth_inds,drop=F],1,sum)
    indsi = which(countT1>thresh)
    if(length(indsi)>0){
    mat1 = apply(mat[indsi,,drop=F],c(1,2),as.numeric)
   # mat1[mat1[,2]==0,2] = 'A'
  #  mat1[mat1[,2]==1,2] = 'C'
  #  mat1[mat1[,2]==2,2] = 'G'
  #  mat1[mat1[,2]==3,2] = 'T'
   
    IDv = apply(cbind(rep(chr, length(indsi)),mat1[,1:2,drop=F]),1,paste,collapse=".")
    array_nmes[[2]] = IDv
    mat2 = array(dim = c(2,length(indsi), length(dinds)), dimnames=array_nmes)
    mat2[1,,] = mat1[,dinds,drop=F]
    mat2[2,,] =mat1[,einds,drop=F]
   
    
    if(dim(mat2)[2]>0){
        mat_all= abind(mat_all,mat2, along=2)
      #  mat_all_e = rbind(mat_all_e,mat2_e)
      
    }
    }
  } 
  }
mat_all
  #print(head(mat_all[mat_all$pv2<1e-100,]))
  
}

# cumulative 
readH5_c<-function(h5file, df, filenames, thresh=0, log=F,tokeepi = NULL, seqlen=30000){
  IDS = df$ID
  chrs = df$chrs
  header = h5read(h5file,"header")
  if(header[2]!="base") header = c(header[1], "base", header[-1])
  depth_inds = grep("depth[0-9]", header)
  error_inds = grep("errors[0-9]", header)
  dinds  = grep("depth", header)
  einds  = grep("error", header)
  header[dinds] = paste("depth",filenames,sep="_")
  header[einds] =  paste("error",filenames,sep="_")
  #countT = apply(df[,grep("count[0-9]", names(df)),drop=F],1,min)
  pos_ind = which(header=="pos")
  names = h5ls(h5file)$name
  inds = which( df$ID %in% names)
  array_nmes = list(c("depth","error"), 1:seqlen, c(sub("depth_","",header[dinds]) ))
  mat_all  = array(0,dim = c(2,seqlen, length(dinds)),dimnames=array_nmes)
  
  if(length(inds)>0){
    for(i in 1:length(inds)){
      print(paste(i, length(inds), sep=" of "))
      ID = as.character(IDS[inds[i]])
      chr =  as.character(chrs[inds[i]])
      if(is.null(chr) || length(chr)==0) chr = as.character(df$ORFs[inds[i]])
      mat = t(h5read(h5file,as.character(ID)))
      mat_all[1,mat[,1],] =  mat_all[1,mat[,1],]+mat[,dinds]
      mat_all[2,mat[,1],] =  mat_all[2,mat[,1],]+mat[,einds]
    } 
  }
  mat_all
  #print(head(mat_all[mat_all$pv2<1e-100,]))
  
}

# split into transcripts
readH5_s<-function(h5file, df, filenames, thresh=0, log=F, src_i=1 ,seqlen=30000){
  IDS = df$ID
  chrs = df$chrs
  header = h5read(h5file,"header")
  if(header[2]!="base") header = c(header[1], "base", header[-1])
  depth_inds = grep("depth[0-9]", header)
  error_inds = grep("errors[0-9]", header)
#  tokeepi = grep(grep("depth_",grep(tokeep, header,v=T),v=T),header)
  dinds  = grep("depth", header)
  einds  = grep("error", header)
 # header[dinds] = paste("depth",filenames,sep="_")
#  header[einds] =  paste("error",filenames,sep="_")
  #countT = apply(df[,grep("count[0-9]", names(df)),drop=F],1,min)
  pos_ind = which(header=="pos")
  names = h5ls(h5file)$name
  inds = which( df$ID %in% names)
  array_nmes = list(c("depth","error"), 1:seqlen, df$ORFs )
  num_trans = dim(df)[[1]]
  mat_all  = array(0,dim = c(2,seqlen, num_trans),dimnames=array_nmes)
  
  if(length(inds)>0){
    for(i in 1:length(inds)){
      print(paste(i, length(inds), sep=" of "))
      ID = as.character(IDS[inds[i]])
      chr =  as.character(chrs[inds[i]])
      if(is.null(chr) || length(chr)==0) chr = as.character(df$ORFs[inds[i]])
      mat = t(h5read(h5file,as.character(ID)))
      mat_all[1,mat[,1],i] =  mat[,dinds[src_i]]
      mat_all[2,mat[,1],i] =  mat[,einds[src_i]]
    } 
  }
  mat_all
  #print(head(mat_all[mat_all$pv2<1e-100,]))
  
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
  ratio1 = rep(NA, length(control))
  ratio2 = rep(NA, length(control))
  
  for(j in 1:length(control)){
    indsj = c(indsDepth[control[j]], indsError[control[j]], indsDepth[case[j]], indsError[case[j]])
    #indsj2 = c(indsDepth[case[j]], indsError[case[j]], indsDepth[control[j]], indsError[control[j]])
    pv1[j] = betaBinomialP1(v,indsj , binom=binom, lower.tail=T, log=log)
    pv2[j] = betaBinomialP1(v, indsj, binom=binom, lower.tail=F, log=log)
    r =.ratios(v,indsj )
    ratio1[j] = r[2]*100
    ratio2[j] = r[1]*100
  }
 # pv1m =chisqCombine(pv1,log=log)
#  pv2m = chisqCombine(pv2,log=log)
  c(pv1,pv2,ratio1, ratio2)
 # c(pv1m,pv2m)
#  2*min(pv1m, pv2m, na.rm=T)
}

.ratios<-function(v,ord){
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
  
  c(x1/sizex,y/size)
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
   # print(paste(x1,sizex, shape1, shape2))
    pv_res = pbetabinom.ab(x1,size = sizex,shape1 = shape1[[1]],shape2 =shape2[[1]],log=log)
    
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

chisqCombine<-function(pv,log=F){
  #if(TRUE) return(pv[1])
  nonNA = which(!is.na(pv))
  
  if(length(nonNA)==0) return(NaN) else if(length(nonNA)==1) return (pv[nonNA]);
  if(log){
    resp = pchisq(sum(unlist(lapply(exp(pv[nonNA]),qchisq,lower.tail=F, df=1 ))),df=length(nonNA), lower.tail=F,log=T)
    
  }else{
  resp = pchisq(sum(unlist(lapply(pv[nonNA],qchisq,lower.tail=F, df=1 ))),df=length(nonNA), lower.tail=F)
  }
}

.write<-function(DE1, resdir, filename="results.csv", numeric_inds=1:9){
  if(length(numeric_inds)>0){
  DE1[,numeric_inds] = apply(DE1[,numeric_inds], c(1,2), function(x) sub(' ' , '', sprintf("%5.3g",x)))
  }
  write.table(DE1[attr(DE1,"order"),],file=paste(resdir,filename,sep="/") , quote=F, row.names=F, sep="\t", col.names=T)
}
.exact<-function(v, nrow=2, ncol=2){
  if(sum(v)==0) return (c(NA,NA,NA,NA))
 # print(v)
  f = fisher.test(matrix(v,nrow=nrow,ncol=ncol))
  c(f$p.value,f$estimate, v[2,1]/v[1,1], v[2,2]/v[1,2])
}
.chisq<-function(v, nrow=2,ncol=2){
  if(sum(v)==0) return (c(NA,NA,NA,NA))
  f = chisq.test(matrix(v,nrow=nrow, ncol=ncol))
  c(f$p.value, NA,v[2,1]/v[1,1], v[2,2]/v[1,2])
}
.extractFromDepth<-function(dpth, indsc=1:2,ORF=NULL,pos=NULL){
  dnmes = dimnames(dpth)[[2]]
  if(!is.null(ORF)){
    dpth =dpth[,grep(ORF,dnmes),,drop=F]
  }
  if(!is.null(pos)){
    dpth =dpth[,grep(pos,dimnames(dpth)[[2]]),,drop=F]
  }
  if(!is.null(indsc)){
   dpth = dpth[,,indsc,drop=F]
  }
  res = dpth
  ratio = apply(res,c(2,3),function(x) x[2]/x[1])
  res1 = abind(res,ratio, along=1)
  dimnames(res1)[[1]][3] = "ratio"
  res1
}
.processDM<-function(depth, filenames, control_names, infected_names, thresh_min =100, method=.chisq, plot=F, adjust=getOption("np.adjustMethod","BH")){
  if(is.null(depth)) return(NULL)
  inf = dimnames(depth)[[3]]
  DE = list()
  method1 = if(method=="chisq.test") .chisq else if(method=="fisher.test") .exact else method
  inds = which(filenames %in% c(control_names, infected_names))

  row_inds = apply(depth[1,,inds],1,min)>thresh_min
 # df1 = df[,row_inds,inds,drop=F]
  if(length(which(row_inds))<=0) return (DE)
  for(i in 1:length(control_names)){
    DE[[i]] = DEdepth(depth[,row_inds,,drop=F], control_names[i], infected_names[i], tojoin=1:3, method=method1, adjust=adjust)
  }
  nmes=apply(cbind(infected_names, control_names),1,function(x).largestPrefix(x[1],x[2]))
  
  if(length(control_names)>1){
   DE[[length(DE)+1]] = DEdepth(depth[,row_inds,], control_names, infected_names, tojoin=1:3,method=method1, adjust=adjust)
    nmes = c(nmes,"meta")
  }
  names(DE) =  nmes
for(i in 1:length(DE))  attr(DE[[i]],"nme") = names(DE)[i]
 # if(length(DE)==1) return(DE[[1]])
    invisible(rev(DE))
}

.printOptions<-function(out, prefix="np."){
  nme=grep(prefix,names(options()),v=T)
  vals = unlist(lapply(nme,function(x) getOption(x)))
  names(vals) = nme
  write.table(as.matrix(vals),out,row.names=T,col.names=F)
}
DEdepth<-function(df1,control_names, infected_names,tojoin=1:3, maxLogFC=10, method=chisq.test, adjust  = "none"){
 inds = which(dimnames(df1)[[3]] %in%  c(control_names, infected_names))

  nme = dimnames(df1)[[3]]
  pvs =data.frame(matrix(NA, ncol = 4*length(control_names), nrow = dim(df1)[[2]]))
  st_col = 1
  inds_all = c()
 # print(length(row_inds))
  #if(length(which(row_inds))>0){
  for(i in 1:length(control_names)){
    nmes_i = c(grep(control_names[i], nme,v=T),grep(infected_names[i], nme,v=T))
    inds_i = which(nme %in% nmes_i)
    inds_all = c(inds_all, inds_i)
    dfi =  df1[,,inds_i,drop=F]
   
    indsk = st_col:(st_col+4-1)
    st_col = st_col+4
  #print(dim(dfi))
    pv1 = t(apply(dfi,2,method))
    pvs[,indsk] =pv1
      #t(apply(dfi, 1, betaBinomialP2, depth_inds_i, error_inds_i,1,2,binom=F, log=F))
    names(pvs)[indsk] = paste(control_names[i], c("pv","OR","ratio_control","ratio_case"),sep="_")
  }
  #}
  #print(tojoin)
  control_ratio=apply(pvs[,grep("ratio_control", names(pvs)),drop=F],1, mean)
  infected_ratio=apply(pvs[,grep("ratio_case", names(pvs)),drop=F],1, mean)
  ORFs = dimnames(df1)[[2]]
  logFC = log2(infected_ratio/control_ratio)
  logFC[is.na(logFC)] = 0
  logFC[!is.na(logFC) & logFC< -maxLogFC]=-maxLogFC
  logFC[!is.na(logFC) & logFC>maxLogFC]=maxLogFC
  pvals=apply(pvs[,grep("pv", names(pvs)),drop=F],1, min)
  OR=apply(pvs[,grep("OR", names(pvs)),drop=F],1, min)
  
  p.adj = p.adjust(pvals,method=adjust)
 pvs=data.frame( cbind(logFC,OR, p.adj,pvals))
 # if(length(tojoin)>0)  pvs= cbind(df[row_inds,tojoin,drop=F],pvs)
# df2 = df1[1,, inds_all,drop=F]
# df3 = df1[2,, inds_all,drop=F]
 result = cbind(ORFs,pvs) #, df2,df3)

 return (result)
}

##edgeR
### qlf = DE_egdeR(df, control_inds, infected_inds)
#pvals1 = qlf$table$P
# pvals1 = apply(pvalsM1, 1, chisqCombine,log=log)
#pvals2 = apply(pvalsM2, 1, chisqCombine,log=log)

#which x is significiantly more or less than expected given y
#if(lower.tail=T returns p(x<=y) else p(x>=y)
##ASSUMES MATCHED DATA BETWEEN CONTROL  AND INFECTED
DEgenes<-function(df,control_names,infected_names,  type="lt", binom=F, log=F
                  ){
  lower.tail = T
  ORFs = as.character(df$ORFs)
  if(!is.null(df$Name)){
  dfn_ind = !is.na(df$Name)
  ORFs[dfn_ind]=as.character(df$Name[dfn_ind])
  }
  grp=df$grp
  if(is.null(grp)) grp = df$type_nme

  control_inds = rep(NA, length(control_names))
  infected_inds = rep(NA, length(infected_names))
  results = list();
  for(i  in 1:length(control_inds)){
    print(i)
	    control_inds[i] = which(names(df)==control_names[i])
      infected_inds[i] = which(names(df)==infected_names[i])
      x = df[,control_inds[i]]
      y = df[,infected_inds[i]]
      pvals1 = betaBinomialP(x,y, binom=binom, lower.tail=lower.tail,log=log)
      pvals2 = betaBinomialP(y,x, binom=binom, lower.tail=lower.tail,log=log)
     pvals = apply(cbind(pvals1,pvals2),1,min, na.rm=T)
    p.adj = p.adjust(pvals,method="BH")
    tpm_control = (x/sum(x))*1e6
    tpm_infected = (y/sum(y))*1e6
    probX1 = (x+0.5)/sum(x+.5)
    probY1 = (y+0.5)/sum(y+.5)
    ratio1 = probY1/probX1
    logFC = log(ratio1)/log(2)
    results[[i]] =  data.frame(ORFs=ORFs,grp = grp, pvals, p.adj, pvals1,pvals2, tpm_control, tpm_infected, ratio1,logFC, sum_control=x,sum_infected=y)
  }
  names(results)=apply(cbind(infected_names, control_names),1,paste, collapse=" v ")
                       #function(x).largestPrefix(x[1],x[2]))
 results
}

.meta<-function(DE1){
  DE2 = data.frame(unlist(DE1,recursive=F))
  pvals1 = apply(DE2[,grep(".pvals1",names(DE2))],1,chisqCombine,log=F)
  pvals2 = apply(DE2[,grep(".pvals2",names(DE2))],1,chisqCombine,log=F)
  pvals = apply(cbind(pvals1,pvals2),1,min)
  p.adj = p.adjust(pvals,method="BH")
 x = apply(DE2[,grep("sum_control",names(DE2))],1,sum,log=F)
 y = apply(DE2[,grep("sum_infected",names(DE2))],1,sum,log=F)
 tpm_control = (x/sum(x))*1e6
 tpm_infected = (y/sum(y))*1e6
 probX1 = (x+0.5)/sum(x+.5)
 probY1 = (y+0.5)/sum(y+.5)
 ratio1 = probY1/probX1
 logFC = log(ratio1)/log(2)
 grp = DE1[[1]]$grp
 ORFs=DE1[[1]]$ORFs
 data.frame(ORFs, grp, pvals, p.adj, pvals1,pvals2, tpm_control, tpm_infected, ratio1,logFC, sum_control=x,sum_infected=y)

}

.largestPrefix<-function(str1,str2){
  for(i in 1:nchar(str1)){
    if(substr(str1,1,i)!=substr(str2,1,i)) break;
  }
  substr(str1,1,i-1)
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
  orfs = gsub("[-+]","",transcript$ORFs)
  geneID= as.character(unlist(lapply(strsplit(orfs,";"), function(v) v[1])))
 # rightGene = as.character(unlist(lapply(strsplit(orfs,";"), function(v) v[length(v)])))
  #indsL =  grep(prefix,geneID,inv=T)
  #indsR =  grep(prefix,rightGene)
  #comb = which(indsL %in% indsR)
#if(length(comb)>0) geneID[comb] = rightGene[comb]
  cbind(transcript, geneID)
}

.readCaseControl<-function(f){
  if(!file.exists(f)) return(NULL)
  gfft = read.table(f, sep="\t", header=F, fill=T)
  list(control = gfft[,1], case =gfft[,2])
}

.addAnnotation<-function(annotfile, transcripts,grp, colid="geneID",nmes = c("chr","ID" , "Name" , "Description","biotype")){
  match_ind = which(names(transcripts)==colid)[1]
  gfft = read.table(annotfile, sep="\t", header=F, fill=T, quote='\"')
  chrind = which(nmes=="chr")
  names(gfft) = nmes
  ID_ind = which(nmes=="ID")
  #gfft[,1] = gsub("transcript:", "", as.character(gfft[,1]))
  mi = match(as.character(transcripts[,match_ind]), as.character(gfft$ID))
  gfft = gfft[mi,]
  gfft$ID = transcripts$ID
  transcripts = cbind(transcripts,grp,gfft[,-c(chrind,ID_ind)])
  return(transcripts)
}

.mergeRows<-function(transcripts1,sum_names = c(), append_names= c(), colid='geneID', collapse="_"){
  if(dim(transcripts1)[[1]]==0) return(NULL)
  sum_names = unique(c(sum_names,"countTotal"));
  ind = which(names(transcripts1) %in% colid)
  ind_s = which(names(transcripts1) %in% sum_names)
  ind_app = which(names(transcripts1) %in% append_names)
  lev = getlev(transcripts1[,ind])
  
  todo = as.character(lev[lev$cnts>1,]$lev)
  if(length(todo)==0) return(transcripts1)
  extract_inds = lapply(todo, function(x) which(transcripts1[,ind]==x))
  torem = sort(unlist(extract_inds))
  subt = transcripts1[unlist(lapply(extract_inds,function(x) x[1])),]
 # print(extract_inds)
  for(i in 1:length(extract_inds)){
    subind = extract_inds[[i]]
    matr_join = transcripts1[subind,ind_app,drop=F]
    subt[i,ind_s] =apply(apply(transcripts1[subind,ind_s,drop=F],c(1,2),as.numeric),2,sum)
    joined = apply(matr_join,2,paste, collapse=collapse)
  #  print(matr_join)
    subt[i,ind_app] =joined
  }
  resu = rbind(transcripts1[-torem,,drop=F], subt)
 # subt
  resu
}

###obsolete function
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
.combineTranscripts<-function(dfs, attributes
                              ){
                chroms=attributes[[which(names(attributes)=="chroms")]]
                chrom_names=attributes[[which(names(attributes)=="chrom_names")]]
                
  lengs = unlist(lapply(dfs,function(x) if(is.null(x)) 0 else dim(x)[1]))
  ncol = dim(dfs[[1]])[2]
  #lengs = lengs[inds][ord]
  res = data.frame(matrix(nrow  =sum(lengs), ncol = ncol))
  names(res) = names(dfs[[1]])
  start=1;
  ranges = matrix(nrow = length(dfs), ncol=3)
  for(i in 1:length(dfs)){
  #  print(i)
    lengi =lengs[i]
    for(k in 1:ncol){
      if(is.factor(dfs[[i]][,k])) dfs[[i]][,k] = as.character(dfs[[i]][,k])
    }
    ranges[i,] = c(start, start+lengi-1,chrom_names[i]);
    if(lengi>0) res[ranges[i,1]:ranges[i,2],] = dfs[[i]]
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
 # dimnames(res)[[1]] = res$ID
  
  res
}


.process<-function(transcriptsl, control_names, infected_names, 
                   annotFile = "annotation.csv.gz"){
  transcriptsl = lapply(transcriptsl, .mergeRows,sum_names= c(control_names, infected_names), colid="geneID" )
  transcripts=.combineTranscripts(transcriptsl, attributes)
   dimnames(transcripts)[[1]] = transcripts$ID
  # 
  info = attr(transcripts,'info')
  transcripts=.addAnnotation(annotFile, transcripts, colid="geneID", nmes = c("ID" , "Name" , "Description","biotype"))
  transcripts
}

.processDE1<-function(transcripts, count_names, i1, i2, resdir, top=5, pthresh = 1e-3,plot=F){
  info = attr(transcripts,"info")
  control_names = count_names[i1]
  infected_names = count_names[i2]
  outp = paste("results", info[i1], info[i2], "csv",sep=".")
  type_names = c(control_names[1], infected_names[1])
 .processDE(transcripts,attributes(transcripts), resdir, control_names, infected_names, type_names = type_names, outp= outp, 
                         type=gsub("count_","",paste(type_names,collapse=" vs ")))
}

.volcano<-function(df, logFCthresh = 1.0, top=10, prefix="", useadj=TRUE,exclude=NULL){
  if(dim(df)[[1]]==0) return(NULL)
  if(!is.null(exclude)){
    df = df[!(df$grp %in% exclude),,drop=F]
  }
  if(useadj){
    pthresh = sort(df$p.adj)[top]
ggp<-ggplot(df, aes(x =logFC, y = -log10(p.adj),color = ifelse(abs(logFC)>0.6,"red","grey"))) 
ggp<-ggp+  geom_point() +  xlab(expression("Fold Change, Log"[2]*"")) +  ylab(expression("Adjusted P value, Log"[10]*"")) 
  }else{
    pthresh = sort(df$pvals)[top]
    ggp<-ggplot(df, aes(x =logFC, y = -log10(pvals),color = ifelse(abs(logFC)>0.6,"red","grey"))) 
    ggp<-ggp+  geom_point() +  xlab(expression("Fold Change, Log"[2]*"")) +  ylab(expression(" P value, Log"[10]*"")) 
  }
  print(paste("pthresh",pthresh))
ggp<-ggp+  geom_vline(
    xintercept = c(-0.6,0.6),
    col = "red",
    linetype = "dotted",
    size = 1) 
ggp<-ggp+  geom_hline(
    yintercept = c(-log10(0.01),-log10(0.05)),
    col = "red",
    linetype = "dotted",
    size = 1)
ggp<-ggp+  theme(
  plot.title = element_text(size=8)
)
ggp<-ggp+  theme(legend.position = "none")+
  scale_colour_manual(values = c("grey", "red")) 
if(useadj){
ggp<-ggp+  geom_text_repel(data=subset(df,abs(logFC) >= logFCthresh & p.adj < pthresh),
                  aes(logFC, -log10(p.adj), label = ORFs),size = 3, color="steelblue")
}else{
  ggp<-ggp+  geom_text_repel(data=subset(df,abs(logFC) >= logFCthresh & pvals < pthresh),
                             aes(logFC, -log10(pvals), label = ORFs),size = 3, color="steelblue")
}
ggp<-ggp+ggtitle(sub(" v "," vs ", paste(attr(df,"nme"),prefix)))
ggp
}




.comparisonPlot<-function(DE2, transcriptsl1, inds = c(2,3), countThresh = 100, excl=c("sequins")){
  nme_lfc = grep("logFC",names(DE2),v=T)
  nme_p = grep("p.adj",names(DE2),v=T)
  orf_name = grep("ORF", names(DE2), v=T)
  nmes_lfc =nme_lfc[inds] 
  nmes_p = nme_p[inds]
  inds1_lfc = which(names(DE2) %in% nmes_lfc)
  inds1_p = which(names(DE2) %in% nmes_p)
  countTotal = transcriptsl1$countTotal
  grpname = transcriptsl1$grp
  df = DE2[countTotal>countThresh & !(grpname%in% excl),]
  subset_p =df[head(order(apply(df[,inds1_p[1:2]],1,function(x) max(abs(x))), decreasing=F),10),]
  subset_lfc =df[head(order(apply(df[,inds1_lfc[1:2]],1,function(x) min(abs(x))), decreasing=T),10),]
  df1 = df[apply(df[,grep("p.adj",names(DE2))],1,min) < 1e-2,,drop=F]
  
  print(subset_p[,c(1,inds1_lfc, inds1_p)])
  #print(subset_lfc[,c(1,inds1_lfc, inds1_p)])
  ggp<-ggplot(df1, aes_string(x=nmes_lfc[1], y=nmes_lfc[2]))+geom_point()+ggtitle(paste(nmes_lfc))
  ggp<-ggp+geom_text_repel(data=subset_p,
                           aes_string(x=nmes_lfc[1], y=nmes_lfc[2], label = orf_name[1]),size = 3, color="steelblue")
  
  ggp1<-ggplot(df, aes_string(x=nmes_p[1], y=nmes_p[2]))+geom_point()+ggtitle(paste(nmes_p))
  ggp1<-ggp1+scale_y_continuous(trans='log10')+scale_x_continuous(trans='log10')
  ggp1<-ggp1+geom_text_repel(data=subset_p,
                             aes_string(x=nmes_p[1], y=nmes_p[2], label = orf_name[1]),size = 3, color="steelblue")
  invisible(list(ggp=ggp,ggp1=ggp1))
}

.getAllPairwiseComparisons<-function(info,start=1){
  todo = list()
  if(length(info)==1) return(todo)
  for(i in (start+1):length(info)){
    for(j in start:(i-1)){
      todo[[length(todo)+1]] = c(j,i)
    }
  }
  names(todo) = unlist(lapply(todo, function(x) paste(info[x],collapse=" v ")))
  todo
}

.processDE<-function(transcripts, attributes, resdir, control_names, infected_names,type_names=c("control","infected"), 
                     outp = "results.csv", type=""){
  print(cbind(control_names, infected_names))
   DE1 = try(DEgenes(transcripts, control_names, infected_names));
   if(length(DE1)>1){
     metaDE = .meta(DE1)
     DE1[[length(DE1)+1]]=metaDE
     names(DE1)[[length(DE1)]] = "meta"
     DE1 = rev(DE1)
   }
   for(i in 1:length(DE1)) attr(DE1[[i]], "nme") = names(DE1)[[i]]
  #.write(DE1 ,resdir,outp)
 #if(length(DE1)==1) return (DE1[[1]])
  #.vis(DE1,i=1,min.p=1e-50)
  # .vis(DE1,i=2,min.p=1e-50)
  invisible(DE1)
#  invisible(list(DE1=DE1, transcripts=transcripts))
}
.testIsoformsAll<-function(df, isoforms_i, filenames,control_names, infected_names, n=5, test_func = chisq.test, 
                           adjust=getOption("np.adjustMethod","BH")
                         ){
  ORFs = as.character(df$ORFs)
  dfn_ind = !is.na(df$Name)
  ORFs[dfn_ind]=as.character(df$Name[dfn_ind])
  grp=df$grp
  inds = attr(isoforms_i,"inds")
  if(length(inds)==0) return (NULL)
#  df1 = df[inds,!duplicated(names(df)),drop=F]
  results = list();
  pvm = matrix(NA, ncol = length(control_names), nrow=length(inds))
  lfc=  matrix(NA, ncol = length(control_names), nrow=length(inds))
  for(i  in 1:length(control_names)){
    cols = c(grep(control_names[i], filenames),grep(infected_names[i], filenames))
   dfi =data.frame(t(data.frame(lapply(isoforms_i, .testIsoforms1,cols,  n=2, test_func=test_func))))
   names(dfi)=c("p","logFC")
   p.adj = p.adjust(dfi$p,method=adjust)
   pvm[,i] = dfi$p
   lfc[,i] = dfi$logFC
   results[[i]] = data.frame(ORFs = ORFs[inds], grp=grp[inds],p.adj= p.adj,pvals=dfi$p, "logFC" = dfi$logFC)
  }
  nmes = apply(cbind(infected_names, control_names),1,function(x).largestPrefix(x[1],x[2]))
  if(length(control_names)>1){
  pcomb =  apply(pvm, 1,chisqCombine);
  logfc_comb = apply(lfc, 1,mean, na.rm=T)
  results[[length(results)+1]]= data.frame(ORFs= ORFs[inds], 
                                           grp = grp[inds],pvals = pcomb, p.adj = p.adjust(pcomb, method=adjust), 
                                           logFC=logfc_comb)
  nmes = c(nmes,"meta")
  }
  names(results)=nmes
for(i in 1:length(results)) attr(results[[i]], "nme") =names(results)
  invisible(rev(results))
}

.testIsoforms1<-function(x, cols, n1=2, thresh = 10, test_func=chisq.test){
  tb = x$counts[,cols,drop=F]
  sums = apply(tb,2,sum)
  if(min(sums)<thresh) return (c(NA,NA))
  if(dim(tb)[[1]]>n1){
    tb[n1,] = apply(tb[-(1:(n1-1)),,drop=F],2,sum)
    tb = tb[1:n1,,drop=F]
  }
  ratio=apply(tb,2,function(x) x[1]/sum(x))
  c(test_func(tb)$p.value, log2(ratio[2]/ratio[1]))
}

readIsoformH5<-function(transcripts_, h5file,  depth =1000){
  .ext<-function(x) paste(x[unique(c(0,which(x>0)))],collapse=",")
  header = h5read(h5file,"header")
  names = h5ls(h5file)$name
  inds = which(transcripts_$countTotal>depth)
  if(length(inds)==0) return (NULL)
 
  IDS = lapply(transcripts_[inds,,drop=F]$ID,function(x) strsplit(x,"_")[[1]]);
  IDS = lapply(IDS, function(x) x[which(x %in% names)])
  inds1 = unlist(lapply(IDS, length))>0
  IDS = IDS[inds1]
  trans = vector("list", length = length(IDS))
  if(length(IDS)==0) return (NULL)
  for(i in 1:length(IDS)){
    ID = IDS[[i]]
    mats = lapply(ID, function(x) t(h5read(h5file, as.character(x) )))
    lens = unlist(lapply(mats, function(x)dim(x)[[1]]))
    cnts = array(dim = c(sum(lens),length(header)-1))
    transi = rep("",dim(cnts)[[1]])
    #rown = rep("" , dim(cnts)[[1]])
    #print(ID)
    start = 1;
    for(j in 1:length(ID)){
      mat = mats[[j]]
      rowinds = start:(start+lens[j]-1)
      start = start+lens[j]
      cnts[rowinds,] = mat[,2:length(header),drop=F]
      transi[rowinds] = apply(mat[,-(1:length(header)),drop=F],1,.ext)
    }
    cntT = apply(cnts,1,sum)
    ord = order(cntT, decreasing=T)
    cnts = cnts[ord,,drop=F]
      trans[[i]] = list(breaks = transi[ord], counts = cnts[ord,,drop=F])
  } 
  names(trans) = IDS
  #trans = trans[unlist(lapply(trans, function(x) !is.null(dim(x))))]
  inds2 = unlist(lapply(trans, function(x) dim(x$counts)[1]))>1
  trans = trans[inds2]
  
  attr(trans, "inds")=inds[inds1][inds2] #which(transcripts_$ID %in% names(trans))
  trans
  
}



.readH5All<-function(transcripts, infile, attributes,filenames,  thresh,tokeepi =NULL, readH5_ = readH5_c){
  #filenames=  attributes[which(names(attributes)=="info")][[1]]
  depth_ls = try(readH5_(infile, transcripts, filenames, thresh =thresh,tokeepi = tokeepi, log=F))
  return(.transferAttributes(depth_ls, attributes)[[1]])
}
 

.shorten<-function(str, len=31, split = NULL){
  str = gsub("_pass","",str)
  str = gsub("_","",str)
  str = gsub("vs","v",str)
  str = gsub("infected","inf",str)
  str = gsub("control","ctrl",str)
  
  if(nchar(str)>len){
    if(!is.null(split)){
      strs = paste(unlist(lapply(strsplit(str, split)[[1]], function(x) substr(x,1,14))), collapse=" v ")
      str = strs
    }
    str = substr(str,1,len)
  }
  str
}
.write_xlsx1<-function(sheets, fn, split=" v ", len = 31){
  names(sheets) = unlist(lapply(names(sheets), .shorten, len, split))
  print(names(sheets))
  write_xlsx(sheets, fn)
}
#extracts and rearranges depth
.mergeDepth<-function(i,depths_combined_spliced){
  res =depths_combined_spliced[[1]][,,i,drop=F]
  if(length(depths_combined_spliced)>1){
  for(k in 2:length(depths_combined_spliced)){
    res = abind(res,depths_combined_spliced[[k]][,,i,drop=F],along=3)
  }
  }
  dimnames(res)[[3]] = names(depths_combined_spliced)
  res
}

.readTranscriptsHost<-function(infilesT1, 
                               filter = NULL,
                  target= list(chrom="character", leftGene="character", rightGene="character", start = "numeric", end="numeric", ID="character")
              ,combined_depth_thresh =100                                 
  ){

  header = names(read.table( infilesT1,sep="\t", head=T, nrows = 3, comment.char='#', fill=T))
  ncol = dim(read.table( infilesT1,sep="\t", head=F, nrows = 2,skip = 2, comment.char='#'))[2]
#print(header)
#print(ncol)
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
  
  colClasses = rep(NA, ncol);
  colClasses[header_inds] = target
  #colClass = cbind(rep("numeric", length(extra)), colClasses)
  transcripts = read.table( infilesT1,skip=2,sep="\t", head=F, comment.char='#', colClasses= colClasses)
names(transcripts) = header[1:(dim(transcripts)[2])]
  if(!is.null(filter)){
    for(k in 1:length(filter)){
      nme_ind = grep(names(filter)[k], names(transcripts))
      if(length(nme_ind)==0) stop(paste("not found ",names(filter)[k]))
      indsk=which(transcripts[,nme_ind]==filter[k])
      transcripts = transcripts[indsk,,drop=F]
    }
  }
  header_inds1 = match(names(target),names(transcripts))
  countT = as.numeric(transcripts$countTotal)
  indsk = !is.na(countT) & countT>=combined_depth_thresh
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
  sig =  which(pvs<=thresh)
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


.qqplot1<-function(df, nme,log=F,min.p = 1e-20,main="", add=F, col=2) .qqplot(df[,which(names(df)==nme)],log,min.p,main,add,col)
  

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
.getRanges<-function(matr){
  levs = getlev(matr$chrs)
  levs = levs[order(as.numeric(sub("chr","",as.character(levs[,1])))),]
  ranges = matrix(NA, ncol=3, nrow = dim(levs)[1])
  for(i in 1:(dim(levs)[1])){
    chr = as.character(levs[i,1])
    inds = which( matr$chrs==chr)
    ranges[i,] = c(min(inds), max(inds), chr)
  }
  ranges
}
.vis<-function(dpth, nme,min.p = 1e-50,log=F, chroms=NULL, xlim = NULL, funct = .log10p, usePos=F){
  
  pv_inds = which(names(dpth) %in% nme)
  pvs = dpth[,pv_inds[1]]
  if(!is.null(funct)){
   pvs = funct(pvs, log=log, min.p = min.p)
  }
    ranges = .getRanges(dpth)
    pos=dpth$pos
    if(is.null(pos)) pos = dpth$start
    offset = 0
    chrs = dpth$chrs
    lev= ranges[,3]
  chr_inds = unlist(lapply(chrs, function(x) which(lev==x)))
    for(j in 1:(dim(ranges)[1])){
       r2 = as.numeric(ranges[j,1:2])
       if(usePos){
        pos[r2[1]:r2[2]] = pos[r2[1]:r2[2]]+offset
       }else{
         pos[r2[1]:r2[2]] = (r2[1]:r2[2]) 
       }
        offset =  pos[r2[2]]
        if(length(chroms)==1) offset =0
        chrs[r2[1]:r2[2]] =ranges[j,3]
      
    }
   
    toplot = cbind(pos, pvs, chr_inds)  #, col=0,ylim = c(0,-min(pvs[,i])))
  
    if(is.null(chroms)) toincl  = 1:length(pos) else toincl = which(chrs %in% chroms)
    if(!is.null(chroms)) xlab  =paste(chroms) else xlab = paste(sub("chr","",ranges[,3]))
    
    plot(toplot[toincl,1:2],col=0, xlim = xlim,main=nme, xlab = xlab)
    for(j in 1:(dim(ranges)[1])) {
      if(is.null(chroms) || ranges[j,3] %in% chroms ){
      r2 = as.numeric(ranges[j,1:2])
        print(j)
        print(r2)
       # print(toplot[r2[1]:r2[2],1])
        lines(toplot[r2[1]:r2[2],1:2],type="p", col=j)  #col=0,ylim = c(0,-min(pvs[,i]) ),col=0)
      }
      
    }
    dimnames(toplot)[[2]] = c("pos","pv","chr_inds")
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
