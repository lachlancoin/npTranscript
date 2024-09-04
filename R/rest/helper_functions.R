.mergeC<-function(v){
  res = c()
  if(length(v)>0){
    for(i in 1:length(v)) res = c(res,v[[i]])
  }
  res
}


.merge1<-function(t,num_cols = c(),uniq_cols=c(), addName=NULL, rowNames=F){
  t=t[unlist(lapply(t, nrow))>0]
  t=t[unlist(lapply(t, length))>0]
  if(length(t)==0) return(NULL)
  
  t1 = t[[1]]
  if(!is.null(addName)){
    nme = rep(names(t)[[1]], nrow(t1))
    t1 = cbind(t1, nme)
  }
  #if(length(t)==1) return(t1)
  if(length(t)>1){
    for(i in 2:length(t)){
      t_i=t[[i]]
      if(!is.null(addName)){
        nme = rep(names(t)[[i]], nrow(t_i))
        t_i = cbind(t_i, nme)
      }
      t1 = rbind(t1,t_i)
    }
  }
  if(rowNames){
    dimnames(t1)[[1]] = .mergeC(lapply(t, function(t_)dimnames(t_)[[1]]))
  }else{
    dimnames(t1)[[1]]=1:nrow(t1)
  }
  t2 = data.frame(t1)
  names(t2) = dimnames(t1)[[2]]
  if(!is.null(addName)){
    names(t2)[ncol(t2)]=addName
  }
  
  if(length(uniq_cols)>0){
    facts = factor(apply(t2[,which(names(t2) %in% uniq_cols),drop=F],1,paste, collapse="."))
    t2 = t2[!duplicated(facts),,drop=F]
  }
  num_cols = num_cols[which(num_cols %in% names(t2))]
  lev_cols =names(t2)[!( names(t2)%in% num_cols)]
  for(j in num_cols){
    t2[[j]] = as.numeric(t2[[j]])
  }
  for(j in lev_cols){
    t2[[j]] = factor(unlist(t2[[j]]))
    
  }
  
  t2
}


whichpart1<-function(angle, n=10, one_for_each=F){
  nulli=unlist(lapply(angle, is.null))
  # names(angle) = 1:length(angle)
  wp1=lapply(angle[!nulli], function(x1){
#   x1 = apply(x,2,combine_func)
    wp_ = whichpart(x1,n)
    wp_1 = x1[wp_]
    names(wp_1) = wp_
    wp_1
  })
  wp11 = unlist(wp1)
  wp2 = if(one_for_each)  1:length(wp11) else whichpart(wp11,n=n)
  #names(wp2)
  a1 = lapply(names(wp11[wp2]),function(v1){
    split1 = strsplit(v1,"\\.")[[1]]
    c(paste(split1[-length(split1)],collapse="."), split1[length(split1)])
  })
  
  lapply(a1, function(a) c(which(names(angle)==a[1]),as.numeric(a[2])))
}

whichpart <- function(x, n=10) {
  nonNA = which(!is.na(x))
  if(n==1) return(which(x==min(x[nonNA]))[1])
  if(length(nonNA)>n) {
    #nx <- length(x)
    #nacnt = nx -length(nonNA)
    # p <- nx-n-nacnt
    xp <- sort(x, partial=n)[n]
    inds = which(x <= xp)
    if(length(inds)<n){
      #print(xp)
      #print(paste(nx,nacnt,p,n))
      ##print(head(xp))
      ##print(head(x))
      #print(inds)
      stop("problem in which part")
    }
    inds = inds[1:n]
  }else{
    inds = nonNA	
  }
  inds[order(x[inds])]
  
}



.keepUniq<-function(vars, removeTrans=T){
  if(removeTrans){
    vals1 = unlist(lapply(names(vars), function(st)paste(sort(
      unlist(lapply(strsplit(st,",")[[1]], function(st1) paste(strsplit(st1,"\\.")[[1]][-1], collapse="."))))
      ,collapse=",")))
  }else{
    vals1 = unlist(lapply(names(vars), function(st)paste(sort( strsplit(st,",")[[1]])
                                                         ,collapse=",")))
  }
  vars[!duplicated(vals1)]
}


.findMinRMSV<-function(rmsv_, mult=rep(1, length(levels(rmsv_$data))), fspls.sum=T){
  
  if(fspls.sum){
    beams = levels(rmsv_$beam)
    names(beams)=beams
    rmsv_allp = lapply(beams, function(p1){
      rmsv_sub=subset(rmsv_, beam==p1)
      rmsv_sub1  =pivot_wider(rmsv_sub,id_cols="phens", values_from="value", names_from="data")
      rmsv=as.matrix(rmsv_sub1[,-1,drop=F])
      dimnames(rmsv) [[1]] = rmsv_sub1$phens
      rmsv[is.na(rmsv)]=0 ## set NA to zero for multiplication
      #rmsv1 =   (rmsv %*% mult)[,1]
      #attr(rmsv1,"phen")=p1
      #sum(rmsv1)
      rmsv[,1]
    })
    ord=order(unlist(rmsv_allp))
    r1=unlist(rmsv_allp[ord])
    return(r1)
    #rmsv_allp[which(min(unlist(rmsv_allp)))]
  }else{
    phens1 = levels(rmsv_$phens)
    names(phens1) = phens1
    rmsv_allp = lapply(phens1, function(p1){
      rmsv_sub=subset(rmsv_, phens==p1)
      rmsv_sub1  =pivot_wider(rmsv_sub,id_cols="beam", values_from="value", names_from="data")
      rmsv=as.matrix(rmsv_sub1[,-1,drop=F])
      dimnames(rmsv) [[1]] = rmsv_sub1$beam
      rmsv1 =   (rmsv %*% mult)[,1]
      attr(rmsv1,"phen")=p1
      rmsv1
    })
    rmsv_allp[[which.min(unlist(lapply(rmsv_allp, function(v) min(v))))[1]]]
  }
}
.replot_dist<-function(rdsfile,  min=1,max=4, outpdf1=NULL, typed="datas"){
  fi = grep("rds",dir(rdsfile, full=T),v=T)
  if(length(fi)<min) return(NULL)
  names(fi) = unlist(lapply(fi, function(st)gsub(".rds","",rev(strsplit(st,"/")[[1]])[1])))
  fi = fi[order(as.numeric(names(fi)))]
  fi = fi[as.numeric(names(fi))>0]
  if(length(fi)<min) return(NULL)
  fi = fi[1:min(length(fi),max)]
  
  df0 = .merge1(lapply(fi, function(fi1){
    ar = readRDS(fi1)
    .merge1(lapply(ar$datas[[typed]], function(d){
      y1=d$y[,1]
      yp =d$ypreds_all$ypreds[[1]][[1]][,1]
      roc1 =roc(y1,yp,plot=F)
      df2=data.frame(cbind(roc1$sensitivities, roc1$specificities))
      names(df2) = c("sens","spec")
      nme = names(ar$datas[[1]]$train[[1]]$prev)[1]
      print(paste(nme, roc1$auc))
      nme=gsub(",","\n",nme)
      title=paste(nme, collapse="\n")
      title=nme
     
      cbind(df2,title)
      #  if(type=="area") return(plotAreas(yp,y1,title=nme)) else return(ggroc(roc1))
    }), num_cols = c("sens","spec"), addName="lineage")
  }), num_cols = c("sens","spec"), addName="index")
  
  
  df = .merge1(lapply(fi, function(fi1){
    print(fi1)
    ar = readRDS(fi1)
    .merge1(lapply(ar[[typed]], function(d){
      y1=d$y[,1]
      yp =d$ypreds_all$ypreds[[1]][[1]][,1]
      df2=data.frame(cbind(yp,y1))
      names(df2) = c("value","pheno")
      tab = table(df2)
      df3=data.frame(cbind(as.numeric(dimnames(tab)[[1]]), tab))
      
      names(df3) =c("knots","0","1")
     
      df4 = pivot_longer(df3,names(df3[-1]))
      
      names(df4) = c("knots", "pheno", "counts")
      df4=df4[df4$counts>0,,drop=F]
      #ggplot(df4, aes(x=knots, y=value))+geom_point()
      #ggplot(df4, aes(x=knots, y=name, size=value))+geom_point()
      #ggplot(df2, aes(x=factor(pheno), y=value))+geom_jitter(alpha = 0.9, width=0.1, size=.2)
      roc1 =roc(y1,yp,plot=F)
      nme = names(ar$datas[[1]]$train[[1]]$prev)[1]
      print(paste(nme, roc1$auc))
      nme=gsub(",","\n",nme)
      title=paste(nme, collapse="\n")
      title=nme
#      cbind(df2,title)
      cbind(df4,title)
      #  if(type=="area") return(plotAreas(yp,y1,title=nme)) else return(ggroc(roc1))
    }), num_cols = c("knots","value","counts"), addName="lineage")
  }), num_cols = c("knots","value","counts"), addName="index")
 # df$label = apply(cbind(as.character(df$pheno), round(df$value,2)),1,paste,collapse=" ")
  df = df[!is.na(df$pheno),]
  df$pheno=as.numeric(as.character(df$pheno))
  ggp1=ggplot(df, aes(x=knots, y=pheno))#+geom_jitter(alpha = 0.9, width=0.1, size=.2)+ggtitle(rdsfile)#+geom_violin(aes(x=pheno,y=value))
  ggp1<-ggp1+geom_point(aes(size=counts),color="grey")
  ggp1<-ggp1+geom_text_repel(aes(x=knots,y=pheno, label=counts ),size=2)
  if(!is.null(df0)){
  ggp1<-ggp1+geom_line(data=df0, aes(x=1-spec,y=sens), linetype='dashed',colour="grey")
  }
  #ggp1=ggplot(df, aes(x=knots, y=value, color=pheno, linetype=type))+geom_line()+geom_point(aes(size=counts))
#  ggp1= ggp1+geom_text(data=txt_df,nudge_y=0.02,
#                       inherit.aes =T, 
#                       aes(x=knots,y=value,label=label,color=pheno),size=2)
  ggp1<-ggp1+facet_grid(lineage~title)+ggtitle(rdsfile)+scale_x_continuous(limits = c(-0.1,1.1))+scale_y_continuous(limits = c(-0.1,1.1))
  ggp1
if(!is.null(outpdf1))  try(ggsave(outpdf1, plot=ggp1, width =45, height =45, units = "cm",limitsize=F))
  #attr(ggp1,"text_df")=txt_df
  invisible(ggp1)
}

.replot<-function(ar,  typed="datas"){
    #ar = readRDS(rdsfile)
    df=.merge1(lapply(ar[[typed]], function(d){
      y1=d$y[,1]
      yp =d$ypreds_all$ypreds[[1]][[1]][,1]
      roc1 =roc(y1,yp,plot=F)
      nme = names(ar$datas[[1]]$train[[1]]$prev)[1]
      print(paste(nme, roc1$auc))
      nme=gsub(",","\n",nme)
      title=paste(nme, collapse="\n")
      title=nme
      area_plot=getAreaPlot(yp, y1,nme)#, title=title)
      area_plot = addExtraLines(area_plot,yp,y1,nme)
      mat = d$train[[length(d$train)]]$prev[[1]]$betas[[1]]
      if(is.matrix(mat) ){
        beta= try(mat[nrow(mat),1])
        beta = sign(beta)*log10(abs(beta))
        area_plot1 = data.frame(rbind(c(0,beta,0,"beta","summary",nme), c(1,beta,0,"beta","summary",nme)))
        names(area_plot1) = names(area_plot)
        area_plot=rbind(area_plot, area_plot1)
      }
      area_plot
      #  if(type=="area") return(plotAreas(yp,y1,title=nme)) else return(ggroc(roc1))
    }), num_cols = c("knots","value"), addName="lineage")
  df$label = apply(cbind(as.character(df$pheno), round(df$value,2)),1,paste,collapse=" ")
  beta_inds = df$pheno=='beta'
  beta_scale = max(abs(df$value[beta_inds]))
  df$value[beta_inds] = df$value[beta_inds]/(2*beta_scale) + .5
 df
}
.replot_plot<-function(df, title="", outpdf1=NULL){
  txt_df = subset(df ,type=="summary" & knots==0)
  txt_df$knots=0.02
  txt_df$knots[txt_df$pheno=="beta"] = 0.95
  txt_df$label1 = unlist(lapply(as.character(txt_df$label), function(st)strsplit(st," ")[[1]][1]))
  df$label = as.character(df$label)
  ggp1=ggplot(df, aes(x=knots, y=value, color=pheno, linetype=type))+geom_line()+geom_point(aes(size=counts))+ggtitle(title)
  if(!is.infinite(beta_scale)){
  ggp1=ggp1+scale_y_continuous(
    "value", 
    sec.axis = sec_axis(~ (.+.5)*(2*beta_scale), name = "betas")
  )  
  }
  ggp1= ggp1+geom_text(data=txt_df,nudge_y=0.02,
                       inherit.aes =T, 
                       aes(x=knots,y=value,label=label,color=pheno),size=2)
  ggp1<-ggp1+facet_grid(lineage~drug)
  attr(ggp1,"text_df")=txt_df
  if(!is.null(outpdf1))  try(ggsave(outpdf1, plot=ggp1, width =45, height =45, units = "cm",limitsize=F))
  
  ggp2=ggplot(txt_df, aes(x=drug, y=value, color=lineage))+geom_point()+facet_grid(label1~.)
  
  list(ggp1 = ggp1,ggp2=ggp2)
}


addExtraLines<-function(df,yp,y1,title="",input=list()){
  auc = (.calcAUCW(yp,y1))[2]
  youden = (.youden(as.matrix(yp),y1)[2])
  area = attr(df,"area")
  #area =  (.areaBetween(yp,y1))[[2]]
  
  extra = c(auc, youden,area)
  nme_e = c("auc","youden","area")
  extra = data.frame(rbind(cbind(extra,  0), cbind(extra,1)))
  
  extra = cbind(c(nme_e,nme_e), extra)
  extra = cbind(extra, apply(extra[,c(1,2)],1,paste,collapse="="))
  names(extra) = c("pheno","value","knots","label")
  
  extra = cbind(extra,"summary")
  names(extra)[ncol(extra)]="type"
  extra = cbind(extra,0)
  names(extra)[ncol(extra)]="counts"
  
  df=rbind(df, extra[,match(names(df), names(extra))])
  # print(paste(i1,length(kn)))
  cbind(df,title)
}
getAreaPlot<-function(yp, y1,title = "", input = list()){
  levs = sort( unique(y1[!is.na(y1)]))
  names(levs)=levs
  knots=sort(unique(unlist(lapply(levs, function(t){
    inds = y1==t
    cdf=ecdf(yp[inds])
    kn=stats::knots(cdf)
  }))))
  if(!is.null(input$range)){
    knots = c(input$range[0], knots, input$range[1])
  }else{
    knots = c(0,knots,1)
  }
  df=.merge1(lapply(levs, function(t){
    inds = y1==t
    tab1 = table(yp[inds])
    cdf=ecdf(yp[inds])
    mi2=match(knots, names(tab1))
    no_dupl = !duplicated(mi2)
    #mi2 = mi2[!duplicated(mi2)]
    tab1_col = tab1[mi2[no_dupl]]
    tab1_col[is.na(tab1_col)]=0
    df1 = data.frame(cbind(knots[no_dupl],cdf(knots[no_dupl]), tab1_col))
    names(df1)= c("knots","value","counts")
    df1
  }), addName="pheno",num_cols=c("knots","value"))
  pw = pivot_wider(df, names_from="pheno", id_cols="knots")
  diff=0
  for(k in 2:nrow(pw)){
    delta = (pw[k,1]-pw[k-1,1])
    diff = diff+((pw[k-1,2] + pw[k,2] - (pw[k-1,3] + pw[k,3]))/2) * delta
  }
  area=diff
  
  df = cbind(df, "cumulative")
  
  
  
  names(df)[ncol(df)] = "type"
 
 
  attr(df,"area")=area
  df
}
.combinePv<-function(pvs){
  return(min(pvs))
  chisq = qchisq(pvs, df=1, lower.tail=F)
  pchisq(sum(chisq), df = length(pvs), lower.tail=F)
}

.cntNA<-function(vec) length(which(is.na(vec)))

.rep<-function(v,n) t(matrix(rep(v,n),ncol=n))


.getRMSPrev=function(rmsv_){
  beam_levs = levels(rmsv_$beam)
  names(beam_levs) = beam_levs
  rmsv1=unlist(lapply(beam_levs, function(beam){
    ab=subset(rmsv_, beam==beam)
    sum(ab$value,na.rm=T)  #could also try min ? 
  }))
  best_beam = names(which.min(rmsv1))[1]
  rmsv2 = subset(rmsv_, beam==best_beam)
  rmsv_prev1 = rmsv2$value
  names(rmsv_prev1) = apply(rmsv2[,names(rmsv2) %in% c("subpheno","pheno","data")],1,paste,collapse=".")
  rmsv_prev1
}


.addColumn<-function(phenot,json1){
  json_nme = names(json1)
  names(json_nme) = json_nme
  res_all=lapply(json_nme,function(nme){
    disease_class = phenot[[nme]]
    res=data.frame(lapply(json1[[nme]],function(tblkk){
#      tblkk = tbl[kk,]
       dc = rep(NA, length(disease_class))
       dc[(disease_class %in% tblkk)] = disease_class[disease_class %in% tblkk]
       tbls1 = names(sort(table(dc), decr=T))
       factor(dc, levels = tbls1)
      #factor(bacviral, levels = tblkk)
    }))
    names(res) = names(json1[[nme]])
    res
  })
  res_all1 = data.frame(unlist(res_all,rec=F))
  cbind(phenot, res_all1)
}
.inferFamily<-function(y){
  unlist(lapply(y, function(y1) {
    ty = table(y1)
    if(is.factor(y1)) {
      return(if(length(ty)>2) "multinomial" else "binomial")
    }else if(length(ty)==2 && names(ty)[[1]]=="0" && names(ty)[[2]] == "1"){
     return("binomial")      
    }else if( max(apply(cbind(y1,round(y1) ),1,function(x) abs(diff(x))),na.rm=T)==0){
      return("ordinal")
    }else{
      return("gaussian")
    }
  }))
}
.getSubset<-function(phenotypes, col){
  levs = levels(phenotypes[[col]])
  names(levs) = levs
  l1=lapply(levs,function(lev){
    r1 = (phenotypes[[col]]==lev)
    r1[is.na(r1)]=FALSE
    r1
  })
  l1[!unlist(lapply(l1, function(xx) length(which(xx))))==0]
}
.subset<-function(dat, subset){
  dat1 = dat$clone(deep=T)
  dat1$subset = subset
  dat1
}


.appendPredictedPheno<-function(phenotypes, data, nmes =c("Bacterial","Viral","NonInfectious") , thresh =0.5){
  ypred_all1 = data$ypreds_all$ypreds[[1]][[1]]
  dimnames(ypred_all1)[[1]] = dimnames(data$y)[[1]]
  names(ypred_all1 ) = nmes

  disease_class_predicted = 
    factor(apply(ypred_all1,1,function(xx) names(ypred_all1)[which.max(xx)]), levels = names(ypred_all1))
  disease_class_predicted [apply(ypred_all1,1, max)<thresh] = NA
  
 print( table(cbind(disease_class_predicted),data$y[,1]))
  phenotypes[['disease_class_predicted']] = disease_class_predicted
#  return(cbind(phenotypes, disease_class_predicted))
  phenotypes   
}


.splitByPheno<-function(datas, phenotypes, nme='disease_class_predicted'){
  subsets = .getSubset(phenotypes, nme)
  datas_=unlist(lapply(subsets, function(subset){
    lapply(datas, function(d) {
      .subset(d, subset)
    })
  }),rec=F)
  c(datas, datas_)  
}
#datas1 is independent validation
.calcRMSVAll<-function(datas,datas1,cv, label){
  if(!is.null(datas1)){
    mi1 = match(names(datas), names(datas1))
    if(length(which(!is.na(mi1)))>0){
    datas = datas[!is.na(mi1)]
    datas1 = datas1[mi1[!is.na(mi1)]]
    print(cbind(names(datas), names(datas1)))
    }
  }
  
#  cnts_all=unlist(lapply(models, function(m) (m$cnt)))
#  print(table(unlist(lapply(models, function(m) (names(m$prev))))))
#  if(length(models)==1){
#    return(NULL)
#  }
  numv = -1 #c(0:min(cnts_all),-1)
  names(numv) = numv
  names(numv)[which(numv<0)]="full"
  if(cv){
  rmsval = .merge1(lapply( numv, function(kk){
    #kk=0 means evaluate the full models
    rms_cv = .merge1(lapply(datas,function(d) d$getRMSVAll(numvar=if(kk<0) NULL else kk)),num_cols="value",addName="dataset")
    rms_cv
  }),num_cols="value",addName="numvars")
 
  }else{
      k = length(datas[[1]]$train)
      data1_inds = 1:length(datas1); names(data1_inds) = names(datas1)## because we kept the full model here
      rmsval =.merge1(lapply( numv, function(kk){
        print(kk)
        #kk=-1 means evaluate the full models
        .merge1(

          lapply(data1_inds, function(ij) {
            ij2 = if(length(datas)==1) 1 else ij
            prev1 = datas[[ij2]]$train[[k]]$prev
            prev2 = if(label=="discovery") prev1 else datas1[[ij]]$translate(prev1)
            datas1[[ij]]$initY()
            datas1[[ij]]$importModel(prev2,
                                                                         numvar=if(kk<0) NULL else kk,
                                                                 datas[[ij2]]$family)
            })
                ,num_cols ="value", addName = "dataset")
      })
      ,num_cols="value",addName="numvars")
      
      
     
     
  }
  rmsval$value[rmsval$measure=="AUC_all"] =  rmsval$value[rmsval$measure=="AUC_all"]-0.5
  return(cbind(rmsval, label))
}

.rewind<-function(models, datas){
  for(k in 1:length(models))models[[k]]$unwind(datas,k) 
}
.trim<-function(models, datas,to_keep , inds_m=1:length(models)){
  for(k in inds_m){
    models[[k]]$keep(to_keep)
    lapply(datas, function(d) d$train[[k]]$keep(to_keep))
  }
}
.updateModels<-function(models, datas, prev,inds_m, useInternalCVAsStopping){
  for(k in inds_m){
    model = models[[k]]
    if(!useInternalCVAsStopping) model$prev_old=model$prev
    model$prev=prev
    model$cnt = length(prev)
    lapply(datas, function(d) d$updateModel1(k,prev, useOld = useInternalCVAsStopping ))
  }
}

.reorder<-function(models, datas, o, inds_m = 1:length(models)){
  for(k in inds_m){
    model = models[[k]]
    model$reorder(o)
    lapply(datas, function(d) d$reorder(o,k))
  }
}
.summarise<-function(rms_cv, dim1="beam",avg=F){
  #dim1 = 'pheno'
  select1 =   apply(rms_cv[,which(names(rms_cv) %in% dim1),drop=F],1,paste, collapse=".")
  select=factor(select1, lev = select1[!duplicated(select1)])
  
  levs = levels(select)
  .merge1(lapply(levs, function(l1){
    rms_cv1 =rms_cv[select==l1,,drop=F]
  rms_cv2 = rms_cv1[1,,drop=F]
  rms_cv2$value = if(avg) mean(rms_cv1$value,na.rm=T) else sum(rms_cv1$value,na.rm=T)
  rms_cv2$subpheno="comb"
  rms_cv2
  }), num_cols="value")
}
.runAnalysisAll<-function(datas, datas1, params,rdsdir=NULL){
  options(params)
  if(!is.null(rdsdir)) dir.create(rdsdir, rec=T)
  {  #this is setup
      genes_incls=getOption("genes_incls",NULL)
      cv_skips_allowed=getOption("cv_skip",0) ;  signatures = getOption("signatures",NULL)
      incls = getOption("incls",list(names(datas$all$data)))
      ##this is for cross validation evaluation 
      if(length(datas1)>0){
        mi1 = match(names(datas), names(datas1))
        if(length(which(!is.na(mi1)))==0) stop("!!")
        datas = datas[!is.na(mi1)]
        datas1 = datas1[mi1[!is.na(mi1)]]
        # print(cbind(names(datas), names(datas1)))
      }
      num_var =sum(unlist( lapply(datas[[1]]$data, ncol)))
      max_vars = min(num_var,getOption("maxvars",100))
      data1_inds = 1:length(datas1)
      names(data1_inds) = names(datas1)
      types_all = names(datas[[1]]$data)
      names(types_all) = types_all
      var_threshs = lapply(types_all, function(x) getOption("var_thresh",1e-8))
      print("initialising")
      for(ij in 1:length(datas)){print(ij);datas[[ij]]$init1(var_threshs, genes_incls=genes_incls)}
      nreps_all =unlist(lapply(datas, function(d) ncol(d$looc$incl)))
      nreps = table(nreps_all)
      if(length(nreps)>1) {
        ##NEED TO REDO BATCHING
        if(getOption("fspls.batch",0)==0) stop("!!")
        print(paste("new nrep",min(nreps_all)-1 ))
        for(ij in 1:length(datas))datas[[ij]]$init1(var_threshs, nrep = min(nreps_all)-1, batch=0)
        nreps_all =unlist(lapply(datas, function(d) ncol(d$looc$incl)))
        nreps = table(nreps_all)
        if(length(nreps)>1) stop("could not resolve problems with nrep and batch")
      }
      nrep1  = as.numeric(names(nreps))
      nrep=if(nrep1==1)1 else nrep1-1 
      print("building models")
      models=lapply(1:nrep1, function(k) {
        ##following line could be applied on yPreds above
        rmsv_=.merge1(lapply(datas, function(d)d$getRMSV(k)),num_cols="value",addName="data")  #this is call 3 to datas
        rmsv_prev1 =   .getRMSPrev(rmsv_)
        modelObj$new(names(datas), names(datas[[1]]$data),  rmsv_prev1)
        
      })
      names(models) = 1:length(models)
      # jks=1; k=1;
      rms_prev=sum(models[[1]]$rmsv_prev, na.rm=T)
      skipped=0
      useInternalCVAsStopping= getOption("useInternalCVAsStopping",F)
      orderCV = getOption("orderCV",T)
      pv_only=getOption("fspls.pval_only",F)
      beam=getOption("fspls.beam",c(1,1))
      logpthresh=log(getOption("fspls.pthresh1",0.05))
      replaceModel = beam[1] *beam[2] >1  || !useInternalCVAsStopping  ## can only not replace models if no beam
  }
  for(jks in 1:length(incls)){
    print(paste("training", jks))
    for(model in models) model$finished=F
    finished = FALSE
    to_keep = 1  ## which samples to use in next round
    maxvars = getOption("maxvars",100)
    cnt1=1
    maxvars_jks = maxvars[min(jks, length(maxvars))]
    while(!finished && models[[length(models)]]$cnt<sum(maxvars)){
      if(cnt1>maxvars_jks){
        break;
      }else{
        cnt1 = cnt1+1
      }
      incl=incls[[jks]]
      inds_m = if(useInternalCVAsStopping) 1:length(models) else length(models)
      pvslist_all =lapply(inds_m, function(k){
        model = models[[k]]
        nxt_var = if(models[[k]]$cnt<length(signatures) ) signatures[[models[[k]]$cnt+1]] else NULL
        model$simpleForwardTrain(datas,k,  exclude=exclude,nxt_var = nxt_var, incl = incl ,  weights=NULL, to_keep=to_keep)
      })
      pvslist_all1 =pvslist_all[[length(pvslist_all)]]$pv
      pvslist_all2 = pvslist_all[[length(pvslist_all)]]$cumpv
      print(pvslist_all1)
      pvtk1=  pvslist_all1<logpthresh
      if(length(which(pvtk1))==0){
        .trim(models, datas,c(),inds_m  )
        print("breaking on pvalue")
        finished=T
        break;
      }
      
      if(!pv_only &useInternalCVAsStopping & nrep>1){
        rms_cv1 = .merge1(lapply(datas,function(d) .summarise(d$getRMSVAll())),num_cols="value",addName="dataset")
        rms_prev1=min(rms_cv1$value)
        if(rms_prev1>=rms_prev){
          .trim(models, datas,c(),inds_m )
          finished=T
          break;
        }
      }
      vars = names(pvslist_all2)
      prev = models[[length(models)]]$prev
      if(length(models)>1 && replaceModel) {
        .updateModels(models, datas, prev,1:(length(models)-1),useInternalCVAsStopping)
      }
      if(!pv_only){
          if(length(models)==1){
              rms_cv = .merge1(lapply(datas,function(d) .summarise(d$getRMSV(length(models)))),num_cols="value",addName="dataset")
        }else{
              rms_cv = .merge1(lapply(datas,function(d) .summarise(d$getRMSVAll())),num_cols="value",addName="dataset")
          }
          rms_cv = rms_cv[match(vars, rms_cv$beam),,drop=F]
          o = if(orderCV) order(rms_cv$value) else order(pvslist_all2) 
          names(o) = if(orderCV)  rms_cv$beam[o]  else names(pvslist_all2)[o]
      }else{
        o = order(pvslist_all2)
        names(o) = names(pvslist_all2)[o]
      }
         to_keep =o[1:min(length(vars),beam[2])]
         #attr(models,"to_keep")=to_keep
          pvtk=   pvslist_all1[to_keep]<logpthresh 
          if(length(which(pvtk))==0){
            .trim(models, datas,c()  )
            finished=T
            break;
          }
          #if(length(which(!pvtk))>0) pvtk = pvtk[1:which(!pvtk)[1]]  # only keep until first non-sig
          #to_keep = to_keep[pvtk]
#          tokeep1 = rms_cv$beam[to_keep]
          if(!pv_only){
            if(!useInternalCVAsStopping || length(models)==1) rms_prev1=min(rms_cv$value[to_keep],na.rm=T)
            if(rms_prev1>rms_prev){
              skipped = skipped+1
              print(paste("skipped ",skipped))
            }else{
              #reset skipped and reset bar
              print("reset skipped")
              skipped=0
              rms_prev = rms_prev1
            }
          
          if(skipped> cv_skips_allowed){
            if(length(to_keep)>0){ ## need to wind back since rms got worse
              for(k in 1:length(models)){
                models[[k]]$keep(c())
                lapply(datas, function(d) d$train[[k]]$keep(c()))
              }
            }
            finished=T
            print(paste("breaking here because cv detetiorating", models[[1]]$cnt))
            
          #  for(k in which(!finished)){
          #    models[[k]]$unwind(datas,k) 
          #  }
            break;
          }
          }
          if(length(to_keep)==0){
            finished=T
            break
          }
          to_keep =1:min(length(vars),beam[2])  ## really only comes into play if beam>1
          attr(models,"to_keep")=to_keep
          .reorder(models, datas, o)          
#          attr(models,"to_keep")=to_keep
      if(!is.null(rdsdir)){
        outdir1 =  paste(rdsdir,paste(models[[length(models)]]$cnt,"rds",sep="."),sep="/")
        xls_file = paste(rdsdir,paste(models[[length(models)]]$cnt,"xlsx",sep="."),sep="/")
        print(paste("saving data to outdir1", outdir1))
        #sr =  .saveDatas(datas,datas1, params, outdir1,xls_file, types=c("AUC","area","sens_spec"))
        sr =  .saveDatas(datas,datas1, outdir1,xls_file, types=c("AUC","area","youden_sens","youden_spec"))
        if(is.null(sr$rms$validation)){
          if(is.null(sr$rms$crossvalidation)){
            print(head(.summarise(sr$rms$discovery, c("beam","pheno","measure"),avg=T)))
          }else{
            print(head(.summarise(sr$rms$crossvalidation, c("beam","pheno","measure"),avg=T)))
            
          }
        }else{
            print(head(.summarise(sr$rms$validation, c("beam","pheno","measure"),avg=T)))
        }
        
      }
    }
  }
  
  
  invisible(models)
}


.logLik<-function(y,ypred=mean(y), family="binomial"){
  if(family=="binomial"){
    yp1=.logistic(ypred)
    ll1=sum(apply(cbind(y,yp1), 1,function(v) if(v[1]==1) log(v[2]) else log(1-v[2])))
  }else{
    ll1=sum(apply(cbind(y,ypred), 1,function(v) dnorm(v[2], mean = v[1],log=T)))
  }
  
}
.lrt<-function(ll1,ll2,  df1 =  attr(ll1,"df")[1] ,  df2 = attr(ll2,"df")[1],totweight=1,log.p=F){
  if(df2>df1) stop(paste("df2 should be bigger than df1 ",df2,df1))
  pchisq((2*(ll1 - ll2))/(totweight),df1-df2,lower.tail=FALSE,log.p=log.p)
}
.norm1<-function(g){
  sqrt(sum((g - mean(g, na.rm = T))^2, na.rm = T))
}
#.plotAUC<-function(rdsf){
#  
#}


.plotRMS2<-function(tEnvs,type="cv",fam=tEnvs[[1]][[1]]$datas[[1]]$family,
                    measures = NULL){
  #print(tEnv$rms_list_all)
  a7=.merge1(lapply(tEnvs, function(tEnvs1){
    .merge1(lapply(tEnvs1, function(tEnv){
      names(tEnv$rms_list_all)=1:length(tEnv$rms_list_all)
      a1 = .merge1(list(
        cv=.merge1(lapply(tEnv$rms_list_all, function(t) .mod(t$cross)), addName="lens1",num_cols=c("value","lens")),
        val=.merge1(lapply(tEnv$rms_list_all, function(t) .mod(t$val)), addName="lens1",num_cols=c("value","lens")),
        disc=.merge1(lapply(tEnv$rms_list_all, function(t) .mod(t$disc)), addName="lens1",num_cols=c("value","lens"))),
        num_cols = c("value","lens","lens1"),addName="type")
      o = order(a1$dataset)
      #a1[o,]
      a1$value = -1*(a1$value)
      a2 = subset(a1,beam1==1)
      if(fam=="multinomial"){
        return(a2)
      }
      a3 =  a2[grep("\\.mid$",a2$subpheno),,drop=F]
      a4 =  a2[grep("\\.low$",a2$subpheno),,drop=F]
      a5 =  a2[grep("\\.high$",a2$subpheno),,drop=F]
      low = a4$value
      high = a5$value
      a6 = cbind(a3,low, high)
      a6
    }),addName="dataset1",num_cols=c("value","lens","lens1","low","high"))
  }),addName="drug",num_cols=c("value","lens","low","lens1","high"))
  a71 = a7[a7$type==type,,drop=F]
  if(!is.null(measures)){
    a71 = subset(a71, measure %in% measures)
  }
  if(fam %in% c("multinomial","binomial")){
    ggp<-ggplot(a71,aes(x=lens1, y=value, color=subpheno, shape=dataset,linetype=dataset))+geom_line()+geom_point()
    ggp<-ggp+facet_grid(measure ~dataset)
  }else{
    
    ggp<-ggplot(a71,aes(x=lens1, y=value, color=dataset1, linetype=beam1))+geom_line()+geom_point()
    ggp<-ggp+geom_ribbon(aes(ymin=low, ymax=high, fill=dataset),alpha = 0.05)
    ggp<-ggp+facet_grid(measure~drug,scales="free_y")+ggtitle(type) # +ggtitle(paste(a1$beam[a1$len==max(a1$lens)][1],"cross validation",sep="\n"))
  }
  return(ggp)
}



.replot2<-function(tEnvs, type1="cv", type="area", outpdf1=NULL, typed="datas"){
  df = .merge1(lapply(tEnvs, function(tEnvs1){
    .merge1(lapply(tEnvs1, function(tEnv){
      names(tEnv$rms_list_all)=1:length(tEnv$rms_list_all)
      
      .merge1(lapply(tEnv[[typed]], function(d){
        y1=d$y[,1]
        yp =d$ypreds_all$ypreds[[1]][[1]][,1]
        roc1 =roc(y1,yp,plot=F)
        nme = names(ar$datas[[1]]$train[[1]]$prev)[1]
        print(paste(nme, roc1$auc))
        nme=gsub(",","\n",nme)
        title=paste(nme, collapse="\n")
        title=nme
        area_plot=getAreaPlot(yp, y1,nme)#, title=title)
        area_plot = addExtraLines(area_plot,yp,y1,nme)
        mat = d$train[[length(d$train)]]$prev[[1]]$betas[[1]]
        if(is.matrix(mat) ){
          beta= try(mat[nrow(mat),1])
          beta = sign(beta)*log10(abs(beta))
          area_plot1 = data.frame(rbind(c(0,beta,0,"beta","summary",nme), c(1,beta,0,"beta","summary",nme)))
          names(area_plot1) = names(area_plot)
          area_plot=rbind(area_plot, area_plot1)
        }
        area_plot
        #  if(type=="area") return(plotAreas(yp,y1,title=nme)) else return(ggroc(roc1))
      }), num_cols = c("knots","value"), addName="lineage")
    }), num_cols = c("knots","value","counts"), addName="index")
  }),num_cols = c("knots","value","counts"), addName="drug")
  df$label = apply(cbind(as.character(df$pheno), round(df$value,2)),1,paste,collapse=" ")
  beta_inds = df$pheno=='beta'
  beta_scale = max(abs(df$value[beta_inds]))
  df$value[beta_inds] = df$value[beta_inds]/(2*beta_scale) + .5
  txt_df = subset(df ,type=="summary" & knots==0)
  txt_df$knots=0.02
  txt_df$knots[txt_df$pheno=="beta"] = 0.95
  
  ggp1=ggplot(df, aes(x=knots, y=value, color=pheno, linetype=type))+geom_line()+geom_point(aes(size=counts))+ggtitle(rdsfile)
  if(!is.infinite(beta_scale)){
    ggp1=ggp1+scale_y_continuous(
      "value", 
      sec.axis = sec_axis(~ (.+.5)*(2*beta_scale), name = "betas")
    )  
  }
  ggp1= ggp1+geom_text(data=txt_df,nudge_y=0.02,
                       inherit.aes =T, 
                       aes(x=knots,y=value,label=label,color=pheno),size=2)
  ggp1<-ggp1+facet_grid(lineage~drug)
  attr(ggp1,"text_df")=txt_df
  if(!is.null(outpdf1))  try(ggsave(outpdf1, plot=ggp1, width =45, height =45, units = "cm",limitsize=F))
  
  invisible(ggp1)
}
