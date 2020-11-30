.calcPropCI<-function(x, conf.int=0.95,  method="prop.test"){
  if(is.na(x[2])) return (c(NA,NA,NA))
  if(x[1]>x[2]) return (c(NA,NA,NA))
  v = binom.confint(x[1], x[2], method=method, conf.int=conf.int)[4:6] *(1e6)
unlist(v)
}
.getHeaderH5<-function(datafile, toreplace=list()){
  header = h5read(datafile,"header")
  if(length(toreplace)>0){
    for(i in 1:length(toreplace)){  
      header[header==names(toreplace)[i]] =toreplace[[i]]
      
    }
  }
  header
}
.getIsoInfo<-function(datafile, h5file,toreplace=list()){
 
 header=.getHeaderH5(datafile, toreplace)
 mat =h5ls(datafile)
 count_entries= apply(mat[grep("/counts",mat$group),1:2],1,paste,collapse="/")
 print(h5file)
 if(file.exists(h5file)){
   mat1 =h5ls(h5file)
   orf_entries = mat1[grep("/depth",mat1$group),]$name
   
 }else{
  orf_entries = mat[grep("/trans",mat$group),]$name
 }
  counts=lapply(count_entries, function(x) h5read(datafile,x))
  
  total_reads = apply(data.frame(counts),1,sum)
  #if(sum(total_reads)==0){
  #  total_reads = rep(1e6,length(total_reads))
  #}
  names(total_reads) = header
  experiments <- data.frame(experiments=header)
 exps= separate(experiments,1, c('molecule_type', 'cell', 'time'), sep='_', remove = T) %>%
    transform(  molecule_type = factor(molecule_type), cell = factor(cell), time = factor(time)) 
  num_breaks=unlist(lapply(strsplit(orf_entries,","),function(x) length(x)-1))
  vals = list(
     grep("end",grep("start|leader", orf_entries,v=T),v=T),
     grep("end",grep("start|leader", orf_entries,v=T),v=T,inv=T),
     grep("end",grep("start|leader", orf_entries,v=T,inv=T),v=T),
     grep("end",grep("start|leader", orf_entries,v=T,inv=T),v=T,inv=T)
  )
  names(vals) = c("5_3", "5_non3","non5_3","non5_non3") 
 vals_ind= lapply(vals,function(x) which(orf_entries %in% x))
  type_name=rep("NA", length(num_breaks))
  for(i in 1:length(vals_ind)){
    type_name[vals_ind[[i]]] = names(vals_ind)[i]
  }
  
  list(total_reads = total_reads, experiments=exps, 
       orfs=data.frame(ORFs=orf_entries, type_name=factor(type_name),num_breaks = factor(num_breaks) )
                )
}

.processInfo<-function(isoInfo){
  ORFs = isoInfo$orfs
  exps = isoInfo$experiments
 # choices = lapply(levels(ORFs$num_breaks), function(x) sort(as.character(ORFs$ORFs[ORFs$num_breaks==x])))
  choices1 = lapply(levels(ORFs$type_name),function(x) as.character(ORFs$ORFs[ORFs$type_name==x]))
  nmes1 = c("5_3","5_non3","non5_3","non5_non3")
#  choices1= vector("list", length(nmes1))
  cmax = max(as.numeric(levels(ORFs$num_breaks)))
#  c1 = vector("list",cmax+1)
  choices = lapply( 0:cmax, function(x) as.character(ORFs$ORFs[ORFs$num_breaks==x]))
  names(choices) = paste(0:cmax,"junctions")
  if(length(choices1)<4){
    
    choices1 = c(choices1,rep("",4-length(choices1)))
    names(choices1) = nmes1
  }else{
    names(choices1) = levels(ORFs$type_name)
  }
  names(choices) = paste(as.numeric(levels(ORFs$num_breaks)),"junctions")
#  if(length(choices)>4){
#  choices[[4]] = as.character(unlist(choices[-(1:3)]))
  
#}
 # names(choices)[[1]] = "zero junctions"
#  names(choices)[[2]] = "one junction"
#  names(choices)[[3]] = "two or more junctions"
  molecules=levels(exps$molecule_type)
  times = levels(exps$time)
  cells = levels(exps$cell)
  
  
  list(choices=choices, choices1=choices1,molecules=molecules, times=times, cells = cells)
}


.readTPM<-function(datafile){
  transcripts <- .readTranscripts(datafile)
  experiments <- attr(transcripts, 'info')
  experiments
  count_idx <- which(grepl(pattern = 'count[0-9]', x = colnames(transcripts)))
  total_reads = apply(transcripts[,count_idx],2,sum)
  names(total_reads)= experiments
  #if(!is.null(method)){}
   tpm= transcripts[,count_idx]
  
  .processTPM(tpm, experiments, transcripts$ORFs)
}
 .processTPM<-function(mat, experiments, ID,levels =NULL,split=T,
                       toreplace1 = c("leader_ORF1ab,S_ORF1ab,ORF10_end","leader_leader,S_ORF1ab,ORF10_end"),
                       toreplace2= c("leader_ORF1ab,ORF1ab_ORF1ab,ORF10_end","leader_leader,ORF1ab_ORF1ab,ORF10_end")
                       ){ 
   toplot=ID
   inds = if(is.null(levels)) 1:dim(mat)[2] else which(experiments %in% levels)
   tpm = mat[,inds,drop=F]
   experiments = experiments[inds]
   colnames(tpm) <- experiments
  #props <- prop.table(x = as.matrix(transcripts[,count_idx]), margin = 2)
  # tpm=props*1e6
  #for use in 'split_by', if I can get it to work
  time_vec <- c('0hpi','2hpi','24hpi','48hpi')
  exp_vec <- c('virion','vero','calu','caco')
  #calculate tpm
  #tpm <- props*1e6
  tpm <- cbind(ID=ID,tpm)
  #prep tpm_df
  #separate(TPM1, c('TPM', 'lower', 'upper'), sep=',', remove = T) %>%
  if(split){
  as.data.frame(tpm) %>% melt(id.vars='ID', measure.vars=experiments, value.name = 'count') %>%
    separate(variable, c('molecule_type', 'cell', 'time'), sep='_', remove = T) %>%
    transform( count=as.numeric(count), molecule_type = factor(molecule_type), cell = factor(cell), time = factor(time, levels = time_vec)) -> tpm_df
  }else{
    as.data.frame(tpm) %>% melt(id.vars='ID', measure.vars=experiments, value.name = 'count', variable.name='sample') %>%
      transform(count=as.numeric(count),ID=factor(ID, levels=toplot),sample = factor(sample, levels =levels)) -> tpm_df
    
  } #attr(tpm_df,"total_reads") = total_reads
  if(length(toreplace1)>0){
        tpm_df$ID= factor(.fix(as.character(tpm_df$ID),toreplace1,toreplace2), levels = .fix(toplot,toreplace1,toreplace2))
  }
  tpm_df
 }

  .fix<-function(v ,nme,toreplace){
   for(i in 1:length(toreplace)){
    v[v== nme[i]] = toreplace[[i]]
   }
    .inner<-function(x){
      paste(unique(strsplit(x,"_")[[1]]),collapse="_")
      
    }
   v2 =  lapply(v,function(x) strsplit(x,",")[[1]])
   for(i in 1:length(v2))v2[[i]] =paste( unlist(lapply(v2[[i]], .inner)),collapse=",")
   unlist(v2)
 }
 
.readIso<-function(x, isofile,header, group="/trans"){
  len  = length(header)
  cnts = h5read(isofile, paste(group, x,sep="/"))
  c(cnts, rep(0,len -length(cnts)))
}

.mergeGroups<-function(matname,merge_by_all, ord = NULL, max_trans=10){
  x=matname
  merge_all =strsplit(merge_by_all,"\\|")[[1]]
  tomerge = lapply(merge_all,.findEntries1, matname)
  names(tomerge) = merge_all
 
    indsx = !(x %in% unlist(tomerge))
    if(!is.null(ord))  ord = ord[indsx]
    others = x[indsx]
  if(!is.null(ord)){
    if(length(others)>max_trans){
      others = others[ord[1:max_trans]]
    }
  }
  
  names(others) = others
  l = c(tomerge,as.list(others))
  l
}

.getGroups<-function(x1, group_bys){
  group_l = unlist(strsplit(group_bys,":")[[1]])
  l1 = list(x1)
  for(i in 1:length(group_l)){
    l1 = unlist(lapply(l1,.getGroupsInner,group_l[i]),recursive=F)
  }
  l1
}

.getGroupsInner<-function(x1,group_by){
  l = list()
  if(group_by=="all"){
    l = list("all"=x1)    
  }else if(group_by=="type"){
    l = list(
      grep("end",grep("start|leader",x1,v=T),v=T),
      grep("end",grep("start|leader",x1,v=T),v=T,inv=T),
      grep("end",grep("start|leader", x1,v=T,inv=T),v=T),
      grep("end",grep("start|leader", x1,v=T,inv=T),v=T,inv=T)
    )
    #l = vals #lapply(vals,function(x) which(x1 %in% x))
    names(l) = c("5_3", "5_non3","non5_3","non5_non3") 
  }else if(group_by=="juncts"){
    juncts = factor( unlist(lapply(x1,function(x)-1+length(strsplit(x,",")[[1]]))))
    junctlev = levels(juncts)
    #  print(juncts)
    #  print(junctlev)
    l = list()
    for(k in 1:length(junctlev)){
      l[[k]] = x1[which(juncts==junctlev[k])]
    }
    names(l) = junctlev
  }else{
    l[[1]] = grep(group_by,x1,v=T)
    l[[2]] = grep(group_by, x1,inv=T,v=T )
    names(l) = c(group_by,paste("!",group_by))
  }
  l[ unlist(lapply(l, length))>0]
}





  
.subsetFCFile<-function(mat, plot_params){
  print("getting subset")
  print(plot_params)
  toplot=plot_params$toplot5
  toplot = toplot[toplot!="-"]
  matname = dimnames(mat)[[1]]
  toplot2 = plot_params$toplot2
  toplot2 = toplot2[unlist(lapply(toplot2,nchar))>2]
  group_by = plot_params$group_by
  merge_by = plot_params$merge_by
  if(length(toplot)==0 ){
    if(length(toplot2)==0) toplot2 = "all"
    toplot = toplot2
    if(toplot!="all"){
      mat= mat[matname %in% .findEntries1(toincl,matname,tojoin=plot_params$tojoin),,drop=F]
    }
    if(merge_by!=""){
      ncol = dim(mat)[[2]]
      print("MERGING!!!")
      print(merge_by)
      
      grps=.mergeGroups(matname,merge_by)
    #  print(grps[[1]])
    #  print(names(grps))
      mat1 = matrix(nrow = length(grps), ncol = ncol)
      for(i in 1:length(grps)){
        mat1[i,] = apply(mat[grps[[i]],,drop=F],2,sum)
      }
      dimnames(mat1) = list(names(grps), dimnames(mat)[[2]])
      mat = mat1
    }else if(group_by != "No grouping"){
      ncol = dim(mat)[[2]]
      grps = .getGroups(matname,group_by)
      mat1 = matrix(nrow = length(grps), ncol = ncol)
      for(i in 1:length(grps)){
        mat1[i,] = apply(mat[grps[[i]],,drop=F],2,sum)
      }
      dimnames(mat1) = list(names(grps), dimnames(mat)[[2]])
      mat = mat1
    }
  }else if(length(toplot>0)){
    mat = mat[matname %in% toplot,,drop=F]
  }
 # print(dimnames(mat))
  mat
}

.findEntries<-function(x,isofile, group ,tojoin="OR", merge=F){
  mat = h5ls(isofile)
  mat = mat[mat$group==group,,drop=F]
  .findEntries1(x,mat$name, tojoin) 
}

.findEntries1<-function(x,matname ,tojoin="OR"){
  print(paste("tojoin",tojoin))
  join_and = tojoin=="AND"
  join_not = tojoin=="AND NOT"
if(join_and || join_not)  x2 =  matname else   x2 =  c()
  for(j in 1:length(x)){
    
    if(x[j]=="all"){
      x1 = matname
    }else if(x[j]=="non3"){
      x1=grep("end",matname,inv=T,v=T)
    }else if(x[j]=="non5"){
      x1=grep("leader",matname,inv=T,v=T)
    }else if(x[j]=="non5_non3"){
      x1 = grep("end",grep("leader",matname,inv=T,v=T),inv=T,v=T)
    }else if(x[j]=="non5_3"){
      x1 = grep("end",grep("leader",matname,inv=T,v=T),inv=F,v=T)
    }else if(x[j]=="5_non3"){
      x1 = grep("end",grep("leader",matname,inv=F,v=T),inv=T,v=T)
    }else if(x[j]=="5_3"){
      x1 = grep("end",grep("leader",matname,inv=F,v=T),inv=F,v=T)
    }else if(length(grep("juncts",x[j]))>0){
      if(length(grep("=",x[j]))>0){
         num = as.numeric(strsplit(x[j],"=")[[1]][2])
        x1=matname[unlist(lapply(strsplit(matname,","),length))==(num+1)]
      }else if(length(grep("<",x[j]))>0){
        num = as.numeric(strsplit(x[j],"<")[[1]][2])
        x1=matname[unlist(lapply(strsplit(matname,","),length))<(num+1)]
      }else if(length(grep(">",x[j]))>0){
        num = as.numeric(strsplit(x[j],">")[[1]][2])
        x1=matname[unlist(lapply(strsplit(matname,","),length))>(num+1)]
      }
    }else{
     x1=grep(x[j],matname,v=T)
    }
    if(join_and || (join_not && j=1)) x2 = x2[x2 %in% x1] else if (join_not) x2 = x2[!(x2 %in% x1)] else x2 = c(x2, x1[!(x1 %in% x2)])
}
  x2
}
#.readIsoGrep<-function(x, isofile,header, group="/trans"){
# x1 = .findEntries(x,isofile,group);
#  mat1 = t(data.frame( lapply(x1, .readIso, isofile, header, group)))
#  mat2 = matrix(apply(mat1,2,sum),nrow=1,ncol=dim(mat1)[2])
#  mat2
#}
.readTotalIso<-function(isofile, group="/trans", trans=NULL){
  names = h5ls(isofile)
  isoheader = h5read(isofile,"header")
  #total_reads = apply(allcnts,2,sum)
  if(is.null(trans))  trans = names[grep("/trans",names$group),]$name
  allcnts = t(data.frame(lapply(trans, .readIso, isofile, isoheader)))
  dimnames(allcnts) = list(trans, isoheader)
  cnts=apply(allcnts,1,sum)
  names(cnts) = trans
  order(cnts, decreasing=T)
}


getKmer<-function(base, pos,v = c(-1,0,1)){
  res = matrix("-", nrow = length(pos), ncol = length(v))
  for(k in 1:length(v)){
    
    
    vect = pos+v[k]
    inds = vect>0 & vect <= length(base)
    res[inds,k]  = base[vect[inds]]
    
  }
  apply(res,1,paste, collapse="")
}
.mergeMatr<-function(cols, cols1){
  pos_c = unique(sort(c(cols$pos, cols1$pos)))
  cols2 = data.frame(array(dim = c(length(pos_c), 4)))
  names(cols2) = names(cols)
  cols2$pos = pos_c
  cols2$s_e = as.factor(rep("end", length(pos_c)))
  cols2$type = as.factor(rep(1, length(pos_c)))
  cols2$depth = rep(0, length(pos_c))
  inds1 = match(cols1$pos, pos_c)
  cols2$depth[inds1] = cols2$depth[inds1]  + cols1$depth 
  inds = match(cols$pos, pos_c)
  cols2$depth[inds] = cols2$depth[inds]  + cols$depth 
  attr(cols2,"inds") = inds
  attr(cols2,"inds1") = inds1
  cols2
}

.mergeBreak<-function(brP){
  heatm = brP[[1]]$heatm
  rows = brP[[1]]$rows
  cols = brP[[1]]$cols
  if(length(brP)>1){
    for(i in 2:length(brP)){
       heatm1 = brP[[i]]$heatm
       rows1 = brP[[i]]$rows
       cols1 = brP[[i]]$cols
      cols2 = .mergeMatr(cols, cols1)
      rows2 = .mergeMatr(rows, rows1)
       heatm2 = array(0,dim = c(dim(rows2)[[1]], dim(cols2)[[1]]), dimnames = list(rows2$pos, cols2$pos))
      heatm2[attr(rows2,"inds"), attr(cols2,"inds")] =  heatm2[attr(rows2,"inds"), attr(cols2,"inds")] + heatm 
      heatm2[attr(rows2,"inds1"), attr(cols2,"inds1")] =  heatm2[attr(rows2,"inds1"), attr(cols2,"inds1")] + heatm1 
      heatm = heatm2
      rows = rows2
      cols=cols2
    }
  }
  list(heatm =heatm, rows = rows, cols = cols)
}

.readAnnotFile<-function(fi, total_reads, norm=F , conf.level=0.95,  levels=NULL, 
                         orfs="E,N"){
  annot = read.table(fi, head=T)
  type_name = names(total_reads)
  type_nme= names(total_reads)
if(is.null(levels)){
  levels=type_nme
}
  #if(!norm){
  #  total_reads = rep(1,length(total_reads))
  #  names(total_reads) = type_name
  #}else{
  #  total_reads = total_reads/1e6
  #}
  toinclude=which(type_nme %in% levels)
  spi  = grep("Spliced", names(annot))
  uspi = grep("Unspliced", names(annot))
  ratio = data.frame(matrix(nrow = dim(annot)[1], ncol = length(spi)*3))
  mixt = data.frame(matrix(nrow = dim(annot)[1], ncol = length(spi)*3))
  
 nme_r = c("ORF", "Start", "Ratio","lower","upper", "type","spliced","unspliced","logdiff",
           "logtotal")
  
  offset = 0
  if(is.null(type_nme)) type_nme = 1:length(spi)
  nme = c()
  if(is.null(toinclude)) toinclude = 1:length(spi)
  ratio1 = data.frame(matrix(nrow = dim(annot)[1]*length(toinclude), ncol = length(nme_r)))
  
  names(ratio1) = nme_r
  for(i in toinclude){
    sumr = annot[,spi[i]] + annot[,uspi[i]]
    confints =   binom.confint(annot[,spi[i]],sumr ,method="logit",conf.level=conf.level)
  
    confints1 = confints#  binom.confint(sumr ,rep(total_reads[i], length(sumr)), method="logit",conf.level=conf.level)
    
    #confints[,4:6]*sumr
    nme = c(nme,paste(c("mean","lower","upper"), type_nme[i],sep="_"))
   # print(nme)
    indsi = ((i-1)*3+1):((i)*3)
  #  ratio[,indsi ]= confints[,4:6]
    ranges = offset + 1:length(annot$Gene)
    ratio1[ranges,1] = as.character(annot$Gene)
    ratio1[ranges,2] = as.numeric(as.character(annot$Start))
    ratio1[ranges,3:5] = confints[,4:6]
    ratio1[ranges,6] = rep(type_nme[i], length(ranges))
    ratio1[ranges,7] = annot[,spi[i]]; #/total_reads[i]
    ratio1[ranges,8] =  annot[,uspi[i]];#/total_reads[i]
  
    ratio1[ranges,9] = apply(-log2(confints[,4:6]),1,paste,collapse=":")
    ratio1[ranges,10] = apply(log2(confints[,4:6]*1e6),1,paste,collapse=":")
    #ratio1[ranges,7:9] = mixt[,indsi]
    offset = offset + length(ranges)
  }
  #names(ratio) =nme
 
  orfs_include=unlist(strsplit(orfs,","))
  ratio1 = ratio1[ratio1$ORF %in% orfs_include,,drop=F]
 ratio1
}
#  ggp= NULL
.plotAnnotFile<-function(ratio1, levels=NULL,barchart=F,showEB = F,facet = NULL, showSecondAxis=F,coeff=1,diff=0, y_text="Ratio", size=20,textsize=20, linesize=1.5){
   if(!barchart){
    timevec=c("0hpi","2hpi","24hpi","48hpi")
    ratio3= separate(ratio1,6, c('molecule_type', 'cell', 'time'), sep='_', remove = T) %>%
      transform(  molecule_type = factor(molecule_type), cell = factor(cell), time = factor(time,levels=timevec))
    
    ratio5=
      separate(ratio3,logdiff, c('value', 'lower', 'upper'), sep=':', remove = T) %>%
      transform( lower=as.numeric(lower),upper=as.numeric(upper),value=as.numeric(value)  ) 
    ratio5 = ratio5[!( is.infinite((ratio5$value)) | is.na(ratio5$value)),]
    
    # if(is.null(orfs_include)) orfs_include = levels(factor(ratio3$ORF))
    lims = c(0,max(ratio5$value[is.finite(ratio5$value)],na.rm=T))
    
    if(!showSecondAxis){
      ggp1<-ggplot(ratio5, aes(x=time))
      ggp1<-ggp1+geom_line(position=position_dodge(width=0.1),aes(y=value ,group=interaction(molecule_type, cell, ORF), color = cell), size=linesize)
      
    if(showEB) {
      ggp1<-ggp1+ geom_errorbar(aes(ymin=lower, ymax=upper, group=interaction(molecule_type, cell, ORF), color = cell),
                                position=position_dodge(width=0.1)) #,colour="black")
      ggp1<-ggp1+geom_point(position=position_dodge(width=0.1),aes(y=value ,group=interaction(molecule_type, cell, ORF), color = cell, shape=ORF,size=10))
      
     # ggp1<-ggp1+geom_ribbon(aes(ymin=lower, ymax=upper, group=interaction(molecule_type, cell, ORF), color = cell), linetype=1, alpha=0.1)
    }else{
      ggp1<-ggp1+geom_point(aes(y=value ,group=interaction(molecule_type, cell, ORF), color = cell, shape=ORF,size=10))
      
    }
      ggp1<-ggp1+ scale_y_continuous( name = "Log2 (total - sub-genomics)", limits=lims)
      data  = ratio5
        }else{
  #    ratio3$logtotal = (ratio3$logtotal-diff)/coeff
      ratio4 = melt(ratio3,id.vars=c("ORF","molecule_type","cell","time"), measure.vars=c("logdiff","logtotal")) %>%
        separate(value, c('value', 'lower', 'upper'), sep=':', remove = T) %>%
        transform( lower=as.numeric(lower),upper=as.numeric(upper),value=as.numeric(value)  ) 
     ratio4 = ratio4[!( is.infinite((ratio4$value)) | is.na(ratio4$value)),]
      indst = ratio4$variable=="logtotal"
      ratio4$value[indst] = (ratio4$value[indst]-diff)/coeff
      ratio4$lower[indst] = (ratio4$lower[indst]-diff)/coeff
      ratio4$upper[indst] = (ratio4$upper[indst]-diff)/coeff
      
      data=ratio4
      #print(head(data))
      ggp1<-ggplot(ratio4, aes(x=time))
      ggp1<-ggp1+geom_line(aes(y=value ,group=interaction(molecule_type, cell, ORF,variable), color = cell, linetype=variable))
      if(showEB) ggp1<-ggp1+geom_ribbon(aes(ymin=lower, ymax=upper, group=interaction(molecule_type, cell, ORF,variable), color = cell, linetype=variable),  alpha=0.1)

       #ggp1<-ggp1+geom_line(aes(y=(logtotal-diff)/coeff ,group=interaction(molecule_type, cell, ORF), color = cell, linetype="dashed"))
       #ggp1<-ggp1+geom_point(aes(y=(logtotal-diff)/coeff ,group=interaction(molecule_type, cell, ORF), color = cell, shape=ORF,size=10))
       
       ggp1<-ggp1+ scale_y_continuous(
      name = "Log2 (total - sub-genomics)",
      sec.axis = sec_axis(~.*coeff+diff, name="Log2 TPM")
    
    )
        }
    
    ggp1<-ggp1+theme_bw()+theme(text = element_text(size=textsize), axis.text.x = element_text(size = rel(1.0),angle = 25, hjust=1.0))
    if(facet=="molecules"){
      ggp1<-ggp1+facet_grid(rows=vars(molecule_type))
    }else if(facet=="cells"){
      ggp1<-ggp1+facet_grid(rows=vars(cell))
    }else if(facet=="times"){
      ggp1<-ggp1+facet_grid(rows=vars(time))
    }
    else if(facet=="ORF"){
      ggp1<-ggp1+facet_grid(rows=vars(ORF))
    }
    return(list(ggp=ggp1, data=data))
  }else{
   
    ratio1$type=factor(ratio1$type,levels = levels)
    ORF="ORF"
    ord="Start"
    x1 =  paste("reorder(", ORF, ",", ord,")", sep="") 
    ratio1 = .expand(ratio1,nme="type")
    
    ggp<-ggplot()
    fill="type"
    if(facet=="molecules"){
      fill="cell_time"; #type
      ggp<-ggp+facet_grid(rows=vars(molecule_type))
    }else if(facet=="cells"){
      fill="mol_time"
      ggp<-ggp+facet_grid(rows=vars(cell))
    }else if(facet=="times"){
      fill="samp_cell"
      ggp<-ggp+facet_grid(rows=vars(time))
    }
    y_text = "Ratio"
    ggp<-ggp+geom_bar(data=ratio1, 
                      position=position_dodge(),aes_string(x=x1,y=y_text,fill=fill, colour=fill,ymin="lower" ,ymax="upper"),stat="identity")
  
   # ggp<-ggp+ geom_bar(position=position_dodge(), aes_string(y="Ratio"),stat="identity")
      #geom_bar(aes_string(x=x1, y="Ratio", fill = "type", colour = "type"),stat="identity", position = "dodge")
    if(showEB){
      ggp<-ggp+geom_errorbar(data=ratio1,position=position_dodge(width=0.9),
                             aes_string(x=x1,y=y_text,fill=fill, colour=fill,ymin="lower" ,ymax="upper"),colour="black")
    } #ggp<-ggp+geom_errorbar(aes_string(x=x1,ymin="lower", ymax="upper"), width=.2)#, position="dodge")
   ggp<-ggp+scale_y_continuous(limits = c(0,1))
    
     ggp<-ggp+ggtitle("Percentage of ORF covering reads which include leader")
     ggp<-ggp+xlab("ORF")
     ggp<-ggp+theme_bw()+theme(text = element_text(size=textsize), axis.text.x = element_text(size = rel(1.0),angle = 25, hjust=0.75))
    return(list(ggp=ggp, data =ratio1) )
  }
}


.calcConf<-function(v, conf.level=0.95, method="bayes"){
  #print(conf.level)
  if(conf.level<0.0001) return (c(NA, v[2]/v[1],NA))
  res = binom.confint(x=v[2],n=v[1],conf.level=conf.level, method=method)
  as.numeric( c(res$lower,res$mean,res$upper))
  
}
.overl<-function(x,y){
  min(x[2]-y[1], y[2] - x[1])
}
.makeCombinedArray<-function(clusters_, errors_, xlim, downsample=F, thresh = 1000, alpha = 1.0, ci = 0.995, max_num= 10, t=NULL,motifpos=list(), fisher=F){
  
  dm = dim(clusters_)
  ij =2
  if(dm[2]<=3){
    clusters_ = pivot_wider(clusters_, id_cols=pos, names_from="clusterID", values_from=names(clusters_)[3], values_fill=0)
    errors_= pivot_wider(errors_, id_cols=pos, names_from="clusterID", values_from=names(errors_)[3], values_fill=0)
    ij=1
  }
  range = which(clusters_[,ij] <= xlim[2] & clusters_[,ij]>=xlim[1])
  if(length(range)==0) return(ggplot())
  dn1 = list(c("depth","error"),unlist(clusters_[range,ij]),names(clusters_)[-(1:ij)])
  depth = array(dim = unlist(lapply(dn1,length)),  dimnames = dn1)
               
  depth[1,,] = as.matrix(clusters_[range,-(1:ij)])
  depth[2,,] = as.matrix(errors_[range,-(1:ij)])
 # print(dim(depth))
  .plotError(depth,  thresh = 1000, downsample = downsample, alpha = alpha, 
             ci = ci, max_num =max_num,t=t, xlim =xlim, pvAsSize=T, logy=T, fisher=fisher, motifpos=motifpos)
}
.restrict<-function(x,thresh){
 # print("h")
#  print(x)
  if(x[1]>thresh){
    ratio=thresh/x[1]
    x = floor(x*ratio)
  }
 # print("h1")
#  print(x)
  x
  
}
.plotError<-function(depth,range = 1:dim(depth)[[2]],alpha = 1.0, downsample = F, t1=NULL, method="logit",  
                     max_num = 20, pval_thresh = 1e-3,
                     ci=0.95, thresh = 1000, extend=T, log=F, adj=F, 
                     xlim = NULL,pvAsSize=T, logy=T, fisher=F, motifpos = list()){
  inds1 = apply(depth[1,range,],1,min)>thresh
  range = range[inds1]
  if(length(range)==0) return(ggplot())
 
  inds1 = apply(depth[1,range,],1,min)>thresh
  range = range[inds1]
  dfs = list()
  pos = as.numeric(dimnames(depth)[[2]][range])
  
  ##x dimension is 1  total vs error
  .testStatistic<-function(x) chisq.test(x)$p.value
  if(fisher){
  .testStatistic<-function(x) fisher.test(x)$p.value
  }
  d1 = depth[,range,]
  #print(d1[,1:5,])
  if(downsample){
  d1 =   apply(d1,c(2,3),.restrict,thresh)
  #print(d1[,1:5,])
  }
  pval = apply(d1,2, .testStatistic)
  if(adj) pval = p.adjust(pval, method="BH")
  if(max_num<length(pval)){
    ord  = order(pval)
    pval_thresh =  min(pval_thresh,pval[ord[max_num]])
  }
  for(i in 1:dim(depth)[[3]]){
    print(i)
    dfc = t(apply(depth[,range,i,drop=F],2,.calcConf, method=method,conf.level=ci))
    
    dimnames(dfc) = list(pos,c("lower","mean","upper"))
    type = rep(dimnames(depth)[[3]][i],length(pos))
    if(i==1){
      diffc = rep(0, length(pos))
    }else{
      dj = cbind(dfc[,2], dfs[[i-1]][,2])
      diffc = apply(dj,1,function(x) x[1]-x[2]) 
      
    }
    log10pv = -log10(pval)
    dfs[[i]] = data.frame(dfc,type,pos,pval,log10pv, diff=diffc)
    if(i==2){
      dfs[[1]]$diff = -1*dfs[[2]]$diff
    }
  }
  
  df = NULL
  
  for(i in 1:length(dfs)){
    df = if(is.null(df)) dfs[[i]] else rbind(df, dfs[[i]])
    
  }
  return(df)
}
.plotError1<-function(df,  t1, xlim, motifpos,pval_thresh=1e-3, logy=T, pvAsSize=T,ci=0.95, alpha=1.0){
  posy = max(df$upper)
 df$type=as.factor(df$type)
 #print(head(subset(df,pval <= pval_thresh & diff>0)))
  ty = levels(df$type)
  ggp<-ggplot(df, aes(x=pos,y=mean, colour=type,ymin=lower ,ymax=upper))
  if(pvAsSize){
  ggp<-ggp+ geom_point(position=position_dodge(), aes(y=mean, size=log10pv),alpha = alpha, stat="identity")
  }else{
    ggp<-ggp+ geom_point(position=position_dodge(), aes(y=mean),alpha = alpha,stat="identity")
    
  }
  if(ci>0.01){
 ggp<-ggp+geom_errorbar(position=position_dodge(width=0.0),colour="black")
  }
  ggp<-ggp+ggtitle(paste("Error rate by position (",ci*100,"% binomial CI) ",sep=""))
  ggp<-ggp+ylab("error rate")
  #if(log) ggp<-ggp+scale_y_continuous(trans='log10')
  ggp<-ggp+geom_text_repel(data=subset(df,pval <= pval_thresh & diff>0),
                  aes(pos,upper , label = sprintf("%3.2g" ,log10pv)),size = 3)
                    #if(type==ty[1]) "red" else "steelblue")
  ##can use pos instead of pval
  #if(!extend) ggp<-ggp+xlim(min(range), max(range))
  if(!is.null(t1) ){
    if( !is.null(t1$sideCols)){
    ggp<-ggp+geom_vline(xintercept = t1$Minimum, linetype="solid", color=t1$sideCols)
    ggp<-ggp+geom_vline(xintercept = t1$Maximum, linetype="dashed", color=t1$sideCols)
    }
  }
  if(!is.null(xlim)) ggp<-ggp+xlim(xlim)
  if(logy) ggp<-ggp+scale_y_continuous(trans='log10')
  if(length(motifpos)>0){
    for(jk in 1:length(motifpos)){
      ggp<-ggp+geom_vline(xintercept = motifpos[[jk]], linetype=jk+1, color="black", alpha=0.5)
    }
  }
  
  ggp+theme_bw()
}

.is_equal<-function(in1, in2){
    nmes = unique(c(names(in1),names(in2)))
   # print(in1)
  #  print(in2)
    for(nme_i in nmes){
     i1 = which(names(in1)==nme_i)
     i2 = which(names(in2)==nme_i)
   
     if(length(i1)==0 && length(i2)==0) return(TRUE)
     if(length(i1)!=1 || length(i2)!=1) return(FALSE)
     if(is.null(in1[[i1]]) && is.null(in2[[i2]])) return(TRUE)
     if(is.null(in1[[i1]]) || is.null(in2[[i2]])) return(FALSE)
     
      if(!identical(in1[[i1]],in2[[i2]])) return (FALSE)
    }
    return(TRUE)
  
}
writeFasta<-function(kmer10, file){
  write(">seq1", file)
  write(kmer10[1], file, append=T)
  for(i in 2:length(kmer10)){
    write(paste(">seq", i, sep=""), file, append=T)
    write(kmer10[i], file, append=T)
  }
}

getKmerPerc<-function(res) res$errorsK[[2]][2,]

matchKmerPerc<-function(kmerPerc_i, lev_i){
  prop = rep(0, length(kmerPerc_i))
  
  prop[match(lev_i[,1], names(kmerPerc_i))] = lev_i[,2]
  prop = as.numeric(prop)
  counts = prop
  prop = prop/sum(prop)
  names(prop) = names(kmerPerc_i)
  
  ratio = prop/kmerPerc_i
  names(ratio) = names(kmerPerc_i)
  #ratio
  cbind(counts, prop, kmerPerc_i, ratio)
}

getErrorComb<-function(error, sum, base_, bases_ = levels(as.factor(base_))){
  errorP = rep(0,length(bases_))
  sumP = rep(0,length(bases_))
  baseComp = rep(0, length(bases_))
  su  = sum(error)
  su1 = sum(sum)
  sul = length(error)
  for(k in 1:(length(bases_))){
    baseComp[k] = length(which(base_==bases_[k]))/sul
    errorP[k] =  sum(error[base_==bases_[k]])/su
    sumP[k] =  sum(sum[base_==bases_[k]])/su1
  }
  #errorP[5] =su
  #	baseComp[5] = sul
  errorComb = rbind(errorP, sumP, baseComp, errorP/sumP)
  
  errorComb = apply(errorComb,c(1,2), function(s) as.numeric(sprintf("%5.3g",s)))
  dimnames(errorComb) = list(c("error", "sum", 'comp', 'ratio') , as.character(bases_))
  errorComb
}

readTables<-function(files, files1){
  tables = list()
  for( i in 1:length(files)){
    tables[[i]] = read.table(files[i], header=T, sep=",", as.is=T)
  }
  names(tables) = files1
  tables
}

processMods<-function(tables,depth_thresh = 10, frac_thresh =0.4, depth_thresh1 = 10, skip = c(), kmer_v = seq(-1,1,1)){
  all_levs = list()
  files1 = names(tables)
  all_pos = matrix(nrow=0, ncol=2)
  nmes =  unlist(lapply(files1, split1))
  long_t = matrix(nrow = 0, ncol =6)
  bases = NULL
  
  errorsK = list()
  kmer = NULL
  kmer_bases = NULL
  kmer_ratios = NULL
  base_ratios = NULL
  kmer_abs = NULL
  for( i in 1:length(tables)){
    a = tables[[i]]
    errorP = rep(0,length(bases))
    baseComp = rep(0, length(bases))
    
    sum = a$match + a$mismatch  + a$baseDel#+a$baseIns #+ a$refClipped
    if(is.null(kmer)) {
      
      kmer = getKmer(a$base, 1:length(a$base), kmer_v)
      kmer_bases = levels(as.factor(kmer))
      bases = levels(as.factor(a$base))
      base_ratios = matrix(nrow = length(bases), ncol = length(files))
      
      kmer_ratios = matrix(nrow = length(kmer_bases), ncol = length(files))
      kmer_abs = matrix(nrow = length(kmer_bases), ncol = length(files))
      dimnames(kmer_ratios) = list(as.character(kmer_bases), nmes)
      dimnames(kmer_abs) = list(as.character(kmer_bases), nmes)
      dimnames(base_ratios) = list(as.character(bases), nmes)
    }
    error =  a$mismatch  + a$baseDel#+a$baseIns #+ a$refClippedR
    frac = (error /sum)
    
    
    
    
    
    errorsK[[i]] = getErrorComb(error,sum,  kmer, kmer_bases)
    
    kmer_ratios[,i] = errorsK[[i]][4,]
    kmer_abs[,i] = errorsK[[i]][1,]
    frac[sum < depth_thresh1] = NA
    a1= cbind(a[,1:4], sum, frac)
    type =  rep(nmes[i], length(sum))
    long_t = rbind(long_t, cbind(a1, type))
    #all[[i]] = a1
    all_levs[[i]] = getlev((a1[which(a1$frac > frac_thresh & a1$sum>depth_thresh),,drop=F]$base))
    #inds = which(a1$sum>depth_thresh)
    inds1 = which(a1$frac > frac_thresh & a1$sum>depth_thresh)
    if(i %in% skip){
      print(paste("skipping", i))
    }else{
      all_pos = rbind(all_pos, cbind(inds1, a1$sum[inds1]))
    }
    #all_mod[[i]] = a1[inds1,,drop=F]
    if(i==1){
      
      big_t = a1[,c(1,2,5,6)]
      names(big_t)[3:4] = paste(names(big_t)[3:4], nmes[i], sep=".")
    }else{
      big_t1 = a1[,c(5,6)]
      names(big_t1) = paste(names(big_t1), nmes[i], sep=".")
      big_t = cbind(big_t, big_t1)	
      
    }	#plot(a1$pos[inds], a1$frac[inds])
  }
  names(all_levs) = nmes
  all_pos1 = unique(sort(all_pos[,1]))
  depth1 = rep(NA, length(all_pos1))
  for(i in 1:length(all_pos1)) depth1[i] = max(all_pos[all_pos[,1] == all_pos1[i],2], na.rm=T)
  big_t1 = cbind(big_t[,grep('frac', names(big_t), invert=T)], big_t[,grep('frac', names(big_t))])
  big_t1 = cbind(big_t1[,c(1,2)], apply(big_t1[,-c(1,2)], c(1,2), function(s) sprintf("%5.3g", s)))
  
  names(errorsK) = nmes
  or_r = order(kmer_ratios[,1], decreasing = T)
  or_abs = order(kmer_abs[,1], decreasing = T)
  head_ratios =kmer_ratios[or_r,]
  head_abs = kmer_abs[or_abs,]
  list(big_t = big_t1, long_t = long_t, all_pos = cbind(all_pos1, depth1), all_levs = all_levs,  errorsK = errorsK, kmer_ratios = kmer_ratios, kmer_abs = kmer_abs, kmer = kmer, kmer_bases = kmer_bases, or_r = or_r, or_abs = or_abs, head_ratios = head_ratios, head_abs = head_abs)
}

printRes<-function(results, v = 1:length(results), top = 5, offsets = 0:10, ks = 1:10, bases = c("A", "C", "T", "G"), remove_repeat = T){
  out_res = list()
  nmes = c()
  k = 1;
  for(i in v){
    if(nchar(results[[i]]$kmer[1]) %in% ks){
      offset = which(as.numeric(strsplit(names(results)[i],"_")[[1]])==0)-1
      if(offset %in% offsets){
        nmes[k] = names(results)[i]
        print(nmes[k])
        ratios = results[[i]]$head_ratios;
        
        nonblank = !unlist(lapply(dimnames(ratios)[[1]], function(x) length(grep('-', x))))
        rownmes = data.frame(strsplit(dimnames(ratios)[[1]], ""))
        ba = as.character(unlist(rownmes[offset+1,]))
        prior =if(offset<=0) rep("N", length(ba)) else as.character(unlist(rownmes[offset,])) 
        subs =if(offset+2>dim(rownmes)[1]) rep("N", length(ba)) else as.character(unlist(rownmes[offset+2,])) 
        indsr = if(remove_repeat)  which(nonblank & ba!=prior & ba != subs & ba %in% bases) else which(nonblank & ba %in% bases)
        out_res[[k]] = ratios[indsr,,drop=F]
        print("top")		
        print(head(out_res[[k]], top))
        print("bottom")	
        print(tail(out_res[[k]], top))
        k = k+1
        
        
      }
    }
  }
  names(out_res) = nmes
  invisible(out_res)
}

motifAnalysis<-function(results, ks, offset, base, top=1, remove_repeat=T){
  resu = printRes(results, top = 4, ks = ks, offsets = offset, bases = base, remove_repeat=remove_repeat)
  dn = getDimNames(resu)
  subsetTop = findKmers(results, kmer=head(dn[[1]],top), offsets = offset)
  subsetBottom = findKmers(results, kmer=tail(dn[[1]],top), offsets = offset)
  
  list(dimTop = lapply(subsetTop, dim), top = lapply(subsetTop, getAvg),dimBottom = lapply(subsetBottom, dim),  bottom = lapply(subsetBottom, getAvg))
  
}

.expandLev<-function(lev_){
	m1 = apply(t(data.frame(lapply(lev_[,1], function(x) strsplit(as.character(x), "\\.")[[1]]))), c(1,2), as.numeric)
	mat  = cbind(m1,lev_[,2])
	dimnames(mat)= NULL
	mat
}
.stEndKmer<-function(lev_a,right= 0:8, left = -8:-1){
	res = data.frame(
	breakSt = lev_a[,1],
	breakStart_l = getKmer(fasta[[1]], lev_a[,1],right),
	breakStart_r = getKmer(fasta[[1]], lev_a[,1],left),
	breakEnd = lev_a[,2],
	breakEnd_l = getKmer(fasta[[1]], lev_a[,2],right),
	breakEnd_r = getKmer(fasta[[1]], lev_a[,2],left),
	cnt = lev_a[,3]
	)
	res
}


getDimNames<-function(out_res, type = 1){
  l = list()
  for(i in 1:length(out_res)){
    l[[i]] = dimnames(out_res[[i]])[[type]]	
  }
  names(l) = names(out_res)
  l
}

getKmerOffsets<-function(max = 5){
  kmers = list()
  kmers[[1]] = 0
  k = 2
  for(i in 1:max){
    km_i = 0:i
    for(j  in 1:length(km_i)) {
      kmers[[k]] = km_i-j+1
      print(kmers[[k]])
      k =k+1;
    } 
  }
  
  names(kmers) = unlist(lapply(kmers, paste, collapse="_"))
  kmers
}

getAvg<-function(matr){
  apply(apply(matr[,grep('frac', names(matr))], c(1,2), as.numeric),2,mean, na.rm=T)
}

findKmers<-function(results, kmers="CTC", offsets = -5:5, col_inds = c(1,2,12,13,14,15,16,17,18,20)){
  
  out_res = list()
  k = 1
  nmes = c()
  inds_res = list()
  for(kmer in kmers){
    for(i in 1:length(results)){
      offset = which(as.numeric(strsplit(names(results)[i],"_")[[1]])==0)-1
      if(offset %in% offsets){
        res = results[[i]]
        if(kmer %in% res$kmer){
          nmes = c(nmes,paste(names(results)[i], kmer,sep="_"))
          inds = which(res$kmer==kmer)
          inds_res[[k]] = inds
          out_res[[k]] = res$big_t[inds,col_inds]
          k = k+1
        } 
      }
    }
  }
  names(inds_res) = nmes
  names(out_res) = nmes
  out_res
  #list(out_res = out_res, inds_res = inds_res)
}

filterT<-function(big_t, range){
  big_t[big_t$pos <= range[2] & big_t$pos>=range[1],,drop=F]
}

plotF<-function(long_t, base="C"){
  
  lt1 = long_t[which(long_t$base==base),]
  ggp<-ggplot(lt1,aes(x=pos,fill = type, colour = type, y=frac))+geom_point() + theme_bw()+ggtitle(paste(base, " %ESB"))
  #ggp<-ggp+scale_y_continuous(limits = c(0, 1))
  invisible(ggp)
}



paste1<-function(st, end) paste(st, end, sep="/")
plotFAll<-function(lt1, range = c(20000,30000), exclude = "allreads"){
  lt1$frac = 100*lt1$frac
  if(!is.null(exclude)) {
    
    lt1 = lt1[ grep(exclude, lt1$type, inv=T),]
  }
  long_t = lt1[lt1$pos < range[2] & lt1$pos>=range[1],]
  
  pos  = long_t[!is.na(long_t$frac),]$pos
  min = min(pos)
  max = max(pos)
  
  #long_t = long_t[long_t$pos <= max(pos) & long_t$pos>=min(pos),, drop=F]
  #miny = min(long_t$frac, na.rm=T)
  maxy = max(long_t$frac, na.rm=T)
  ggps = list()
  charc = c("A", "C", "G", "T")
  for(i in 1:length(charc)){
    ggps[[i]] = plotF(long_t, charc[i])+scale_y_continuous(limits = c(0, maxy))+scale_x_continuous(limits = c(min, max))
    
  }
  names(ggps) = charc
  
  ml<-marrangeGrob(ggps,nrow = 2, ncol = 2) 
  invisible(ml)
}


getlev<-function(x1, todo = NULL){ 
  x = if(is.null(dim(x1))) x1 else apply(x1, 1,paste, collapse=".")
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


getSideCols<-function(rn, t, default = "#000000"){
  rown = as.numeric(rn)
  sideC =  as.character(t$sideCols)
  cols = rep(NA, length(rown))
  #t$Maximum = floor(t$Maximum/10)*10+1
  #t$Minimum = floor(t$Minimum/10)*10+1
  for(i in 1:length(rown)){
    mtch = which(t$Minimum<=rown[i] & t$Maximum>=rown[i])
    
    if(length(mtch)==0) cols[i] = default
    else cols[i] = sideC[mtch[1]]
  }
  cols
}


plot1<-function(x,y, t, log = ""){
  plot(x,y, log=log)
  abline(v = t$Minimum, col=2)
  abline(v = t$Maximum, col=3)
  
}
split1<-function(fi) strsplit(fi,"_")[[1]][1]




plotClusters<-function(df,seq_df, k1, totalReadCount, t, motifpos, peptides,rawdepth = T, 
                       linetype="sampID", colour="clusterID", title = "", ylab=if(rawdepth)  "depth" else "TPM", logy=F, 
                       leg_size = 6, xlim  = NULL, show=F, updatenmes = F, fill = F, alpha=0.5,
                       showSeqThresh=500, size=20, linesize=0.1, textsize=20,zoom=T){
  if(!is.factor(df$clusterID)) df$clusterID = as.factor(df$clusterID)  #types[df$type]
 # names(df)[3] = "depth"
  #ids =  as.character(rel_count$ID)
  #if(max(rel_count[,2])>1) stop('this should not happen')
  
  if(!rawdepth){
    
    #for(i in 1:length(ids)){
    #  df[df$type== ids[i],2] = (df[df$type== ids[i],2] * rel_count[i,2])
    #}
    df[,k1]  = df[,k1]*(1e6/totalReadCount)
  }
 # if(logy){
#    df[df[,k1]<=0.01,k1]=0.01
#  }
  if(is.null(xlim)){
    xlim = c(min(df$pos), max(df$pos))

    ylim = c(max(0,min(df[,k1], na.rm=T)),max(1,max(df[,k1], na.rm=T)))
  }else{
    xincl = df$pos <= xlim[2] & df$pos >=xlim[1]
    ylim = c(max(0,min(df[xincl,k1], na.rm=T)),max(1,max(df[xincl,k1], na.rm=T)))
  }
  ggp<-ggplot()
 # print(paste("linesize",linesize))
  ggp<-ggp+geom_line(data=df, aes_string(x="pos", fill="clusterID", colour = colour, linetype=linetype, y = names(df)[k1], alpha=alpha)) +theme_bw()
  if(!is.null(seq_df) && dim(seq_df)[1]<showSeqThresh){
   seqy = seq_df$seqy
    seq_df$seqy = rep(ylim[2],length(seqy))
  #  seq_df$seqy = seqy 
    print(seq_df)
    ggp<-ggp+geom_text(data=seq_df,aes(x=pos, y=seqy, color=sequence, label=sequence))
  }
   # ggp<-ggp+geom_line(inherit.aes=T, aes(alpha=alpha)) 
if(fill) ggp<-ggp+geom_area()
  ggp<-ggp+ggtitle(title)
  if(logy) ggp<-ggp+scale_y_continuous(name=ylab,trans='log10') #, limits=ylim)
  else ggp<-ggp+labs(y= ylab)+ylim(ylim)
  if(leg_size==0  || length(levels(as.factor(as.character(df$type))))>20){
    ggp<-ggp+theme(legend.position="none",text = element_text(size=textsize))
  }else{
    ggp<-ggp+theme(text = element_text(size=textsize),plot.title = element_text(size = rel(1.0), face = "bold"),legend.position="right",
legend.title=element_text(size=leg_size), legend.text=element_text(size=leg_size))
  }  
  #ggp<-ggp+theme(, axis.text.x = element_text(size = rel(1.0)))
  
  #
  if(updatenmes){
    
    br = apply(clusters[range,,drop=F], 1, findBreaks)
    br_names = apply(br, 2, function(x) paste(x[1:2], collapse=":"))
    
    ggp<-ggp+ scale_colour_discrete(name="cluster",
                                    breaks = levels(df$type),
                                    labels=br_names[as.numeric(levels(df$type))])
  }
  
  if(!is.null(t$sideCols)){
    ggp<-ggp+geom_vline(xintercept = t$Minimum, linetype="solid", color=t$sideCols, alpha=0.5)
    ggp<-ggp+geom_vline(xintercept = t$Maximum, linetype="dashed", color=t$sideCols, alpha=0.5)
  
  }
  if(length(motifpos)>0){
    for(jk in 1:length(motifpos)){
    ggp<-ggp+geom_vline(xintercept = motifpos[[jk]], linetype=jk+1, color="black", alpha=0.5)
    }
  }

if(!is.null(peptides)){
  ggp<-ggp+geom_vline(xintercept = peptides[,1], linetype="dotted", color="blue", alpha=0.5)
  ggp<-ggp+geom_vline(xintercept = peptides[,2], linetype="dashed", color="blue", alpha=0.5)
  #ggp<-ggp+geom_vline(xintercept = fimo$start[fimo$strand=="-"], linetype="dotted", color="grey")
}
  #abline(v = t$Maximum, col=3)
  if(!is.null(xlim)) {
    if(zoom) ggp<-ggp+facet_zoom(xlim=xlim) else ggp<-ggp+xlim(xlim)
  }
 # ggp<-ggp
  if(show) print(ggp)
  invisible(ggp)
}

#ma <- function(x, n = 5){
#	filter(x, rep(1 / n, n), sides = 2)

#}
distbin<-function(x) dist(x, method="binary")

.loess<-function(clusters, span){
  cars.lo <- stats::loess(depth ~ pos, clusters, span = span)
        y = try(stats::predict(cars.lo, clusters, se = FALSE))
 y
}

loess_smooth<-function(clusters,  inds,span = 0.05){
  if(span==0) return (clusters)
  nmes = names(clusters)
#print(nmes)
#print(dim(clusters))
    for(j in inds){	
      names(clusters)[j] = "depth"
     
      if(max(clusters$depth,na.rm=T)>0){
        y = try(.loess(clusters, span))
	if(inherits(y,"try-error")) {
		print("error loess")
		#print(which(is.na(clusters$pos)))
		#plot(clusters$pos,clusters$depth)
		#print()
        }else{
		#print("is ok ")
        	y[clusters$depth==0] =0
       	 	clusters[,j]=  y
	}
      } 
      names(clusters)[j] = nmes[j]
    }
  
  clusters
}




makeCovPlots<-function(clusters_,  header, total_reads, t, type_nme,  leg_size=6, logy=F, rawdepth = T, xlim = list(NULL),   show=F, fill=F
){
 dinds  = grep("depth", header)
  ggps = list()
  for(j in 1:length(xlim)){
	leg_size1 = if(j==1) leg_size else  0
    for( k in 1:length(type_nme)){
      ggps[[length(ggps)+1]] = plotClusters(clusters_, dinds[k],  total_reads[k], t, fimo, xlim = xlim[[j]], rawdepth = rawdepth,  
                                            title =type_nme[k], logy=logy, leg_size =leg_size1, show=show, fill =fill)
    }
  }
  
  ml<-marrangeGrob(ggps,nrow = length(type_nme), ncol=length(xlim)) 
  
  invisible(ml)
}
#prop0 how far back to go 
getGeneBP<-function(t, genes, left, right, left_buffer = 10){
  endcs = matrix(ncol = 2, nrow = length(genes))
  dimnames(endcs) = list(genes, c("start","end"))
  for(i in 1:length(genes)){
    ind  = which(t$gene==genes[i])
   # t[(ind-1):(ind+1),]
    endcs[i,1] = max(t$Minimum[ind]-left, t$Minimum[ind-1]+left_buffer)
    
   #   t$Minimum[(ind-1)]+prop0 * (t$Maximum[(ind-1)] - t$Minimum[(ind-1)])
    #endcs[i,1] = max(t$Minimum[ind] - maxl, min( t$Minimum[ind] - 100, endcs[i,1]))
    endcs[i,2] = min(t$Minimum[ind] + right, t$Minimum[ind+1])
      #t$Minimum[ind]+prop1 * (t$Maximum[(ind)] - t$Minimum[(ind)])
    
  }
  endcs
}
# brP1 = readBreakPointsH5(h5file,"chrMT007544", "", 0)
readBreakPointsH5<-function(h5file, chrom, type,j, prefix=""){
  id = paste(prefix,chrom,type,j, sep="/")
  mat=h5read(h5file,id)
  heatm = mat[-(1:2),-(1:2)]
  rowlen = nrow(heatm)
  collen= ncol(heatm)
 rows =cbind( data.frame(mat[-(1:2),1:2]), rep("start", rowlen), rep(type,rowlen ))
  cols = cbind(data.frame(t(mat[1:2,-(1:2)])), rep("end", collen), rep(type,collen ))
 
  
  names(rows) = c("pos", "depth", "s_e", "type")
  names(cols) = c("pos", "depth", "s_e", "type")
  dimnames(heatm) = list(rows[,1], cols[,1])
  res =list(heatm = t(as.matrix(heatm)), rows = cols, cols = rows)
  #checkHeatm(res)
  res
}


readBreakPoints<-function(f, type, addOne = F){
	if (length(scan(f,n=2, what="raw"))==0) return (NULL)
  aa = try(read.table(f, sep=",", skip=2))
  nme  = read.table(f, sep=",",nrows = 2 )
  if(dim(nme)[2]<dim(aa)[2]) nme = cbind(c(NA, NA), nme)
  names(nme) = names(aa)
  rowlen = dim(aa)[1]
  collen = dim(nme)[2]-2
  rows  = data.frame(cbind(aa[,1:2], rep("start", rowlen), rep(type,rowlen )))
  cols = data.frame(t(nme[,-(1:2)]), rep("end", collen), rep(type,collen ))
  
  heatm = aa[,-(1:2)]
  
  
  names(rows) = c("pos", "depth", "s_e", "type")
  names(cols) = c("pos", "depth", "s_e", "type")
  
  if(addOne){
    pos_ind = 1  #which(names(rows)=="pos")
   # print(head(rows))
    rows[,pos_ind] = rows[,pos_ind]+1
    cols[,pos_ind] = cols[,pos_ind]+1
  #  print("after")
   # print( head(rows))
  }
  dimnames(heatm) = list(rows[,1], cols[,1])
  
  
  res =list(heatm = as.matrix(heatm), rows = rows, cols = cols)
  #checkHeatm(res)
  res
}
checkHeatm<-function(breakP){
  dimh = dim(heatm)
  if(dimh[1]>1 && dimh[2]>1){
    rowsums = apply(breakP$heatm,1,sum)
    colsums = apply(breakP$heatm,2,sum)	
    if(sum(abs(breakP$rows$depth - rowsums))>1) stop('row error')
    if(sum(abs(breakP$cols$depth - colsums))>1) stop('col error')
  }
  
}

findOrf<-function(fasta){
  seqlen  = length(fasta[[1]])
  mat = as.vector(gregexpr("atg", paste(fasta[[1]],collapse=""))[[1]])
  prots = list()
  lens = c()
  for(i in 1:length(mat)){
    start = mat[i]
    prot_a = translate(fasta[[1]][start:seqlen])
    end_p = grep("\\*",prot_a)[1]-1
    
    end = start + end_p*3-1 + 3
    full_prot = ""
    
    if(!is.na(end) & !is.na(start)){
    if(end-start>3)full_prot = translate(fasta[[1]][start:end])
    }
    prots[[i]] = paste(full_prot,collapse="")
    lens = c(lens, length(full_prot)-1)
    
  }
  names(prots) = mat
  names(lens) = mat
  list(prots = prots, length = lens)
}

subsetBr<-function(breakP, region){
  startc = region[1:3]
  endc = region[4:6]
  rowinds = which(breakP$rows$pos>=startc[1] & breakP$rows$pos<=startc[2])
  colinds = which(breakP$cols$pos>=endc[1] & breakP$cols$pos<=endc[2])
  rows = breakP$rows[rowinds,,drop=F]
  cols = breakP$cols[colinds,,drop=F]
  heatm = breakP$heatm[rowinds,colinds,drop=F]
  rows$depth = apply(heatm,1,sum)
  cols$depth = apply(heatm,2,sum)
  colinds = cols$depth>0
  rowinds = rows$depth>0
  res = list(heatm = heatm[rowinds,colinds,drop=F],rows = rows[rowinds,,drop=F] ,cols = cols[colinds,,drop=F] )
  #checkHeatm(res)
  res
}
.findInd<-function(v,br) min(which(br>=v[1]))

.getMoments<-function(startP){
  maxv = startP[which(startP[,2]==max(startP[,2],na.rm=T))[1],]
  a=apply(startP,1,function(v) v[1]*v[2])
  avg = sum(a)/sum(startP[,2])
  
  a=apply(startP,1,function(v) ((v[1]-avg)^2)*v[2])
  var = sum(a)/sum(startP[,2])
  stderr = sqrt(var)
  c(avg,stderr,maxv[[1]], maxv[[2]])
}




#tolerance is how many base pairs past the break we allow
.protR<-function(start, fasta,t1,protL,
                 frames = 0:2, tolerance = 5, seqlen=length(fasta[[1]]) ){
 # print(paste(start,frame))
 # frame = protL$phase
  if(start>seqlen) stop("start cant be more than seqlen")
  t1_ind =  which(t1$Minimum>=start- tolerance)[1]
  #if this is NA then there no downstream gene
  left_rna = protL$rna
  if(!is.na(protL$start) ){
    phase=(protL$end-protL$start+1)%%3
    frame =(3-phase)%%3
    t1_ind_a =  which(t1$Maximum>=start & t1$Minimum<=start)
    if(length(t1_ind_a)>0) t1_ind  = t1_ind_a 
  }else{
    frame =0
   
      if(!is.na(t1_ind)){
       start = t1$Minimum[t1_ind]
       frames = 0
      }
    
  }
  prots = list()
  rnas = list()
  for(kk in 1:length(frames)){
    frame1 = (frames[kk]+frame)

    rna = c(left_rna,rep("N",frames[kk]), fasta[[1]][start:seqlen])
    full_prot = "NA"
    prot_a = try(translate(rna))
     if(inherits(prot_a,"try-error")) {
	print("translation error")
	print("left")
	print(left_rna)
	print("right")
	
	print(fasta[[1]][start:seqlen])

     }else{
   	 end_p = grep("\\*",prot_a)[1]-1
         end = start + end_p*3-1 + 3 #   ##note plus 3 includes the stop
   	 full_prot = ""
    	if(!is.na(end)){
    	   if(end-start>3){
		rna1 = rna[1:(end-start+1)]
		full_prot = try(translate(rna1))
		 if(inherits(full_prot,"try-error")) {
			print("translation error1")
			print(rna1);
		}
	   }
        }
     }
    #prot_a = translate(fasta[[1]][start:seqlen],frame = frame1%%3)
   
    t1_ind1 =  which(t1$Maximum>=end & t1$Minimum<=end)
    gene = as.character(t1$gene[c(t1_ind, t1_ind1)])
    #rna = c(left_rna,rep("N",frames[kk]), fasta[[1]][start:end])
  
    
     stops = which(full_prot=="*")/length(full_prot)
    if(length(stops)>0 && stops[1] < 0.9999){
      print(c(protL$end, protR$start))
      print(kk)
      print(stops)
      stop(" error")
    }
    
    rnas[[kk]] = rna
    prots[[kk]] = full_prot
    #translate(c(leftRNA, rna))
  }
    
  
  list(prots = prots, rnas = rnas, frame = frame,gene =unique(c(protL$gene, gene)), rna_br = length(protL$rna), prot_br = length(protL$prot))
}



.combinedProt<-function(m1, fasta,t1,seqlen=length(fasta[[1]])){
  col_offset = dim(m1)[2]
  nrow = dim(m1)[1]
 
  ncol1 = col_offset+5
  nmes = c("prot0","prot1","prot2","genes","prot_br" , "rna-br","phaseL","numgenes")
  stringcols = 1:4
  res = matrix("" ,nrow = nrow, ncol = length(nmes))
  
  #names(res)[1:col_offset] =dimnames(m1)[[2]]
  dimnames(res)[[2]] = nmes
  if(nrow>0){
    for(i in 1:nrow){
      
      protL = .protL(m1[i,1] ,fasta,t1)
     # if(is.na(m1[i,2])) stop("err")
      protJoin = .protR(m1[i,2],fasta,t1,protL, seqlen=seqlen )
    #  protJoin = .joinProts(protL, protR)
      
      
      res[i, ] =
        c("","","",
          paste(protJoin$gene, collapse="_")[[1]],
          protJoin$prot_br,
          protJoin$rna_br,
          protJoin$frame,
          length(protJoin$gene))
      for(jk in 1:length(protJoin$prots)){
        res[i,jk] =paste(protJoin$prots[[jk]], collapse="")[[1]]
      }
    }
  }
  res2 = cbind(res,m1)
  res2= apply(res2,c(1,2), as.character)
  #print(res2)
  res2
}
#end is end of first section 
#phase is -1 if there is no ORF here
.protL<-function(end, fasta, t1){
  t1_ind = which(t1$Minimum<end)[1]
  
  if(is.na(t1_ind)){
    rna = c()
    phase=-1
    prot=c()
    gene = c()
    st = NA
  }else{
    st = t1$Minimum[t1_ind]
    rna = fasta[[1]][st:end]
    gene = unique(as.character(t1$gene[t1_ind]))
    if(is.na(gene)) stop("err")
  }
  list(rna = rna, gene = gene, end = end, start = st)
}



.getDepth<-function(breakP_,m1){
  dep = rep(0, dim(m1)[1])
  if(dim(m1)[1]==0) return (dep)
  for(i in 1:length(dep)){
  rowind = which(breakP_$rows$pos==m1[i,1])
  colind = which(breakP_$cols$pos==m1[i,2])
  if(length(rowind==1) & length(colind)==1)  dep[i] = breakP_$heatm[rowind,colind]
  }
  dep
}

.getMaxInds<-function(breakP_, mindepth, todo, ind,dnmes){
  heatm = breakP_$heatm
  dm = dim(heatm)
  ord = order(heatm, decreasing=T)
  vals = heatm[ord]
  max_i = which(vals>=mindepth)
  ord = ord[max_i]
 # dnmes = c("start", "end", paste("depth", todo,sep=""))
  m1 = array(dim = c(length(ord),length(dnmes)), dimnames = list(NULL,dnmes ))
  
  k1 = 1
  k = 1
  while(k1 <= length(ord)){
    matr_inds = which(heatm == heatm[ord[k]], arr.ind = TRUE)
   # print(matr_inds)
    max_k2 = min(dim(matr_inds)[1], length(ord)-k1)
    for(k2 in 1:max_k2){
      if(k2>dim(matr_inds)[1]) stop(" err" )
      if(k1<=length(ord)){
      #  print(paste("k1" ,k1))
        row_k1 = c(breakP_$rows[matr_inds[k2,1],1], breakP_$cols[matr_inds[k2,2],1],ind)
       # if(is.na(row_k1)[2]) {
      #    print((matr_inds))
       #   print(k2)
        #  stop(" err1" )
        #}
       # print(row_k1)
        m1[k1,1:3] = row_k1
        m1[k1,ind+3]= heatm[ord[k1]]
        k1 = k1 +1
      }
    }
    k= k+1
  }
  #if(length(which(is.na(m1[,2])))>0) stop(" err")
  m1
}

.mergeProts<-function(m1s_j, nme="prot0", ind = which(dimnames(m1s_j)[[2]] %in% nme)){
  if(length(ind)<1){
    print(dimnames(m1s_j)[[2]])
    print(nme)
    stop("err")
  }
  joins =  apply(m1s_j[,ind,drop=F],1, function(x) paste(as.character(x), collapse="."))
  lev = unique(joins)
# print(head(joins));
  #print(lev)
  if(dim(m1s_j)[[1]]==0) return (NULL)
#  print(lev)
  m1_merged = m1s_j[1:length(lev),,drop=F]
  depth_inds = grep("depth" , dimnames(m1s_j)[[2]])
  for(kk in 1:length(lev)){
    indskk = which(joins==lev[kk])  
    m1_merged[kk,] = m1s_j[indskk[1],]
    m1_merged[kk,depth_inds] =apply( m1s_j[indskk,depth_inds,drop=F],2,function(x) sum(as.numeric(x)))
  }
  m1_merged
}

getHMClust<-function(breakPs, regions, fasta,t,mind=c(1,1), seqlen=length(fasta[[1]])){
  if(is.null(names(breakPs))) stop("error with breakPs names")
  if(is.null(names(regions))) stop("error with regions names")
  
  t1 = t[t$gene!="none" & t$gene!="leader",]
  dimm = c(length(regions),length(breakPs),  2,2)
  res_df = matrix(nrow =prod(dimm[1:3]) ,ncol = length(dimm)+3 )
  df = data.frame(res_df)
  k=1
  m1s_ = list();
 
  todo = 1:length(breakPs)
  dnmes = c("start", "end", "index", paste("depth", todo,sep=""))
#  prot0 prot1 prot2 genes prot_br rna.br phaseL numgenes start  end index
  #dnmes1 = c("prot0","prot1","prot2","genes","prot_br","rna_br","phaseL", "numgenes",dnmes);
  m1s_merged = NULL;
  m1s_merged1 = NULL;# array(dim = c(0,length(dnmes1)), dimnames = list(NULL,dnmes1 ))
  m1s_merged2 = NULL
   for(j in 1:length(regions)){
    if(TRUE){
    region = regions[[j]]
    m1 = array(dim = c(0,length(dnmes)), dimnames = list(NULL,dnmes ))
    
    breakPs_ = list()
    for(i in todo){
      breakPs_[[i]] = subsetBr(breakPs[[i]],region)
    }
    for(i in todo){
      breakP_ = breakPs_[[i]]
      heatm = breakP_$heatm
      m2 = .getMaxInds(breakP_, mind[i], todo,i,dnmes)
      for(j1 in todo[-i]){
        m2[,3+j1]   =.getDepth(breakPs_[[j1]], m2)
      }
      m1 = rbind(m1,m2)
      
      startP = breakP_$rows[breakP_$rows$s_e=="start",1:2,drop=F]
      endP = breakP_$cols[breakP_$cols$s_e=="end",1:2,drop=F]
      momentsS = .getMoments(startP)
      momentsE = .getMoments(endP)
      df[k,] =c(j,i,1,momentsS)
      df[k+1,] =c(j,i,2, momentsE )
      k =k+2;
    }
    m1s_j = .combinedProt(m1, fasta,t1, seqlen = seqlen)
      m1s_[[j]] =m1s_j
      if(dim(m1s_j)[1]>0){
      m1s_merged_j = .mergeProts(m1s_j,nme="prot0")
      m1s_merged_j1 = .mergeProts(m1s_j,nme="genes")
      m1s_merged_j2 = .mergeProts(m1s_j,nme=c("start","end"))
       # print(head(m1s_merged_j))
        if(is.null(m1s_merged)) m1s_merged = m1s_merged_j
        else  m1s_merged =rbind(m1s_merged,m1s_merged_j)
      if(is.null(m1s_merged1)) m1s_merged1 = m1s_merged_j1
      else  m1s_merged1 =rbind(m1s_merged1,m1s_merged_j1)
      if(is.null(m1s_merged2)) m1s_merged2 = m1s_merged_j2
      else  m1s_merged2 =rbind(m1s_merged2,m1s_merged_j2)
      }
    }
  }
    names(df) = c("region" ,"type" , "s_e", "mean", "sd","maxpos", "maxdepth")
  names(m1s_) = names(regions)
 # names(m1s_merged) = names(regions)
 for(i in 1:length(breakPs)) df[which(df$type==i),2] = names(breakPs)[i]
  for(i in 1:length(regions)) df[which(df$region==i),1] = names(regions)[i]
  for(i in 1:2) df[which(df$s_e==i),3] =c("start","end")[i]
  df$type = factor(df$type, levels = names(breakPs))
  df$region  = factor(df$region, levels = names(regions))
  df$s_e  = factor(df$s_e, levels = c("start","end"))
  m1s_merged = data.frame(m1s_merged)
  m1s_merged[,-(1:4)] = apply(m1s_merged[,-(1:4)],c(1,2),as.numeric)
  m1s_merged1 = data.frame(m1s_merged1)
  m1s_merged2 = data.frame(m1s_merged2)
  m1s_merged1[,-(1:4)] = apply(m1s_merged1[,-(1:4)],c(1,2),as.numeric)
  m1s_merged2[,-(1:4)] = apply(m1s_merged2[,-(1:4)],c(1,2),as.numeric)
  list(df = df,breaks = m1s_, merged = m1s_merged,merged_gene = m1s_merged1,merged_se = m1s_merged2)
 # invisible(ggp)
}

plotHMClust1<-function(hmClust_b,total_reads,type_name, size=3, logT = T, nudge_x = 0.1, nudge_y = 0.1,  plotDepth = FALSE, extra = 0.1){
#shifty = 1.2
#  shiftx = 1.5
  merged = hmClust_b$merged_gene
  depth_inds = grep("depth", names(merged))
  for(k in 1:length(depth_inds)){
    if(!plotDepth) merged[,depth_inds[k]] = merged[,depth_inds[k]]*(1e6/total_reads[k])+extra
    else merged[,depth_inds[k]] = merged[,depth_inds[k]]+extra
  }
 label="TPM + "
 if(plotDepth) label = "depth +"
# merged$depth2 = merged$depth2+0.1
 ggps = list()
 for(k in 2:length(depth_inds)){
   depth1 = names(merged)[depth_inds[k]]
   for(k1 in 1:(k-1)){
     depth2 = names(merged)[depth_inds[k1]]
     print(paste(depth1, depth2))
     nme_old =  names(merged)[depth_inds[c(k,k1)]]
     names(merged)[depth_inds[c(k,k1)]] = c("deptha","depthb")
  ggp<-ggplot(merged, aes(x=deptha, y=depthb))+geom_point()
 # ggp<-ggp+theme_bw()#plot.title = element_text(size = 8, face = "bold"),legend.title=element_text(size=8), legend.text=element_text(size=8))
 ggp<-ggp+	theme(plot.title = element_text(size = 8, face = "bold"),legend.title=element_text(size=8), legend.text=element_text(size=8))

   ggp<-ggp+geom_text(  size=size,nudge_x = nudge_x,nudge_y= nudge_y, aes(x=deptha, y=depthb, label=genes))
   if(logT){
  ggp<-ggp+scale_y_continuous(trans='log10',  limits = c(extra,10*max(merged[,depth_inds[k1]],na.rm=T)))+ylab(paste(type_nme[k1], label, extra))
  ggp<-ggp+scale_x_continuous(trans='log10' , limits = c(extra,10*max(merged[,depth_inds[k]],na.rm=T)))+xlab(paste(type_nme[k], label,extra))
   }
   ggps[[length(ggps)+1]] = ggp
   names(merged)[depth_inds[c(k,k1)]] = nme_old
   }
 }
 ml<-marrangeGrob(ggps,nrow = length(ggps), ncol = 1) 	
 
 ml
}

plotHMClust<-function(dfs){
  ggp = list()
  for(i in 1:length(dfs)){
  p <- ggplot(df, aes(x=region, y=mean, fill=s_e)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(.9))
  
  p <- p+ scale_fill_brewer(palette="Paired") + theme_minimal()+scale_y_continuous(trans='log10')+ggtitle(names(breakPs))
  ggp[[i]] <-p
}
   ml<-marrangeGrob(ggp,nrow = length(breakPs), ncol = 1) 	
  ml
  
}

expandBr<-function(breakP_, region){
  startc = region[1:3]
  endc = region[4:6]
  rows = breakP_$rows
  cols = breakP_$cols
  
  rowl = floor(seq(startc[1], startc[2], by = startc[3]))
  coll = floor(seq(endc[1], endc[2], by = endc[3]))
  if(rowl[length(rowl)]<region[2]) rowl = c(rowl, rowl[length(rowl)] +startc[3])
  if(coll[length(coll)]<region[5]) coll = c(coll, coll[length(coll)]+endc[3])
  rowToBin = apply(rows[,1,drop=F], 1, function(v, br) min(which(br>=v[1])), rowl)
  colToBin = apply(cols[,1,drop=F], 1, function(v, br) min(which(br>=v[1])), coll)
  heatm = breakP_$heatm
  #	
  heatm1 = array(0,dim = c(length(rowl), length(coll)), dimnames = list(rowl, coll))
  if(dim(rows)[1]>1 & dim(cols)[1]>1){
  for(i in 1:dim(rows)[1]){
    i1 = rowToBin[i]
    for(j in 1:dim(cols)[1]){
      j1 = colToBin[j]
      #print(paste(i1,j1))
      heatm1[i1,j1]  = heatm1[i1,j1] + heatm[i,j]
    }
  }
  
}
  vol = startc[3] *endc[3]
  heatm1=apply(heatm1,c(1,2),function(x) x/vol) 
  dimh = dim(heatm1)
  rows1 = data.frame(pos = as.numeric(dimnames(heatm1)[[1]]), depth = apply(heatm1, 1, sum) , s_e = as.factor(rep(as.character(rows$s_e[1]), dimh[1])))
  cols1 = data.frame(pos = as.numeric(dimnames(heatm1)[[2]]), depth = apply(heatm1, 2, sum) , s_e = as.factor(rep(as.character(cols$s_e[1]), dimh[2])))
  res = list(heatm = heatm1, rows = rows1, cols = cols1)
  res
}

#third element of startc and endc is bin
getHMCol<-function(len = 250){
  rev(c(rev(colorRampPalette(brewer.pal(9,"Reds"),bias = 0.5)(len)), colorRampPalette(brewer.pal(9,"Blues"), bias = 0.5)(len)))
}
.findBreaks2<-function(pos,mins){

  res_a = c()
  if(is.null(mins)) return(res_a)
  for(k in 1:dim(mins)[2]){
    min = mins[,k]
    step = pos[2] - pos[1]
    res = rep(NA, length(min))
    for(i in 1:length(min)){
      inds = which(min[i]>=pos & min[i]<pos+step)
      if(length(inds)>=1) res[i] = inds[1]
      
    }
    res[!is.na(res)]
    res_a = c(res_a, res)
  }
  unique(res_a)
}

blankGraph<-function(xlim, xax, yax, title = "" ) {
  ggplot(data.frame()) + geom_point()  + ylim(0, 100) + theme_bw()+xlab(xax)+ylab(yax)+ggtitle(title) +scale_x_continuous(limits = xlim[1:2])
 
}


plotBreakPByGenes<-function(breakP1, t, genesLeft="leader", genesRight="N",
                            left_range=NULL, right_range=NULL, step_left = 10, step_right=100, plotHM=T,
                            logT=T, title="", subtitle="", mult=1, fimo=NULL){
  region =  c(1,100,step_left,28000,30000,step_right)
  if(!is.null(genesLeft)){
    i1 = which(t$gene==genesLeft)
    region[1] = t$Minimum[i1]
    region[2] = t$Maximum[i1]
  }
  if(!is.null(genesRight)){
    i2 = which(t$gene==genesRight)
    region[4] = t$Minimum[i2-1]
    region[5] = t$Minimum[i2]
  }
  if(!is.null(left_range)){
    region[1] = left_range[1]
    region[2] = left_range[2]
  }
  if(!is.null(right_range)){
    region[4] = right_range[1]
    region[5] = right_range[2]
  }
  plotBreakPIntrons(breakP1,t,fimo=fimo,region=region, mult=mult, plotHM=plotHM, logT=logT, title=title, subtitle=subtitle)
                    
}

#plotBreakPIntrons<-function(breakP1, t=NULL, fimo=NULL,
 #                           region =  c(1,100,10,28000,30000,100), mult = 1, plotHM=T, logT = T, title = "", subtitle = ""){


findMaxSeqs<-function(breakP1,region=c(60,80,1,28240,28260,1), fasta, nme=""){
	left=c(-5,5); right = c(-5,5)
	if(nme=="scores33") {
		left =c(0,10); right = c(0,10)
	}
	if(nme=="scores55") {
		left =c(-10,0); right = c(-10,0)
	}
if(nme=="scores53") {
		left =c(-10,0); right = c(0,10)
	}
if(nme=="scores35") {
		left =c(0,10); right = c(-10,0)
	}
	breakP_ = subsetBr(breakP1, region)
 	breakP = expandBr(breakP_, region)  
	mm = breakP$heatm
maxv = max(mm)
	rowcol = which(mm == maxv, arr.ind = TRUE)
	res = list()
	nmes = rep("a", dim(rowcol)[1])
	for(k in 1:dim(rowcol)[1]){
		colind = rowcol[k,2]
		rowind = rowcol[k,1]
		row = breakP$rows[rowind,1]
		col = breakP$cols[colind,1]
		seq1 = fasta[[1]][(col+right[1]):(col+right[2])]
		seq2 = fasta[[1]][(row+left[1]):(row+left[2])]
		nmes[k] = paste(breakP$rows[rowind,1], breakP$cols[colind,1],maxv, sep="_")
		res[[k]] = list(seq1= seq1, seq2 = seq2, rowcol, row=breakP$rows[rowind,], col=breakP$cols[colind,], left=left, right = right, nme=nme)
        }
	names(res) =nmes
	res
}

plotBreakPIntrons<-function(breakP1, t=NULL, fimo=NULL, region =  c(1,5000,100,25000,30000,100), mult = 1, plotHM=T, logT = F, title = "", subtitle = ""){
  breakP_ = subsetBr(breakP1, region)
  col = getHMCol()
  breakP = expandBr(breakP_, region)  #bin = c(startc[3], endc[3]))
  maxv = max(breakP$heatm)
  if(maxv>0){
    heat_breaks = seq(0,maxv, length.out = length(col ))
    if(logT)   heat_breaks = seq(0,log10(maxv), length.out = length(col ))
  }
  title2  = paste(title, subtitle)
  print(title2)
  if(length(breakP$rows)>=2 && length(breakP$cols)>=2 && plotHM){
    if(logT){
      print("log trans")
      heatm = apply(breakP$heatm,c(1,2), function(x) if(x==0) NA else log10(x*mult))
    }else{
      heatm = apply(breakP$heatm,c(1,2), function(x)  x*mult)
    }
    
    rowSideCols = getSideCols(breakP$rows$pos, t)
    colSideCols = getSideCols(breakP$cols$pos, t)

    colsep = NULL
    rowsep = NULL
    inds = 1:length(rowSideCols);
    
     if(!is.null(fimo)){
     
      start_ind = grep("start",names(fimo))
      stop_ind = grep("stop",names(fimo))
      plus_strand = fimo$strand=="+"
      rowsep = .findBreaks2(breakP$rows$pos, fimo[plus_strand,start_ind,drop=F])

      colsep = .findBreaks2(breakP$cols$pos,  fimo[plus_strand,start_ind,drop=F])
     }
	nonNA = apply(heatm, 1, function(x) length(which(!is.na(x))))
	
      if(length(which(nonNA>1))>2){
        heatmap.2(heatm, scale="none" ,col = col,  na.color="LIGHTGREY", rowsep=rowsep,colsep = colsep, dendrogram = "none", trace ="none", Rowv = NULL, Colv = NULL, RowSideColors = rowSideCols, ColSideColors = colSideCols, main = title2)
      }
        
    
  }
  df = rbind( breakP$rows, breakP$cols)
  #if(dim(df)[1]==0) return (NULL)
  depth_str = "depth"
  depth_ind_ = which(names(df)=="depth")
#if(logT && FALSE){
 # depth_ind = which(names(df)=="depth")
  
	# df[,depth_ind]=df[,depth_ind]+0.1
	
  #}
  if(abs(mult-1)>1e-5){
    depth_str = "tpm"
    df$depth = df$depth *mult
    names(df)[depth_ind_] = depth_str;
  }
  #df$s_e = as.factor(df$s_e)
#  print(dim(df))
 # if(dim(df)[1]==0) return (list())
  startc = region[1:3]
  endc = region[4:6]
  #print(df)

  if(length(which(df$s_e=="start"))==0 || max(df[df$s_e=="start", depth_ind_], na.rm=T)<=0 ){
    ggp1<-blankGraph(startc[1:2],"pos", depth_str, title = title2);
  }else{
    ggp1<-ggplot(df,aes_string(x="pos",y=depth_str,fill =  "s_e" , colour = "s_e" ))+geom_point() + theme_bw()+ggtitle(title2)
    ggp1<-ggp1+geom_vline(xintercept = t$Minimum, linetype="solid", color=t$sideCols, alpha =0.5)
    ggp1<-ggp1+geom_vline(xintercept = t$Maximum, linetype="dashed", color=t$sideCols, alpha =0.5)
    if(!is.null(fimo)){
    	ggp1<-ggp1+geom_vline(xintercept = fimo$start[(fimo$strand=="+") & (fimo$motif_id == 'TRS_short')], linetype="dotted", color="black", alpha = 0.5)
    	ggp1<-ggp1+geom_vline(xintercept = fimo$stop[(fimo$strand=="+") & (fimo$motif_id == 'TRS_short')], linetype="dotted", color="grey" , alpha=0.5)
		
		ggp1<-ggp1+geom_vline(xintercept = fimo$start[(fimo$strand=="+") & (fimo$motif_id == 'TRS_long')], linetype="dotdash", color="black", alpha=0.5)
    	ggp1<-ggp1+geom_vline(xintercept = fimo$stop[(fimo$strand=="+") & (fimo$motif_id == 'TRS_long')], linetype="dotdash", color="grey", alpha=0.5)
		
		
    }
    if(logT) ggp1<-ggp1+scale_y_continuous(trans='log10') 
    ggp1<-ggp1 + scale_x_continuous(limits = startc[1:2])+theme(legend.position = "none")
  }
  if(length(which(df$s_e=="end"))==0 ||max(df[df$s_e=="end",depth_ind_], na.rm=T)<=0  ){
   # print(endc)
    ggp2<-blankGraph(endc[1:2],"pos", depth_str, title = title2);
  }else{
    ggp2<-ggplot(df,aes_string(x="pos",y=depth_str,fill =  "s_e" , colour = "s_e" ))+geom_point() + theme_bw()+ggtitle(title2)
    ggp2<-ggp2+geom_vline(xintercept = t$Minimum, linetype="solid", color=t$sideCols)
    ggp2<-ggp2+geom_vline(xintercept = t$Maximum, linetype="dashed", color=t$sideCols)
 if(!is.null(fimo)){
    ggp2<-ggp2+geom_vline(xintercept = fimo$start[(fimo$strand=="+") & (fimo$motif_id == 'TRS_short')], linetype="dotted", color="black", alpha=0.5)
    ggp2<-ggp2+geom_vline(xintercept = fimo$stop[(fimo$strand=="+") & (fimo$motif_id == 'TRS_short')], linetype="dotted", color="grey", alpha=0.5)
	
    ggp2<-ggp2+geom_vline(xintercept = fimo$start[(fimo$strand=="+") & (fimo$motif_id == 'TRS_long')], linetype="dotdash", color="black", alpha=0.5)
    ggp2<-ggp2+geom_vline(xintercept = fimo$stop[(fimo$strand=="+") & (fimo$motif_id == 'TRS_long')], linetype="dotdash", color="grey", alpha=0.5)

}
    if(logT) ggp2<-ggp2+scale_y_continuous(trans='log10') 
    ggp2<-ggp2 + scale_x_continuous(limits = endc[1:2])+theme(legend.position="none")
  }
  invisible(list(start=ggp1, end=ggp2))	
  
  #	invisible(list(plots = list(ggp1,ggp2), heatm = heatm))
  
}
.reorder<-function(nmes, split="_"){
nme_matrix = t(data.frame(strsplit(nmes, split=split)))
order(nme_matrix[,2], decreasing=T)
}



#	llHM(special, "special" , resdir,  breakPs, t, fimo, total_reads, type_nme, log=T)	

plotAllHM<-function(special, resname, resdir, breakPs,t,fimo, total_reads, todo = 1:length(breakPs),  type_nme=names(total_reads ), pdf=T, plotHM = T, plotStEnd = T, logT=T, depth=T){
  print(special)
 
  type_nme = type_nme[1:length(breakPs)]
  outfile1 = paste(resdir,"/", resname, "_breaks_heatmap.pdf",sep="")
  outfile2 = paste(resdir,"/", resname, "_start_end_breaks.pdf",sep="")
  if(plotHM && pdf) pdf(outfile1)
  plotsl = list()

  for(i1 in 1:length(special)){	
    #plots = list()
    for(i in todo){
      mult =  1e6/total_reads[i]
      if(depth)  mult=1
#	print(paste("here ", i1, i))
      plots_i = try(plotBreakPIntrons(breakPs[[i]], t,fimo,region =  special[[i1]],  mult =  mult, logT=logT, title= names(special)[i1], subtitle = type_nme[i], plotHM = plotHM))
	#plots_i
      if(inherits(plots_i,"try-error")) {
	print("error with plotBreakPIntrons")
      }else{
     	      names(plots_i) = paste(names(special)[i1], type_nme[i],names(plots_i), sep="_")
	      if(length(plots_i)>0){
			 plotsl = c(plotsl, plots_i)
	      }
}
    }
  }
   # types = unlist(lapply(plots, typeof));
   # print(types)
    #plots = plots[types=="list"]
   # if(length(plotsl)>0)	plotsl = c(plotsl, plots[.reorder(names(plots))])
    
  
  if(pdf && plotHM) dev.off()
  #types = unlist(lapply(plotsl, typeof));
  #plotsl = plotsl[types=="list"]
  ml<-marrangeGrob(plotsl,nrow = length(breakPs), ncol = 2) 	
  if(plotStEnd) try(ggsave(outfile2, plot=ml, width = 30, height = 30, units = "cm"))
  invisible(ml)
}


.addZero<-function(mat, thresh=10, toadd=0.01){
  len = dim(mat)[[1]]
  ncol = dim(mat)[[2]]
  diffs = apply(cbind(mat[-1,1], mat[1:(len-1),1]),1, function(v)v[1]-v[2])
  gaps=which(diffs>thresh)
  if(length(gaps)>0){
  mat[gaps,-1] =rep(0, ncol-1)
  mat[gaps+1,-1] =rep(0, ncol-1)
  }
  mat
}

.mergeDepthMats<-function(mats){
  pos =unique(sort(unlist(lapply(mats,function(x) x[,1]))))
  
 res=matrix(0,nrow=length(pos), ncol = dim(mats[[1]])[[2]])
 res[,1] = pos
 for(i in 1:length(mats)){
   inds=match(mats[[i]][,1],pos)
   res[inds,-1]=res[inds,-1]+mats[[i]][,-1]
 }
 res
}

.processInternal<-function(mat, sumAll,header,total_reads,span, gapthresh,ID, sumID="all", toAdd=0.01){
   mat =  .addZero(mat, thresh=gapthresh,toAdd)
  if(!sumAll && !is.null(total_reads)){
      mat = t(apply(mat,1,function(v)v/c(1,total_reads)))
  }
  mat = data.frame(mat)
  names(mat) = header
  if(sumAll){
    mat[,2] = apply(mat[,-1,drop=F],1,sum)
    mat = mat[,1:2]
    names(mat)[2] = sumID
  }
  mat1 = loess_smooth(mat, 2:dim(mat)[2], span)
  clusterID = rep(ID,dim(mat1)[[1]])
#  print(head(mat1))
  cbind(clusterID,mat1)
}
readH5<-function(h5file, total_reads, header, toplot, path="depth",gapthresh=100,mergeGroups = NULL,
                 sumID="all", pos = NULL,id_cols = c("molecule","cell","time"), toAdd=0,dinds  = 2*(2:length(header)-2)+2,  span =0.0, cumul= if(!is.null(pos)) F else T, sumAll=F){
 pos_ind = 1
 merge=!is.null(mergeGroups)
 ncols = length(id_cols)
 names = h5ls(h5file)$name
 inds = which(toplot %in% names)
 if(length(inds)==0) return (NULL)
 IDS = toplot[inds]
 clusters_ = NULL
 mats=NULL
 new_cols =  c("pos", "depth", "clusterID",'sampleID', id_cols)
 clusters_ = data.frame(matrix(NA, nrow =0, ncol = length(header)))
 names(clusters_) =header
 if(merge){
   mats <- vector("list", length(IDS))
 }
  for(i in 1:length(IDS)){
	ID = IDS[i]
	
	mat = t(h5read(h5file,paste(path,as.character(ID),sep="/")))
	#print(head(mat))
	if(dim(mat)[1]>0){
	mat=mat[,c(1,dinds),drop=F]
	if(merge){
	  mats[[i]] = mat
	}else{
	  mat2 = .processInternal(mat, sumAll, header, total_reads, span, gapthresh,ID, sumID=sumID, toAdd=toAdd)
		clusters_ = rbind(clusters_,mat2)
	}
	}
  } 
 if(merge){
   for(k in 1:length(mergeGroups)){
     indsk=which(IDS %in% mergeGroups[[k]])
   #  print(indsk)
     if(length(indsk)>0){
      
       
        matsk = mats[indsk]
        nonnull=!unlist(lapply(matsk, is.null)) 
        if(length(which(nonnull))>0){
          mat2 = .mergeDepthMats(matsk[nonnull])
          clusters_=rbind(clusters_,.processInternal(mat2,sumAll, header, total_reads, span, gapthresh,names(mergeGroups)[k], sumID=sumID, toAdd=toAdd))
        }
     }
   }
   
 }
  clusters_
}
plotHeatmap<-function(h5file, header,  transcripts1, jk, logHeatmap = F,xlim = list(c(1,max(clusters$pos+1))),featureName = "depth", 
                      max_h=0 ,title = ""){
  #len = max(clusters[,1]+1)
  if(max_h==0) stop("error must set max_h")
  title = paste(title, featureName)
#countdf = grep(featureName, names(transcripts1))
 IDS = transcripts1$ID
cnames =transcripts1$ORFs;
# apply(cbind(as.character(transcripts1$leftGene), as.character(transcripts1$rightGene)), 1, paste, collapse=".")
  for(jj in 1:length(xlim)){
    print(paste('jj',xlim[[jj]]))
    xmin = xlim[[jj]][1]
    xmax = xlim[[jj]][2]
    
  x = xlim[[jj]][1]:xlim[[jj]][2]
  nrows = dim(transcripts1)[1]
 # defaultVal = NA
  if(featureName=="ratios")defaultVal = NA else defaultVal = 0
  heatm =  matrix(defaultVal,nrow = nrows, ncol=  length(x)) 
  dimnames(heatm) = list(cnames, x)
 dinds = grep(featureName, header)

  for(i in 1:length(IDS)){
    matr = try(h5read(h5file,as.character(IDS[i])))
    if(!inherits(matr,"try-error")) {
	    subm = data.frame(t(matr))
	    names(subm) = header
	    subm = subm[subm$pos>=xmin & subm$pos<=xmax,,drop=F]
	    if(length(dinds)==0){
		ratio =subm[,grep("error", header)[jk]] / subm[,grep("depth", header)[jk]]
		ratio[is.na(ratio) | ratio>1000] = defaultVal
	 	heatm[i,subm$pos-xmin+1] = ratio
	    } else{
		 
	      heatm[i,subm$pos-xmin+1] = subm[,dinds[jk]]
	    }
    }else{
	print("could not find ")
	prnt(transcripts1[i,])
    }
  }
  
  ccounts = transcripts1$countTotal
  cols =getHMCol(10)
 
  if(logHeatmap){
    clustersLog =  apply(heatm, c(1,2), function(x) if(x==0) NA else log10(x))
    logcc = log(ccounts)
    max_h = log10(max_h)
  }else{
    clustersLog = heatm
    logcc = ccounts
    if(featureName=="ratios"){
      logcc = log(ccounts)
      max_h = log10(max_h)
      
      
    }
  }
#if(max_h<=0) stop("error with max_h")
  breaks = seq(0, max_h, length.out = length(cols))
  breaks = c(breaks, breaks[2] + breaks[length(breaks)])
 cols1 = hcl.colors(12, "YlOrRd", rev = TRUE)
  br =  seq(0, max_h, length.out = length(cols1))
  colSidecols = rep(NA, length(logcc))
  for(i in 2:length(br)){
    colSidecols[which(logcc>=br[i-1] & logcc<br[i])] = cols1[i-1]
  }
  colSidecols[which(logcc>=br[length(br)-1])]   = cols1[length(br)-1]
  # print(colSidecols)
  colnames(clustersLog)  = x
  rowSideCols = getSideCols(x, t)
  if(nrows>1 && length(which(!is.na(clustersLog)))>2){
    heatmap(clustersLog[nrows:1,,drop=F],scale="none", Colv = NA,margins = c(5,10), Rowv =  NA,distfun = distbin, main = title, ColSideColors = rowSideCols) #RowSideColors = colSidecols[nrows:1])
  }
  #	
  }
  #na.color="GREY",
#  invisible(heatm)
}



findRows<-function(clusters1, start, end){
  x = as.numeric(dimnames(clusters1)[[2]])
  crit1 = apply(clusters1[,which(x<start & x> 20000)],1,max)==0
  crit2= apply(clusters1[,which(x>=start & x<=end)],1,max)>10
  which(crit1 & crit2)
}




checkF<-function(clusters, trans1, i, offl = 0, offr = 0){
  subc = clusters[(trans1$startPos[i]+1-offl) : (trans1$endPos[i]+offr),]
  print(head(subc))
  print(tail(subc))
  
}

findBreaks<-function(row, mult=10){
  zeros = (row==0)
  
  v = which(zeros[-c(1, length(zeros))])
  start = v[1]+1
  en = which(row[start:length(zeros)]>0.1*max(row))[1]+start-2
  
  res = c(start,en)*mult+1
  c(res, row[start], row[end])
}
readCoords<-function(coords_file){
  ###COREF ANALYSIS (heatmap)
  t = read.csv(coords_file)
  levs = levels(as.factor(t$gene))
  mygroup<-1:length(levs)
  cols = brewer.pal(length(mygroup), "Set3")
  sideCols = cols[match(t$gene, levs)]
  
  t = cbind(t, sideCols)
  t
}

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}


.findIndex<-function(x,v){
v1 = abs(v - x)
ind = which(v1==min(v1))[1]
c(ind,v1[ind])
} 

#"5_3"     "5_no3"   "no5_3"   "no5_no3"
plotLengthHist<-function(reads,t1,   seqlen, type_nme, m_thresh = 100, start_thresh = c(T,T,F,F), end_thresh = c(T,F,T,F), min =0, max =10000, logy = F, binwidth = 10){
ggps = list()
len = t1$readlen
perc = rep(c(0.95,1.0), length(len)/2)
for(j in 1:length(type_nme)){
	
	for(k in 1:length(start_thresh)){
		title = paste(type_nme[j], nmes[k])
		subinds = if(start_thresh[k]) reads$startPos<m_thresh else reads$startPos>m_thresh
		subinds= subinds & ( if(end_thresh[k]) reads$endPos> seqlen-m_thresh else reads$endPos<=seqlen-m_thresh)
		subreads = reads[subinds  & reads$source==(type_nme[j]),,]
		ggp<-ggplot(subreads, aes(length)) + geom_histogram(binwidth = binwidth) + theme_bw()+ggtitle(title)+scale_x_continuous(limits = c(min, max))
		if(logy) ggp<-ggp+scale_y_continuous(trans='log10', label=scientific_10) else ggp<-ggp+scale_y_continuous(label=scientific_10)
		if(k==1 || k==3){
		maxy = ggplot_build(ggp)$layout$panel_scales_y[[1]]$range$range[2] 

		for(m in 1:length(len)){
			ggp<-ggp+ annotate(size = 2, "text", x = len[m], y = perc[m]*maxy, label = as.character(t1$gene[m]))
		}
		}
 		#ggp<-ggp+geom_text(  size=2, aes(x=len, y=depthb, label=gene))
		ggps[[length(ggps)+1]]<- ggp
	}

}

  ml<-marrangeGrob(ggps,ncol = length(type_nme), nrow = length(nmes)) 	
invisible(ml)
}


.getRatios<-function(clusters, dinds, einds){
	
     ratios=  data.frame(matrix(nrow = dim(clusters)[1], ncol = length(dinds)))
    for(i in 1:length(dinds)){
      ratios[,i] = clusters[,einds[i]]/clusters[,dinds[i]]
      ratios[ clusters[,dinds[i]]==0,i] = NA
    }
     names(ratios) = paste("ratios",1:length(dinds), sep="")
  ratios
}



.readFastqHeader<-function(f,
	nme =c("readID", "clusterID", "subID", "readSt", "readEnd", "readLen", "refStart", "refEnd", "refLen", "stLeft", "stRight", "endLeft", "endRight" ),
colclasses = c("character", "character", "character", rep("numeric", 6), rep("character", 4)))
{
	p = pipe(paste("sed -n '1~4p'", f))
	t = (read.table(p, sep=" ", header=F, colClasses = colclasses, fill=T))
	t = t[t[,2]!="",]
	names(t) = nme
t
}

.filterT<-function(transcripts, seqlen,  start_end = NULL, lt = c(T,T), leftGene= NULL, rightGene = NULL){
	inds = rep(T, dim(transcripts)[1])
	if(!is.null(start_end)){
		inds = inds & ( if(lt[1]) transcripts$start<start_end[1] else transcripts$start >= start_end[1])
		inds = inds & (if(lt[2]) (seqlen - transcripts$end<start_end[2]) else (seqlen - transcripts$end) >= start_end[2])
	}
	if(!is.null(rightGene)) inds = inds & transcripts$rightGene==rightGene
	if(!is.null(leftGene)) inds = inds & transcripts$leftGene==leftGene
	transcripts[inds,,drop=F]
}

.readTranscripts<-function(infilesT ){
inf = scan(infilesT, nlines=1, what=character())
transcripts = read.table( infilesT,sep="\t", head=T)
#countAll = apply(transcripts[,grep("count", names(transcripts))],1,sum)
#comb_l = apply(cbind(as.character(transcripts$leftGene), as.character(transcripts$rightGene)),1,function(x) paste(x, collapse="."))
#comb_l2 = apply(cbind(as.character(transcripts$leftGene2), as.character(transcripts$rightGene2)),1,function(x) paste(x, collapse="."))
#transcripts = cbind(transcripts, comb_l, comb_l2, countAll)
o = order(transcripts$countT, decreasing=T)
ORFs = as.character(transcripts$ORFs)
ORFS_uniq = unique(ORFs[o])
transcripts$ORFs = factor(ORFs,levels=ORFS_uniq, labels = ORFS_uniq)
transcripts = transcripts[o,]
err_ratio_inds = grep("error_ratio", names(transcripts))
transcripts[,err_ratio_inds] =apply(transcripts[,err_ratio_inds,drop=F], c(1,2), function(x) if(is.na(x)) -0.01 else x)
if(length(grep("#", inf))>0) attr(transcripts,"info") = sub("#", "",inf)
transcripts
}
.expand<-function(subs,nme="sample"){
  i = which(names(subs)==nme)
  molecule_type = factor(unlist(lapply(as.character(subs[,i]), function(x) strsplit(x,"_")[[1]][1])))
  cell = factor(unlist(lapply(as.character(subs[,i]), function(x) strsplit(x,"_")[[1]][2])))
  time = unlist(lapply(as.character(subs[,i]), function(x) strsplit(x,"_")[[1]][3]))
  time =  factor(time,level= paste(sort(as.numeric(unique(sub("hpi","",time)))),"hpi",sep=""))
  mat1 = cbind(as.character(molecule_type),as.character(cell),as.character(time))
  cell_time =.combine(cell,time)
  mol_time = .combine(molecule_type,time)
  cell_mol = .combine(molecule_type,cell)
  #ord = order(as.numeric(factor(types1_$time, levels=c("0hpi", "2hpi","24hpi","48hpi"))),types1_$cell,types1_$molecules)
 # print("h")
#  print(cell_time)
  cbind(subs,cell, molecule_type, time,cell_time,mol_time,cell_mol)
}
.combine<-function(f1,f2, rev=T){
  vs = apply(cbind(as.character(f1), as.character(f2)),1,paste,collapse="_")
  levs = if(rev)  unique(vs[order(f2,f1)]) else unique(vs[order(f1,f2)])
  factor(vs, levels=levs)
}

.splitTranscripts<-function(transcripts, seqlen, nmes,splice = F){
	transcripts_all = list()
	break_ind = grep("num_", names(transcripts));
	num_breaks = transcripts[,break_ind];
	if(splice){
		lead_inds =transcripts$ORFs %in% sort(unique( c(grep(";leader", transcripts$ORFs,v=T), grep("start", transcripts$ORFs,v=T))))
		transcripts_all[[1]] = 	transcripts[num_breaks==0,,drop=F]
		transcripts_all[[2]] = 	transcripts[num_breaks==1 & lead_inds,,drop=F]
		transcripts_all[[3]] = 	transcripts[num_breaks==1 & !lead_inds,,drop=F]
		transcripts_all[[4]] = 	transcripts[num_breaks>=2 & lead_inds,,drop=F]
		transcripts_all[[5]] = 	transcripts[num_breaks>=2 & !lead_inds,,drop=F]
		names(transcripts_all) = c("unspliced","spliced1_leader","spliced1_noleader","spliced2_leader","spliced2_noleader");
	
	}else{
		transcripts_all[[1]] = (transcripts[which(transcripts$start<100 & transcripts$end > seqlen -100),, drop=F])
		transcripts_all[[2]] = (transcripts[which(transcripts$start<100 & transcripts$end <= seqlen -100),, drop=F])
		transcripts_all[[3]] =( transcripts[which(transcripts$start>=100 & transcripts$end > seqlen -100),, drop=F])
		transcripts_all[[4]] = (transcripts[which(transcripts$start>=100 & transcripts$end <= seqlen -100),, drop=F])
		names(transcripts_all) = nmes
	}
num_clusters = lapply(transcripts_all, function(x) dim(x)[1])
transcripts_all = transcripts_all[num_clusters>0]
	transcripts_all
}

.getCompareVec<-function(type_nme, b=1){
	if(length(type_nme)==1){
		tocompare = list(c(1,1))
	}else if(length(type_nme)==4){
		tocompare = list(c(b,1), c(b,2),c(b,3), c(b,4))[-b]
	}else{
		tocompare = list();
		v = 1:length(type_nme)
		v = v[-b]
		for(i in 1:length(v)){
			tocompare[[i]] = c(b,v[i])
		}
	}
	return(tocompare)
}


.sumByLevel<-function(reads, nme1=c(  "ORFs"), target="source", mincount = 0.1, limit = 20){
	
	target_ind = which(names(reads) %in% target)
	inds1 = which(names(reads) %in% nme1)
	fact1 = as.factor(apply(reads[,inds1,drop=F], 1, function(x) paste(as.character(x), collapse=".")))
	levs1 = levels(fact1)
	matr = reads[1:length(levs1),]
	src_lev = levels(reads[,target_ind])
	counts = matrix(nrow =length(levs1), ncol = length(src_lev))
	names(counts) = as.character(src_lev)
	for(i in 1:length(levs1)){
	 subm =  reads[fact1==levs1[i],,drop=F]
	 for(j in 1:length(src_lev)){
		counts[i,j] = length(which(subm[,target_ind]==as.character(src_lev[j])))
	} 
	 matr[i,] = subm[1,]
	}

	counts = data.frame(apply(counts,c(1,2), function(x) x+mincount))
	names(counts) = as.character(src_lev)
	matr = cbind(matr, counts)
	matr[apply(counts,1,sum)>limit,,drop=F]
}

.plotGeneExpression<-function(reads_, nme1 = c("ORFs"), target = "source", limit = 20, mincount = 0.1, removeNA = T){
	reads1 = if(!removeNA) reads_ else reads_[!is.na(reads_$upstream) & !is.na(reads_$downstream),]
	#print(dim(reads1))
	sumByLevel = .sumByLevel(reads1, nme1=nme1 ,target=target,  mincount = mincount, limit = limit)

	extra = paste("+", mincount, sep="")
	lev = levels(reads_$source);
	#print(head(sumByLevel))
	ggp<-ggplot(sumByLevel, aes_string(x=lev[1], y=lev[length(lev)], fill = "ORFs", colour = "ORFs")) + geom_point(size = 4)  #+ theme_bw()
	ggp<-ggp+scale_y_continuous(trans='log10')+scale_x_continuous(trans='log10')
	ggp<-ggp+xlab(paste("Cell depth", extra)) + ylab(paste("Virion depth", extra))
	ggp<-ggp+theme(text = element_text(size=20))
 	#outfile0 = paste(resdir, "/rel_expression1.pdf", sep="");
	invisible(ggp)
}  
.bindall<-function(l){
	a = l[[1]]
	for(i in 2:length(l))	{
		 a = rbind(a, l[[i]])
	}
	data.frame(a)
}

plotHist<-function(vec, breaks, t = NULL, ylog=T, xlim = NULL, title= ""){
h = hist(vec, breaks,plot=F);
f<-function(h1, i, k) {
			 #norm  = if(normalise) sum(h1$count) else 1
			#norm = size[k]/1e6
			cbind(h1$breaks[-1], h1$count, rep(i, length(breaks)-1),rep(k, length(breaks)-1))
		}
a = data.frame(f(h,1,1))
	names(a) = c("breaks", "counts", "type", "source")
ggp <-ggplot(a, aes_string(x="breaks", y="counts", fill="source", color = "source"))+geom_point(size=2)+ggtitle(title)
if(!is.null(t) && !is.null(t$Minimum)) ggp<-ggp+ geom_vline(xintercept = t$Minimum, linetype="solid")  
		if(!is.null(t) && !is.null(t$start)) ggp<-ggp+ geom_vline(xintercept = t$start, linetype="solid")  
		if(!is.null(xlim)) ggp<-ggp + scale_x_continuous(limits = xlim)
		if(ylog)  ggp<-ggp+scale_y_continuous(trans='log10')
res = list(hist=h,plot=ggp)
}

plotJoins<-function(reads, thresh = c(10,100), binsize = 1, t = NULL, ylog = F, xlim = c(1,30000)){
	outfile2= paste("recomb", binsize, t,paste(xlim, collapse="."), "pdf", sep=".")
	print(outfile2)
	breaks = seq(1,seqlen+binsize, binsize)-0.5
	levs = levels(reads$source)
	size = unlist(lapply(levs, function(x) length(which(reads$source==x))))
	a_ = list()
	h_all = list()
	for(k in 1:length(levs)){
		source = levs[k]
	 	h = list(
			hist(reads[which(reads$start_read<thresh[1] & reads$source %in% source),]$startPos, br=breaks, plot = F),
			hist(reads[which(reads$start_read>thresh[2] & reads$source %in% source),]$startPos, br=breaks, plot = F),
			hist(reads[which((reads$length - reads$end_read)<thresh[1] & reads$source==source),]$endPos, br=breaks, plot = F),
			hist(reads[which((reads$length - reads$end_read)>thresh[2] & reads$source==source),]$endPos, br=breaks, plot = F)
		)
		names(h) = c("start_less", "start_more","end_less", "end_more" )
		h_all[[k]] = h	
		f<-function(h1, i, k) {
			 #norm  = if(normalise) sum(h1$count) else 1
			norm = size[k]/1e6
			cbind(h1$breaks[-1], h1$count/norm, rep(i, length(breaks)-1),rep(k, length(breaks)-1))
		}
		
		for(i in 1:length(h)) {
			print(paste(i,k))
			a_[[length(a_)+1]] = f(h[[i]], i,k)
		}

		
	}
	names(h_all) = levs
	maxpos = which(h_all[[1]][[2]]$counts == max(h_all[[1]][[2]]$counts))
	maxbr= h_all[[1]][[2]]$breaks[maxpos]
	indices = (which(reads$start_read > thresh & reads$startPos %in% (floor(maxbr):(1+floor(maxbr)))))
	subm = reads[indices,]
	a = .bindall(a_)	
	names(a) = c("breaks", "counts", "type", "source")
	a$type = factor(a$type, labels= c(paste("start less", thresh[1]), paste("start more", thresh[2]), paste("end less", thresh[1]), paste("end more", thresh[2])))
	a$source = factor(a$source, labels=levs)
	levst = levels(a$type)
	ggps = list()
	
	for(k in 1:length(levst)){
		inds_a = grep(levst[k], a$type)
		
		
		ggp <-ggplot(a[inds_a, ], aes_string(x="breaks", y="counts", fill="source", color = "source"))+geom_point(size=2)+ggtitle(levst[k])
		if(!is.null(t) && !is.null(t$Minimum)) ggp<-ggp+ geom_vline(xintercept = t$Minimum, linetype="solid")  
		if(!is.null(t) && !is.null(t$start)) ggp<-ggp+ geom_vline(xintercept = t$start, linetype="solid")  
		if(!is.null(xlim)) ggp<-ggp + scale_x_continuous(limits = xlim)
		if(ylog)  ggp<-ggp+scale_y_continuous(trans='log10')
#abline(v = t$Minimum, col=2)+abline(v = t$Maximum, col=3)
		ggps[[k]]<-ggp

	}

	  ml<-marrangeGrob(ggps,nrow = 2, ncol = 2, as.table=F) 
#ggsave(outfile2, plot=ml, width = 30, height = 30, units = "cm")
	invisible(ml)
}

plotJoins1<-function(subr, binsize = 1, t = NULL, ylog = F, xlim = c(1,30000), clusterId = NULL, right = FALSE){
if(!is.null(clusterId)) subr = subr[subr$clusterId==clusterId,,drop=F]
	#outfile2= paste("j", binsize, t,paste(xlim, collapse="."), "pdf", sep=".")
	breaks = seq(1,seqlen+binsize, binsize)-0.5
	#print(head(breaks))
	levs = levels(subr$source)
	#size = unlist(lapply(levs, function(x) length(which(subr$source==x))))
	a_ = list()
	h_all = list()
	for(k in 1:length(levs)){
		source = levs[k]
		indsk = if(right) subr$length- subr$end_read>thresh[2]  else subr$start_read>thresh[2] 
		
	 	h = list(
			hist(subr[indsk & subr$source %in% source,]$startPos, br=breaks, plot = F)
		)
		names(h) = c("start_less") #, "start_more","end_less", "end_more" )
		h_all[[k]] = h	
		f<-function(h1, i, k) {
			 #norm  = if(normalise) sum(h1$count) else 1
			norm = size[k]/1e6
			norm = 1
			cbind(h1$breaks[-1], h1$count/norm, rep(i, length(breaks)-1),rep(k, length(breaks)-1))
		}
		
		for(i in 1:length(h)) {
			print(paste(i,k))
			a_[[length(a_)+1]] = f(h[[i]], i,k)
		}

		
	}
	names(h_all) = levs
	#maxpos = which(h_all[[1]][[1]]$counts == max(h_all[[1]][[2]]$counts))
	#maxbr= h_all[[1]][[2]]$breaks[maxpos]
	#indices = (which(reads$start_read > thresh & reads$startPos %in% (floor(maxbr):(1+floor(maxbr)))))
	#subm = reads[indices,]
	a = .bindall(a_)	
	names(a) = c("breaks", "counts", "type", "source")
	a$type = factor(a$type) #, labels= c(paste("start less", thresh[1]), paste("start more", thresh[2]), paste("end less", thresh[1]), paste("end more", thresh[2])))
	a$source = factor(a$source, labels=levs)
	levst = levels(a$type)
	ggps = list()
	
	for(k in 1:length(levst)){
		inds_a = grep(levst[k], a$type)
		
		
		ggp <-ggplot(a[inds_a, ], aes_string(x="breaks", y="counts", fill="source", color = "source"))+geom_point(size=2)+ggtitle(levst[k])
		if(!is.null(t) && !is.null(t$Minimum)) ggp<-ggp+ geom_vline(xintercept = t$Minimum, linetype="solid")  
		if(!is.null(t) && !is.null(t$start)) ggp<-ggp+ geom_vline(xintercept = t$start, linetype="solid")  
		if(!is.null(xlim)) ggp<-ggp + scale_x_continuous(limits = xlim)
		if(ylog)  ggp<-ggp+scale_y_continuous(trans='log10')
#abline(v = t$Minimum, col=2)+abline(v = t$Maximum, col=3)
		ggps[[k]]<-ggp

	}

	#  ml<-marrangeGrob(ggps,nrow = 2, ncol = 2, as.table=F) 
#ggsave(outfile2, plot=ggps[[1]], width = 30, height = 30, units = "cm")
	#invisible(h_all)
	ggps[[1]]
	#h_all
}

.plotGeneExpr<-function(tocomp,transcripts_all, todo =1:length(transcripts_all),  extra = 0.1, mint_ = 2, maxt_ = 10, count_df = grep('count[0-9]', names(transcripts_all[[1]])),
 type_nme = attr(transcripts_all[[1]], "info")){
ggps = list()
maxv = 0
for(k in todo){
  maxk = max(transcripts_all[[k]][,count_df], na.rm=T)
  if(maxk>maxv) maxv = maxk
}
for(k in todo){
	mint = mint_
 	if(dim(transcripts_all[[k]])[1]>maxt_) mint = max(mint, transcripts_all[[k]]$countT[maxt_])
	cco = transcripts_all[[k]]$countT>mint
	tf_k = transcripts_all[[k]][cco,,drop=F]
	
	tf_k[,count_df] = tf_k[,count_df] + extra
	lev = names(tf_k)[count_df][tocomp]
	type_nme1 = type_nme[tocomp]
	ncol = 0
	i =1 
	{
	  st = i+1
	   for(j in st:length(lev)){
		#print(paste(k,i,j))
		ncol = ncol+1
		ggp<-ggplot(tf_k, aes_string(x=lev[i], y=lev[j], fill = "ORFs", colour = "ORFs")) +geom_point() #+ geom_point(size = 4)  #+ theme_bw()
		ggp<-ggp+scale_y_continuous(trans='log10', limits = c(extra,maxv))+scale_x_continuous(trans='log10', limits = c(extra,maxv))
		ggp<-ggp+xlab(paste(type_nme1[i], extra,sep="+")) + ylab(paste(type_nme1[j], extra, sep="+"))+ggtitle(names(transcripts_all)[k])
		ggp<-ggp+ geom_point(size = 1)
		ggps[[length(ggps)+1]] <- ggp
	  }
	}
}
#print(length(ggps))
nrow = floor(length(todo)/2)
  ml<-marrangeGrob(ggps,nrow = nrow, ncol = 2, as.table=F) 
invisible(ml)
}

.appendGenePosition<-function(reads_leader, t1){
genepos = rep(0, dim(reads_leader)[1])
genes = t1$gene
for(i in 1:length(genes)){
	inds = which(reads_leader$downstream==as.character(genes[i]))
	if(length(inds)>0){
		genepos[inds] = t1$Minimum[i]
	}
}
reads_leader = cbind(reads_leader, genepos)
reads_leader
}


plotErrorViolin<-function(reads,reads_no_leader ,  inds1 = NULL,  inds2 = NULL, x = "ORFs",  ord = "genepos", y = "errorRatio", fill = "source"){
	subm = if(is.null(inds1)) reads else reads[inds1 ,]
	if(!is.null(reads_no_leader)){
		reads1 =if(is.null(inds2)) reads_no_leader else reads_no_leader[inds2 ,] 
		levs = levels(as.factor(reads1$type_nme))
		x_ind = which(names(reads1)==x)
		ord_ind = which(names(reads1)==ord)
		fill_ind = which(names(reads1)==fill)
		reads1[,x_ind] = as.character(reads1[,x_ind])
		for(i in 1:length(levs)){
			indsi = which(reads1$type_nme==levs[i])
			if(length(indsi)>0){
				reads1[indsi,x_ind] = levs[i]
				reads1[indsi,ord_ind] = i
		        }
		}
		subm = rbind(subm, reads1)
	}
	err_inds = which(names(reads)==y)
	
	x1 = paste("reorder(", x, ",", ord,")", sep="") #downstream,genepos)"
	ggp<-ggplot(subm, aes_string(x=x1, y=y,  fill = fill))
	ggp<-ggp + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+ggtitle("Error vs transcript")+xlab("ORF")
	#ggp<-ggp+ylim(c(0, max(reads[,err_inds],na.rm=T)))
	ggp<-ggp+theme(text = element_text(size=18), axis.text.x = element_text(size = rel(0.5), angle = 90))
	invisible(ggp)
}
    								      			        		          		        																														                  					                    	          	  															   	  									 		               				      																																      			        												        														 						        						 															    		        																																      			        												        											 						        															 											    							      											     	      			 		    	 											    				      		         				 	 										    	      	        													      			      				    	      											      																				        			      							    				    	                          	                    			 		          		            																									                                                                                                                                                                                                                                                                		                          	    	      		             	    	 						    	    	    	                                                                                                                                                                                                                                 	          				            		  						                                                                                                                                                                                                    	      			 	           	                      	      	                                                                                         	                                                                  				    					    	    	                                                               		                                             		                                                                                                                                                                   		              																					 		      			      			      				      					                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     	                 										 				    	       	   	             	              						                                                                                                                                                                                                                                                                                                                                                                                         	                                                                                                                                                                        	                                                                        	       	 	        			                  									                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    	      									    				                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       