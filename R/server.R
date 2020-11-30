library(shiny)
library(reshape2)
library(ggplot2)
library(tidyr)
library(rhdf5)
library(RColorBrewer)
library(binom)
library(writexl)
library(shinyjs)
library(seqinr)
library(abind)
library(ggrepel)
#library(GGally)

source( "transcript_functions.R")
#source("shiny-DE.R")

basedir="../data"

toreplace1 = c("leader_ORF1ab,S_ORF1ab,ORF10_end","leader_leader,S_ORF1ab,ORF10_end")
toreplace2= c("leader_ORF1ab,ORF1ab_ORF1ab,ORF10_end","leader_leader,ORF1ab_ORF1ab,ORF10_end")





#toreplace=list(virion="RNA_virion_0hpi", whole_genome_mapped="RNA_vero_24hpi")
decodeFile = paste(basedir,"decode.txt",sep='/')
replace=read.table(decodeFile,sep="\t",head=F)
toreplace = replace[,2]
names(toreplace) = replace[,1]



reorder=T

#dirs = list.dirs(basedir,full.names=F, rec=T)
#dirs=dirs[which(unlist(lapply(dirs,function(x) file.exists(paste(basedir,x,"0.isoforms.h5",sep="/")))))]
#seldir=1
#currdir = paste(basedir,dirs[seldir],sep="/")
#datafile=paste(currdir,"0.isoforms.h5",sep="/")
#h5file=paste(currdir,"0.clusters.h5",sep="/")      

#isoInfo = .getIsoInfo(datafile,h5file, toreplace)
#total_reads = isoInfo$total_reads
#header=names(total_reads)
#info=.processInfo(isoInfo)
#t = readCoords(paste(currdir, "Coordinates.csv",sep="/"))
#fimo_file = paste(currdir,"fimo.tsv",sep="/")
#fimo = read.table(fimo_file, sep="\t", head=T)

info=NULL
header=NULL
isoInfo=NULL
total_reads=NULL
t=NULL
fimo=NULL
datafile=NULL
h5file = NULL





.getlevels<-function(type_nme, molecules, cells, times, reverseOrder=T){
  types_=data.frame(t(data.frame(strsplit(type_nme,"_"))))
  names(types_) = c("molecules","cell","time")
  inds1 =  which(types_$molecules %in% molecules & types_$cell %in% cells & types_$time %in% times)
  types1_ = types_[inds1,,drop=F]
  ord = order(as.numeric(factor(types1_$time, levels=c("0hpi", "2hpi","24hpi","48hpi"))),types1_$cell,types1_$molecules)
 if(reverseOrder) ord = order(types1_$molecules,types1_$cell,as.numeric(factor(types1_$time, levels=c("0hpi", "2hpi","24hpi","48hpi"))))
  
  levels1=type_nme[inds1][ord]
  attr(levels1,"inds1") = inds1
  levels1
}
#run_depth(h5file, total_reads=total_reads)


.getCIs<-function(subs,sample, total_reads1,method, showTPM=F,prefix="",suffix="", after=TRUE){
  if(nchar(prefix)>0){
   if(after) prefix = paste(prefix,"_",sep="") else prefix = paste("_",prefix,sep="")
  }
  inds1 =  match(sample,gsub(prefix,"", names(total_reads1)))
  count_head = paste("count",suffix,sep="")
#  print(count_head)
#  print(names(subs))
  count_ind = which(names(subs)==count_head)
#  print(count_ind)
  subs1 = cbind(subs[,count_ind],total_reads1[inds1])
#  print(subs1)
  nrows = dim(subs)[1]
  cis = matrix(NA,ncol=3,nrow =nrows )
  dimnames(cis)[[2]] = paste(c("TPM","lower","upper"),suffix,sep="")
  for(i in 1:nrows){
   # print(i)
    cis[i,1:3] = .calcPropCI(subs1[i,], method=method, conf.int = conf.int)
  }
  if(!showTPM){
    cis = cis *subs1[,2]/1e6
  #  yname="counts"
  }
  cis
}

.extractTPM<-function(datafile ,  total_reads,countsTotal, p_data){
  tr = NULL  
  xy = p_data$xy
  if(p_data$calcTPMFromAll && !is.null(countsTotal)){
    countsHostVirus1 =   countsTotal
    inds1 = match(names(total_reads) ,countsHostVirus1$sample)
    # print(inds1)
    tr = countsHostVirus1$count[inds1]
    tr = as.numeric(as.character(tr))
    names(tr) = names(total_reads)
    print(total_reads)
    print("new total reads")
    print(tr)
    total_reads = tr
  }
  header = names(total_reads)
  toplot = p_data$toplot
  if(p_data$usegrep){
    x1 = .findEntries(toplot,datafile,"/trans",p_data$tojoin);
  }else{
    x1= toplot
  }
  mat = t(data.frame( lapply(x1, .readIso, datafile, header, "/trans")))
  if(is.null(dim(mat))) mat = matrix(mat,nrow=1,ncol=length(header))
  if(p_data$merge){
    mat = matrix(apply(mat,2,sum),nrow=1,ncol=dim(mat)[2])
  }else if(nchar(p_data$merge_by)>0){ 
    ord=order(apply(mat,1,sum),decreasing=T)
    groups=.mergeGroups(x1,p_data$merge_by, ord=ord, max_trans = p_data$max_trans)
    
  }else if(p_data$group_by != 'No grouping'){
    groups = .getGroups(x1,p_data$group_by)
    toplot=names(groups)
    mat1 = matrix(NA, nrow = length(toplot), ncol  =dim(mat)[2])
    for(j in 1:length(groups)){
      indsj = which(x1 %in% groups[[j]])
      mat1[j,]=apply(mat[indsj,,drop=F],2,sum)
    }
    mat = mat1
  } else{
    
    if(length(x1)>p_data$max_trans){
      ord=order(apply(mat,1,sum),decreasing=T)
      mat = mat[ord[1:p_data$max_trans],,drop=F]
      toplot = x1[ord[1:p_data$max_trans]]
    }else{
      toplot=x1
    }
  }
  #print(toplot)
  cells = p_data$cells; times=p_data$times; molecules = p_data$molecules; reverseOrder = p_data$reverseOrder
  splitby_vec= p_data$splitby_vec
  
  if(p_data$xy){
    subs = list()
    for(i in 1:length(splitby_vec)){
      if(splitby=="molecules"){
        levs1=.getlevels(header,molecules[i], cells, times, reverseOrder)
      } else if(splitby=="cells"){
        levs1=.getlevels(header,molecules, cells[i], times, reverseOrder)
      }else{
        levs1=.getlevels(header,molecules, cells, times[i], reverseOrder)
        
      }
      subs[[i]] =.processTPM(mat, header, toplot, levels=levs1,split=T)
     # print(head(subs[[i]]))
      # subs[[i]]$sample=sub(molecules[i],"",subs[[i]]$sample)
    }
    by=c("ID","cell","time")
    if(p_data$splitby=="cells"){
      by=c("ID","molecule_type","time") 
    }else if(p_data$splitby=="times"){
      by=c("ID","molecule_type","cell")
    }
    subs = merge(subs[[1]],subs[[2]],by=by)
    sample=apply(subs[,2:3,drop=F],1,paste,collapse="_")
  } else if(p_data$barchart){
    levs1=.getlevels(header,molecules, cells, times, reverseOrder)
    subs =.processTPM(mat, header, toplot, levels=levs1,split=F)
    sample = subs$sample
  }else{
    tpm_df= .processTPM(mat, header, toplot,split=T)
    subs=subset(tpm_df, molecule_type %in% molecules & cell %in% cells & time %in% times)
    sample=apply(subs[,2:4,drop=F],1,paste,collapse="_")
  }
  
  if(xy){
    if(p_data$splitby=="times") after =FALSE else after=TRUE
    cis.x = .getCIs(subs,sample,total_reads[ grep(splitby_vec[1],names(total_reads))],method, showTPM, prefix=splitby_vec[1],suffix=".x", after=after)
    cis.y = .getCIs(subs,sample,total_reads[ grep(splitby_vec[2],names(total_reads))],method, showTPM, prefix=splitby_vec[2],suffix=".y", after=after)
    cis = cbind(cis.x,cis.y)
  }else{
    #        print(subs)
    #        print(total_reads)
    cis = .getCIs(subs,sample,total_reads,p_data$method, p_data$showTPM)
  }
  if(!is.null(subs$sample)){
    cell = factor(unlist(lapply(as.character(subs$sample), function(x) strsplit(x,"_")[[1]][1])))
    molecule_type = factor(unlist(lapply(as.character(subs$sample), function(x) strsplit(x,"_")[[1]][2])))
    time = unlist(lapply(as.character(subs$sample), function(x) strsplit(x,"_")[[1]][3]))
    time =  factor(time,level= paste(sort(as.numeric(unique(sub("hpi","",time)))),"hpi",sep=""))
    subs = cbind(subs,cell, molecule_type, time)
  }
  subs = cbind(subs,cis)
  subs
}
.plotTPMData<-function(subs,countsHostVirus,p_data, p_plot,yname){
  xy = p_data$xy
  logy=p_plot$logy
  facet = p_plot$facet
  showTPM=p_data$showTPM
  if(p_data$xy){
    colorby=names(subs)[1]
    shapeby=names(subs)[2] 
    fillby = names(subs)[3]
    if(length(levels(subs[[2]]))<length(levels(subs[[3]]))){
      shapeby=names(subs)[3] 
      fillby=names(subs)[2]
      
    }
    ylim = c(min(subs$TPM.x, subs$TPM.y),c(max(subs$TPM.y, subs$TPM.y)))
    ggp<-ggplot(subs, aes_string(x="TPM.x", y="TPM.y",ymin="lower.y", ymax="upper.y", 
                                 xmin="lower.x", xmax="upper.x"))
    ggp<-ggp+ggtitle(yname)+theme_bw()
    ggp<-ggp +geom_point(inherit.aes=T,aes_string(shape = shapeby,fill=fillby,color=colorby,size=10))
    #   ggp<-ggp +geom_point(inherit.aes=T,aes_string(shape = shapeby,color=fillby,size=2))
    
    if(p_plot$showCI){
      ggp<-ggp+geom_errorbar(colour="black")
    } #
    trans="identity"
    if(p_plot$logy){
      trans="log10"
    }
    ggp<-ggp+ scale_y_continuous(trans=trans,name=splitby_vec[2], limits=ylim)+ scale_x_continuous(limits = ylim,trans=trans,name=splitby_vec[1])
    ggp<-ggp+theme_bw()+theme(text = element_text(size=textsize))
  }else if(p_data$barchart){
    ORF="ID"
    y_text="TPM"
    if(p_data$stack){
      sublevs = levels(subs$ID)
      if(length(sublevs)==4 && sublevs[[1]]=="5_3"){
        subs$ID = factor(as.character(subs$ID), levels=rev(c("5_3","non5_3", "5_non3","non5_non3")))
      }
      
      ggp<-ggplot()
      levs_subs = levels(subs$ID)
      cols_subs =  brewer.pal(n = length(levs_subs), name = "Set2")
      names(cols_subs) = levs_subs
    #  print(head(subs))
     
      ggp<-ggp+geom_bar(data=subs,aes(x=time,y=TPM,fill=ID,color=ID),position="stack",stat='identity')
    
      ylim = layer_scales(ggp)$y$range$range
      # print(ylim)
      if(!is.null(countsHostVirus) ){
        scaling_factor = 100/ylim[2]
        countsHostVirus[3] = countsHostVirus[3] /scaling_factor
        types=c("Host","Virus","Sequin")
        linetype=c("dashed","twodash","solid")
        shape = c(1,2,3)
        
        ggp<-ggp+geom_point(data=countsHostVirus, aes(x=time, y=Reads, shape=Type))
        ggp<-ggp+geom_line(data=countsHostVirus, aes(x=time, y=Reads, linetype=Type, group=Type))
        ggp<-ggp+  scale_color_manual(values = cols_subs)
        ggp<-ggp+ scale_linetype_manual(values = linetype)+scale_shape_manual(values=shape)
        ggp<-ggp+xlab("Conditions")
        ggp<-ggp+ scale_y_continuous(
          name = yname,
          sec.axis = sec_axis(~.*scaling_factor, name="Proportion (%)"))
      }
    }else{
      ggp<-ggplot(subs, aes_string(x=ORF,y=y_text,fill="sample", colour='sample',ymin="lower" ,ymax="upper"))
      ggp<-ggp+ geom_bar(position=position_dodge(), aes_string(y="TPM"),stat="identity")
      if(p_plot$showCI){
        ggp<-ggp+geom_errorbar(position=position_dodge(width=0.9),colour="black")
      } #ggp<-ggp+geom_errorbar(aes_string(x=x1,ymin="lower", ymax="upper"), width=.2)#, position="dodge")
    }
    ggp<-ggp+theme_bw()+theme(text = element_text(size=p_plot$textsize), axis.text.x = element_text(size = rel(1.0), angle = p_plot$angle, hjust=1.0))
    
    #geom_bar(aes_string(x=x1, y="Ratio", fill = "type", colour = "type"),stat="identity", position = "dodge")
    
    
    if(p_plot$logy){
      ggp<-ggp+ scale_y_continuous(trans="log10",name=yname)
    }
    
    
    ggp<-ggp+ggtitle(yname)
    #ggp<-ggp+xlab("ORF")
  }else if(!xy){
    if(p_plot$showCI){
      ggp<-ggplot(subs, aes(x=time, y=TPM ,ymin=lower ,ymax=upper,group=interaction(molecule_type, cell, ID), color = cell, linetype=ID))
      if(p_plot$ribbon){
        ggp<-ggp+ geom_line()  + geom_point(inherit.aes=T,aes(shape = molecule_type,size=10))
        ggp<-ggp+geom_ribbon( linetype=2, alpha=0.1)
      }else{
        ggp<-ggp+ geom_line(position=position_dodge(width=0.1))  + geom_point(position=position_dodge(width=0.1),inherit.aes=T,aes(shape = molecule_type,size=10))
        ggp<-ggp+geom_errorbar(position=position_dodge(width=0.1)) #,colour="black")
      }
    }else{
      ggp<-ggplot(subs, aes(x=time, y=TPM ,group=interaction(molecule_type, cell, ID), color = cell, linetype=ID))
      ggp<-ggp+ geom_line()  + geom_point(inherit.aes=T,aes(shape = molecule_type,size=10))
    }
    
    ggp<-ggp+theme_bw()+theme(text = element_text(size=p_plot$textsize))
    # ggp<-ggp+theme_bw();#+ylim(c(min(subs$TPM, na.rm=T), max(subs$TPM, na.rm=T)))
    if(!showTPM)ggp<-ggp+ylab("Counts")
    # ggp<-ggp+ geom_errorbar(aes(linetype=molecule_type))
    if(logy){
      ggp<-ggp+ scale_y_log10()
    }
  }
  if(facet=="molecules_and_cells"){
    ggp<-ggp+facet_grid(molecule_type~cell)
    ggp<-ggp+facet_grid(cell~molecule_type)
  }else if(facet=="molecules"){
    ggp<-ggp+facet_grid(~molecule_type)
    
  }else if(facet=="cells"){
    ggp<-ggp+facet_grid(~cell)
    
  }else if(facet=="times"){
    ggp<-ggp+facet_grid(~time)
  }else if(facet=="molecules_and_times"){
    ggp<-ggp+facet_grid(molecule_type~time)
    ggp<-ggp+facet_grid(time~molecule_type)
  }else if(facet=="times_and_cells"){
    ggp<-ggp+facet_grid(cell~time)
    ggp<-ggp+facet_grid(time~cell)
  }
  
 # "molecules_and_times","times_and_cells"
  
  ggp
  
  
}

.process1<-function(plottype1,info){
  choices1 = info$choices1
  choices = info$choices
  ind1 = which(names(choices1)==plottype1)
  ind = which(names(choices)==plottype1)
  if(length(ind)==0){
    label=paste("Transcript",names(choices1)[ind1])
    ch = c("-",choices1[[ind1]])
  }else{
    label=paste("Transcript",names(choices)[ind])
    ch = c("-",choices[[ind]])
  }
  
  list(label=label, ch=ch)
}

# HERE IS THE SERVER PART
##OUTPUT IS PASSED TO THE UI
##INPUT PASSES IN INFORMATION
shinyServer(function(input, output,session) {
	#output$instructions <- renderPrint({
	#	print("Upload 0.transcripts.txt.gz file produced by npTranscript");
	#})
	
	#source("shiny-DE.R")
	
	source( "transcript_functions.R")
	source("shiny-DE.R")
	
 plot_depth<-function(tpm_df, total_reads=NULL,  toplot=c("leader_leader,N_end", "N_end"),combinedID="combined", 
                      gapthresh=100, mergeGroups=NULL,downsample = F, molecules="RNA",cells="vero",times=c('2hpi','24hpi','48hpi'), 
                      span = 0.01, sumAll=F, xlim=NULL, motifpos=list(),peptides=NULL, alpha=1.0,t= NULL,logy=T, 
                      showORFs = F,showWaterfall=FALSE,waterfallKmer=3,waterfallOffset=0,top10=10,textsize=20,
                      ci = 0.995, depth_thresh = 1000,
                      path="depth",seq_df = NULL, plotCorr=F, linesize=0.1, reverseOrder=F, calcErrors=F, fisher =F){
   if(path!="depth") calcErrors = FALSE
   if(is.null(tpm_df)){
     print(paste("could not read ",toplot))
     return (ggplot())
   }
   if(calcErrors){
     ggp<-.plotError1(tpm_df, pval_thresh = 1e-3, t1=t, xlim=xlim, motifpos=motifpos, logy=logy, ci=ci, alpha=alpha)
     #  ggp<- .makeCombinedArray(clusters_, errors_, xlim, downsample = downsample, thresh = depth_thresh,alpha = alpha,  ci = ci, max_num = 10,t=t, fisher = fisher,motifpos=motifpos)
     return(ggp)
   }else  if(plotCorr){
     ggp<-ggplot(data =tpm_df, aes(x=Var1, y=Var2, fill=value)) +  geom_tile()+ggtitle(path)
     ggp<-ggp+scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                   midpoint = 0, limit = c(-1,1), space = "Lab", 
                                   name="Pearson\nCorrelation") +
       theme_minimal()+ 
       theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                        size = textsize, hjust = 1),
             axis.text.y = element_text( 
               size = textsize))+
       coord_fixed()+xlab("")+ylab("")
     
     return(ggp)
   }
   if(showWaterfall && !is.null(seq_df)){
     cnt_df = tpm_df
     ggp<-ggplot(cnt_df, aes(kmers, fill = kmers))
     ggp<-ggp+ geom_rect(aes(x = kmers,xmin = id - 0.45, xmax = id + 0.45, ymin = end,ymax = start))+ggtitle(path)
     ggp<-ggp+scale_colour_manual(values = rainbow(dim(cnt_df)[1]))
     ggp<-ggp+theme(text = element_text(size=textsize))
   }else{
     if(sumAll) type_nme = "combined"
     rawdepth = T
     leg_size1=textsize
     show=T
     fill=F
     k = 1
     linetype="clusterID"
     colour="sampleID"
     if(sumAll){
       colour="clusterID"
       linetype="sampleID"
     }
     ylab="depth"
     if(!is.null(total_reads)) ylab="depth per million mapped reads"
     ggp<-plotClusters(tpm_df,seq_df, 4,  1, 
                       t,
                       motifpos,peptides,size=20,linesize=linesize,textsize=textsize,
                       rawdepth = T, linetype=linetype, colour=colour,
                       alpha=alpha, xlim = xlim,ylab=ylab , title =path, logy=logy, leg_size =textsize, show=show, fill =fill)
  
   }
   return(ggp)
 }
  
  run_depth<-function(h5file, total_reads=NULL,  toplot=c("leader_leader,N_end", "N_end"),combinedID="combined", 
                      gapthresh=100, mergeGroups=NULL,downsample = F, molecules="RNA",cells="vero",times=c('2hpi','24hpi','48hpi'), 
                      span = 0.01, sumAll=F, xlim=NULL, motifpos=list(),peptides=NULL, alpha=1.0,t= NULL,logy=T,
                      showORFs = F,showWaterfall=FALSE,waterfallKmer=3,waterfallOffset=0,top10=10,textsize=20,
                      ci = 0.995, depth_thresh = 1000,
                      path="depth",seq_df = NULL, plotCorr=F, linesize=0.1, reverseOrder=F, calcErrors=F, fisher =F){
    if(path!="depth") calcErrors = FALSE
    header =.getHeaderH5(h5file,toreplace)
    if(path=="depth"){
      dinds  = 2*(2:(length(header))-2)+2
    }else{
      dinds = 2:(length(header))
    }
    type_nme = header[-1]
    types_=data.frame(t(data.frame(strsplit(type_nme,"_"))))
    names(types_) = c("molecules","cell","time")
    inds1 =  which(types_$molecules %in% molecules & types_$cell %in% cells & types_$time %in% times)
    types1_ = types_[inds1,,drop=F]
    ord = order(as.numeric(factor(types1_$time, levels=c("0hpi", "2hpi","24hpi","48hpi"))),types1_$cell,types1_$molecules)
    levs=type_nme[inds1][ord]
    toAdd=0
    if(logy)toAdd=0.001
    facts =  apply(types1_,2,function(v) levels(factor(v)))
    same_inds = which(unlist(lapply(facts,length))==1)
    sumID='all'
    if(length(same_inds)>0){
      sumID = paste(unlist(facts[same_inds]),collapse="_")
    }
    id_cols = c("molecule","cell","time")
    tot_reads=NULL
    if(!is.null(total_reads)){
      
      tot_reads =  total_reads[inds1]/rep(1e6,length(inds1))
    }
    if(calcErrors){
      tot_reads = rep(1,length(tot_reads))  ## we dont want to correct for read depth for errors 
    }
    clusters_ = readH5(h5file,tot_reads, c("pos",header[inds1+1]),toAdd = toAdd, 
                       mergeGroups=mergeGroups,sumID=sumID, path=path,toplot,id_cols=id_cols, 
                       gapthresh=gapthresh, dinds = dinds[inds1], pos =NULL, span = span, cumul=F, sumAll=sumAll)
    errors_ = NULL
    if(calcErrors){
      offset_ = length(header)-1
      errors_ = readH5(h5file,tot_reads, c("pos",header[inds1+1]),toAdd = toAdd, mergeGroups=mergeGroups,
                       sumID=sumID, path=path,toplot,id_cols=id_cols, gapthresh=gapthresh, 
                       dinds = dinds[inds1]+1, pos =NULL, span = span, cumul=F, sumAll=sumAll)
      tpm_df<- .makeCombinedArray(clusters_, errors_, xlim, downsample = downsample, thresh = depth_thresh,alpha = alpha,  ci = ci, max_num = 10,t=t, fisher = fisher,motifpos=motifpos)
      return(tpm_df)
 
     }else  if(plotCorr){
    indsp = clusters_$pos <=xlim[2] & clusters_$pos >= xlim[1]
    df = clusters_[indsp,2:dim(clusters_)[2],drop=F]
    nrows = length(which(indsp))
   
    sums=apply(df[-1],2,sum)
    df = df[,c(1,1+order(sums, decreasing=T))]
    len = dim(df)[2]-1
    if(len>1){
      print(names(df))
      #reverseOrder=F
      levels1 = .getlevels(names(df)[-1],molecules, cells, times,reverseOrder)
      cor = cor(df[,match( levels1,names(df))])
     
      df2 = melt(cor)
      return(df2)
    }
  }
    if(is.null(clusters_)){
      print(paste("could not read ",toplot))
      return (NULL)
    }
    
    if(sumAll) levs=names(clusters_)[3]
    
    tpm_df = melt(clusters_,id.vars=c("clusterID","pos"), measure.vars=names(clusters_)[-(1:2)], variable.name="sampleID",value.name='count') %>%
      transform(sampleID=factor(sampleID,levels=levs))
 
    tpm_df$clusterID= factor(.fix(as.character(tpm_df$clusterID),toreplace1,toreplace2))
    
    
    if(sumAll) type_nme = "combined"
    rawdepth = T
    leg_size1=textsize
    show=T
    fill=F
    k = 1
    linetype="clusterID"
    colour="sampleID"
    if(sumAll){
      colour="clusterID"
      linetype="sampleID"
    }
    ylab="depth"
    if(!is.null(total_reads)) ylab="depth per million mapped reads"
    #invisible(tpm_df)
   if(showWaterfall && !is.null(seq_df)){
    allpos =  seq_df$pos
  #  print(seq_df$sequence[xlim[1]:xlim[2]])
    kstart = (waterfallKmer-1)/2
    kmers = getKmer(as.character(seq_df$sequence),1:length(seq_df$sequence),v=-kstart:kstart + waterfallOffset)
  #  print(head(kmers))
    levsk = levels(factor(kmers))
   # print(levsk)
    cnts = unlist(lapply(levsk, function(x)sum(tpm_df$count[tpm_df$pos %in% allpos[which(kmers==x)]])  ))
  #  cnts = cnts/sum(cnts)
    print(cnts)
   cnt_df =  data.frame(cnts=cnts, kmers=levsk)
   top10 = min(top10,length(cnts))
  # print(top10)
    cnt_df = cnt_df[order(cnt_df$cnts,decreasing = T)[1:top10],,drop=F]
    #print(cnt_df)
    id = 1:dim(cnt_df)[1]
    end = cumsum(cnt_df$cnts)
    start = c(0,end[1:(length(end)-1)])
    cnt_df = cbind(id,start,end,cnt_df)
    #print(cnt_df)
    cnt_df$kmers = factor(as.character(cnt_df$kmers),levels = as.character(cnt_df$kmers))
   return(cnt_df)
    
  #  ggp<-ggp+theme(text = element_text(size=10), axis.text.x = element_text(size = rel(0.7), angle = 25, hjust=0.75))
    
   }else{
   return(tpm_df)
    
   }
    return(NULL)
    }
  
	readDir <- function() {
    print(input$dir)
	  h5closeAll()
    print(" updating input dir")
    session$userData$dataDepth = list("depth"=data.frame(), "depthStart"=data.frame(), "depthEnd"=data.frame())
    session$userData$dataPlot = list("depth"=NULL, "depthStart"=NULL, "depthEnd"=NULL)
    session$userData$depthOpts = list("depth"=list(), "depthStart"=list(), "depthEnd"=list())
    
    currdir = paste(basedir,input$dir,sep="/")
    datafile = paste(currdir,"0.isoforms.h5",sep="/")
    h5file=paste(currdir,"0.clusters.h5",sep="/")   
    if(file.exists(datafile)){
      isoInfo = .getIsoInfo(datafile, h5file,toreplace)
      total_reads = isoInfo$total_reads
      ##this gets order by counts
      if(reorder){
        order = .readTotalIso(datafile, group="/trans", trans=as.character(isoInfo$orfs$ORFs))
        isoInfo$orfs = isoInfo$orfs[order,,drop=F]
      }
      
      info=.processInfo(isoInfo)
      print(paste("set", datafile))
      ch=c(names(info$choices1), names(info$choices))
      type_nmes=names(total_reads)
      session$userData$isoInfo=isoInfo
      session$userData$info=info
      session$userData$total_reads = total_reads
      session$userData$header=type_nmes
    }else{
      ch = c()
    }
   annot_file = paste(currdir, "annotation.csv.gz",sep="/")
    if(file.exists(annot_file)){
      annots = read.table(annot_file,sep="\t", head=F)
      names(annots) = c("chr","ens","name","ens1","type")
      session$userData$annots=annots
    }
   counts_file = paste(currdir, "Counts_genome1.csv",sep="/")
   options2 = c("showTranscriptPlot","logy","showCI", "TPM_amongst_viral","barchart","ribbonCI","mergeCounts", "stacked", "reverseOrder")
   totick2 = c("showTranscriptPlot","ribbonCI","barchart","TPM_amongst_viral")
   
   if(file.exists(counts_file)){
     countsHostVirus= read.csv(counts_file)
     options2 = c("showTranscriptPlot","logy","showCI", "TPM_amongst_all" ,"TPM_amongst_viral","barchart","ribbonCI","mergeCounts", "stacked", "reverseOrder")
     
     names(countsHostVirus) =  gsub("X.","",names(countsHostVirus))
    
     sample= apply(countsHostVirus[,1:3],1,paste,collapse="_")
     countsHostVirus=cbind(sample,countsHostVirus)
     inds_a = grep("Map.to", names(countsHostVirus))
     sum_inds = names(countsHostVirus) %in% c("Host","Virus","Sequin","Total")
     torem = c()
     v1 = sample
     dupl = lapply(unique(v1[duplicated(v1)]), function(x)which(v1 ==x))
     for(j in dupl){
      # print(countsHostVirus[j,sum_inds,drop=F])
       summed = apply(countsHostVirus[j,sum_inds,drop=F],2,sum)
      # print(j)
      # print(summed)
      # print(countsHostVirus[j[1],sum_inds])
         countsHostVirus[j[1], sum_inds] = summed
       torem = c(torem, j[-1])
     }   
     countsHostVirus = countsHostVirus[-torem,]
     vars = c(names(countsHostVirus[inds_a]), "Total")
     countsHostVirus=
      melt(countsHostVirus,id.vars=c("sample"),
                   measure.vars=which(names(countsHostVirus) %in% vars), variable.name="ID", value.name="count") %>%
       transform(sample=factor(sample), ID=factor(ID))
     session$userData$countsTotal = countsHostVirus[countsHostVirus$ID=="Total",]
     
     session$userData$countsHostVirus = countsHostVirus[countsHostVirus$ID!="Total",]
     
   }
     
    coords_file = paste(currdir, "Coordinates.csv",sep="/")
    motifText = ""
    if(file.exists(coords_file)){
      t = readCoords(coords_file)
      session$userData$t=t
      orfs=paste(t$gene,collapse=",")
   
      fimo_file = paste(currdir,"fimo.tsv",sep="/")
      fimo=read.table(fimo_file, sep="\t", head=T)
      print(names(fimo))
      motifText = paste(levels(factor(fimo$matched_sequence)),collapse="|")
      session$userData$fimo = fimo
    }else{
      orfs = c()
     
    }
    fastafile = grep("extra",grep("fasta.gz",dir(currdir),v=T),v=T,inv=T)
    if(length(fastafile)>=1){
      fastaseq = read.fasta(paste(currdir,fastafile[1],sep="/"))
      session$userData$fasta =fastaseq
      session$userData$fastaseq = lapply(fastaseq,paste,collapse="")
    }else{
      session$userData$fastaseq=NULL
      session$userData$fasta=NULL
    }
    peptide_file =paste(currdir,"peptides.csv",sep="/")
    if(file.exists(peptide_file)){
      
      peptide=read.csv(peptide_file,  head=F, comment.char="#")
      trans_vals= as.numeric(sub("#","",read.csv(peptide_file, head=F,nrow=1)))
      peptide = peptide[!duplicated(peptide[,2]),-1,drop=F]
      
      peptide= apply(peptide,c(1,2),function(x) (x-1)*trans_vals[2]+trans_vals[1])
      
      
     # peptide=read.csv(peptide_file,  head=F)
    #  peptide = peptide[!duplicated(peptide[,2]),-1,drop=F]
      names(peptide)=c("start","end")
      session$userData$peptide = peptide
      
    }
    session$userData$dirname = gsub("/","_",input$dirname)
    session$userData$currdir=currdir
    session$userData$datafile=datafile
    session$userData$h5file=h5file
    session$userData$results = list()
    session$userData$results = list()
	
	
	print(dir.exists(file.path(currdir,'DE')))
	
	updateCheckboxInput(session, "LoadDE", value = FALSE)
	hide("plotDE")
	hide("DE_cell1")
	hide("DE_cell2")
	hide("DE_time1")
	hide("DE_time2")	
	
	
	print('setting DE_countdata to NULL')
	DE$counts = NULL
	
	if (!dir.exists(file.path(currdir,'DE'))) {
		shinyjs::disable('LoadDE') 
		} else {
		session$userData$DE = file.path(currdir, 'DE')
		shinyjs::enable('LoadDE')
	}
	
	#toggleState('ActivateDE', condition = dir.exists(file.path(currdir,'DE')))
    updateSelectInput(session,"plottype", label = "Category 1", choices=ch, selected=input$plottype)
   # updateSelectInput(session,"plottype1", label = "Category 2", choices=ch, selected=input$plottype1)
    updateSelectInput(session, "toplot5",label = paste("Transcript",names(info$choices1)[1]),choices=c("-",info$choices1[[1]]),selected='-')
  #  updateSelectInput(session, "toplot6",label = paste("Transcript",names(info$choices1)[2]),choices=c("-",info$choices1[[2]]),selected=input$toplot6)
    updateCheckboxGroupInput(session,"molecules", label = "Molecule type",  choices =info$molecules, selected = info$molecules)
    updateCheckboxGroupInput(session,"cells", label = "Cell type",  choices = info$cells, selected = info$cells)
    updateCheckboxGroupInput(session,"times", label = "Time points",  choices = info$times, selected = info$times)
    updateTextInput(session,"orfs", label="ORFs to include", value = orfs)
   # updateSelectInput(session, "depth_plot_type", label ="What to plot", choices=plot_type_ch, selected="depth")
    updateTextInput(session,"motif", label="Show motif", value = motifText)
    
    updateCheckboxGroupInput(session,"options2", label = "Middle panel", choices = options2, selected=totick2) 
    
    #CTAAAC|TTAAAC
    #ACGAAC|ACGATC|ATGAAC
	
	
  }
  
  loadDE <- function(thresh) {
  
  # Import data 
  message(paste('importing DE Data'))
    currdir = session$userData$currdir
  count_files <- list.files(path = file.path(currdir, 'DE'), recursive = T, full.names = T)
  cell_types <- lapply(strsplit(x = count_files, split = '/'), rev) %>% sapply('[', 2) 
  print(cell_types)
  
  print(count_files)
  
  countdata <- lapply(count_files, FUN = 
                        function(x) {read.table((gzfile(x)), header = T, row.names = 1) %>%
                              subset(select = -c(1:5)) %>% 
                                as.matrix() } )
  names(countdata) <- cell_types
  
 #countdata=lapply(countdata,
  #                function(x) .subsetFCFile(x,toplot5,toplot7,toplot8,tojoin,group_by)) 
 
 #print(head(countdata[[1]]))
  return(countdata)
}
  
  
  
	extract_sequence=function(){
	  fasta=session$userData$fasta[[1]]
	  xlim =   c(isolate(input$min_x), isolate(input$max_x))
	  fasta[xlim[1]:xlim[2]]
	}
	extract_sequence_name=function(){
	  fasta=names(session$userData$fasta)[[1]]
	  xlim =   c(isolate(input$min_x), isolate(input$max_x))
	  paste(fasta,xlim[1],xlim[2],sep=" ")
	}
	

	
  depthPlot= function(plot_type, reuse=F) {
    #result = loadData();
 reuseData=F
   if(!file.exists(session$userData$h5file)) return(ggplot())
 
   
    showDepth  = "show_depth" %in% input$options3
    logy = "logy" %in% input$options3
    textsize=input$textsize
    ci=0
    if(  "showCI" %in% input$options3){
      ci=input$conf.int
    }
    
    depth_thresh = input$depth_thresh
    group_by=input$group_by
    reverseOrder=F
    merge_by="" #input$merge_by
  #  plot_type=input$depth_plot_type
    motif = isolate(input$motif)
    fastaseq =session$userData$fastaseq[[1]]
    fasta=session$userData$fasta[[1]]
    maxKmers=isolate(input$maxKmers)
    showORFs="showORFs" %in% input$options3
    fisher = input$test =="fisher"
    showErrors = "showErrors" %in% input$options3
    plotCorr = "plotCorr" %in% input$options3
    showSequence="showSequence" %in% input$options3
    if(plot_type=="depth") showWaterfall=F else showWaterfall="showWaterfall" %in% input$options3
    waterfallKmer=isolate(input$waterfallKmer)
   waterfallOffset=isolate(input$waterfallOffset)
    
    motifpos=list()
    if(nchar(motif>0) && !is.null(fastaseq)){
      motifpos= lapply(strsplit(motif,"\\|")[[1]], function(x) gregexpr(x,fastaseq, ignore.case=T)[[1]])
    }
    mergeCounts='mergeCounts' %in% input$options3
    showPeptides="showPeptides" %in% input$options3
    downsample="downsample" %in% input$options3
    tpm = "TPM_amongst_viral" %in% input$options3
    h5file=session$userData$h5file
    total_reads = NULL
   # fimo = NULL
    t = NULL
    peptides=NULL
   # print(showMotifs)
    #if(showMotifs){
    #fimo = session$userData$fimo
    #}
    #print(fimo)
    if(showORFs){
    t = session$userData$t
    }
    if(showPeptides){
      peptides=session$userData$peptide
    }
    if(tpm){
      total_reads = session$userData$total_reads
    }
    #  print(total_reads)
    # print(h5file)
    ggp=ggplot()
    if(length(grep(plot_type,h5ls(h5file)$group))<=0) return (ggp)
  #  span = 0
    span=input$loess
    if(showDepth  && !is.null(h5file)){
      if(file.exists(h5file)){
        toplot = c(isolate(input$toplot5))#,isolate(input$toplot6))#,isolate(input$toplot7),isolate(input$toplot8))
        tojoin=isolate(input$tojoin)
        merge=F
        toplot = toplot[toplot!="-"]
        combinedID="combined"
        sumAll = length(toplot)>1
        if(length(toplot)==0){
          ##need to get list from this
          toplot=c(isolate(input$toplot7),isolate(input$toplot8))
          toplot = toplot[unlist(lapply(toplot,nchar))>2]
          if(length(toplot)>0){
            combinedID=toplot[1]
            toplot = .findEntries(toplot,h5file,"/depth",tojoin)
          }
          
          merge=T
          sumAll=F
        }
        if(length(toplot)>0 && nchar(toplot[1])>2 ){
          if("sumDepth" %in% input$options3) sumAll=TRUE
          mergeGroups=NULL
          if(mergeCounts) group_by="all"
          if(nchar(merge_by)>0){
            mergeGroups=.mergeGroups(toplot,merge_by)
          
          }else if(group_by != 'No grouping'){
            mergeGroups = .getGroups(toplot,group_by)
          }else if(length(toplot)>input$maxtrans){
            ##only show top number if not merging
            isoInfo=session$userData$isoInfo
            inds_k = sort(match(toplot,isoInfo$orfs$ORFs))
          #  print(inds_k)
         #   print(isoInfo$orfs$ORFs[inds_k][1:input$maxtrans])
            toplot = isoInfo$orfs$ORFs[inds_k][1:input$maxtrans]
          }
          molecules=input$molecules
          cells=input$cells
          times = input$times
        xlim =   c(isolate(input$min_x), isolate(input$max_x))
        alpha = isolate(input$alpha)
        linesize=0.001  #isolate(input$linesize)
        if(xlim[2]<=xlim[1]) xlim = NULL
       seq_df= NULL
        if( showSequence || showWaterfall){
          if(xlim[1]<1) xlim[1] = 1
          if(xlim[2]>length(fasta)) xlim[2] = length(fasta)
          pos = xlim[1]:xlim[2]
          sequence =  fasta[pos]
          seqy = rep(1,length(pos))
          seq_df = data.frame(pos,sequence, seqy) %>%
            transform(pos=as.numeric(pos), sequence=factor(sequence, levels=c("a","c","t","g")))
        }
       plotOpts = list(showORFs = showORFs, motifpos=motifpos,peptides=peptides,xlim =xlim, t=t,
                       alpha=alpha,linesize=linesize, textsize=textsize,calcErrors=showErrors,logy=logy )
       depthOpts = 
         list(total_reads=total_reads,toplot=toplot, downsample = downsample, span = span, 
              mergeGroups=mergeGroups,molecules=molecules, combinedID=combinedID, cells=cells, 
              times = times, sumAll = sumAll,
              plotCorr=plotCorr,reverseOrder=reverseOrder,
              fisher=fisher,
              ci = ci, depth_thresh = depth_thresh,
              showWaterfall=showWaterfall,waterfallKmer=waterfallKmer,waterfallOffset=waterfallOffset, top10=maxKmers)
       depthOpts_o = session$userData$depthOpts[[which(names(session$userData$depthOpts)==plot_type)]]
     #  plotOpts_o = session$userData$plotOpts[[which(names(session$userData$plotOpts)==plot_type)]]
       
       if(reuse ){
         ggp = session$userData$dataPlot[[which(names(session$userData$dataDepth)==plot_type)]]
       }else{
          
            # path=plot_type,seq_df = seq_df
            
           
            if(length(depthOpts_o)>0){
            reuseData = identical(depthOpts, depthOpts_o )
            }else{
              reuseData = F 
            }
            session$userData$depthOpts[[which(names(session$userData$depthOpts)==plot_type)]] = depthOpts
           
            if(reuseData){
              print("reusing data")
              tpm_df = session$userData$dataDepth[[which(names(session$userData$dataDepth)==plot_type)]]
            }else{
              print("not reusing data")
             # print("mergeGroups")
            #  print(mergeGroups)
          tpm_df=run_depth(h5file,total_reads,toplot, seq_df=seq_df, downsample = downsample, span = span, 
                           mergeGroups=mergeGroups,molecules=molecules, combinedID=combinedID, cells=cells, 
                           times = times,logy=logy, sumAll = sumAll,
                    showORFs = showORFs, motifpos=motifpos,peptides=peptides,xlim =xlim, t=t,path=plot_type,
                     alpha=alpha,plotCorr=plotCorr,linesize=linesize, reverseOrder=reverseOrder,
                    textsize=textsize, calcErrors=showErrors,fisher=fisher,
                    ci = ci, depth_thresh = depth_thresh,
                    showWaterfall=showWaterfall,waterfallKmer=waterfallKmer,waterfallOffset=waterfallOffset, top10=maxKmers
                    )
            # print(head(tpm_df))
          session$userData$dataDepth[[which(names(session$userData$dataDepth)==plot_type)]] = tpm_df
          }
         ggp =plot_depth(tpm_df,total_reads,toplot, seq_df=seq_df, downsample = downsample, span = span, mergeGroups=mergeGroups,molecules=molecules, combinedID=combinedID, cells=cells, times = times,logy=logy, sumAll = sumAll,
                     showORFs = showORFs, motifpos=motifpos,peptides=peptides,xlim =xlim, t=t,path=plot_type,
                      alpha=alpha,plotCorr=plotCorr,linesize=linesize, reverseOrder=reverseOrder,
                     textsize=textsize, calcErrors=showErrors,fisher=fisher,
                     ci = ci, depth_thresh = depth_thresh,
                     showWaterfall=showWaterfall,waterfallKmer=waterfallKmer,waterfallOffset=waterfallOffset, top10=maxKmers
          )
          session$userData$dataPlot[[which(names(session$userData$dataDepth)==plot_type)]] =ggp
       }
          return(ggp)
        }
        #run_depth(h5file,toplot=c("leader_leader,N_end")) 
      }
    }
  }
    
  
  
  transcriptPlot=function(){
    if(!file.exists(session$userData$datafile)) return(ggplot())
    if(!"showTranscriptPlot" %in% input$options2) return(ggplot())

	p_data= list(molecules=c("RNA","cDNA"),cells=c("calu","vero"),times=c("24hpi","48hpi"),toplot="all",splitby="NA",xy=FALSE,
	             showTPM=T,merge=F,barchart=T,reverseOrder=T,stack=T,calcTPMFromAll=T,group_by="type",
	             tojoin="OR",usergrep=T,merge_by="",method="logit",conf.int=0.95)
	p_plot = list(textsize=20, logy=T, showCI=F,riboon=F,angle=25)
	
	p_plot$textsize=input$textsize
  p_data$molecules <-  input$molecules 
  p_data$cells <- input$cells 
  p_data$times<-input$times
  p_data$splitby=input$splitby
  splitby_vec=NULL
  p_data$xy=FALSE
    if( p_data$splitby=="molecules"){
      p_data$xy=T
      p_data$splitby_vec=molecules
    }else if( p_data$splitby=="cells"){
      p_data$xy=T
      p_data$splitby_vec=cells
    }else if( p_data$splitby=="times"){
      p_data$xy=T
      p_data$splitby_vec=times
    }
    
    p_plot$facet=input$facet
  
    p_plot$logy = "logy" %in% input$options2
    p_plot$showCI = "showCI" %in% input$options2
    p_plot$ribbon="ribbonCI" %in% input$options2
    p_data$showTPM="TPM_amongst_viral" %in% input$options2 || "TPM_amongst_all" %in% input$options2
    p_plot$showCI = "showCI" %in% input$options2
    p_data$merge='mergeCounts' %in% input$options2
    p_data$barchart="barchart" %in% input$options2
    p_data$reverseOrder="reverseOrder" %in% input$options2
    p_data$stack = "stacked" %in% input$options2
    p_data$calcTPMFromAll = "TPM_amongst_all"  %in% input$options2
    p_data$group_by=input$group_by
    p_plot$angle = input$angle
    p_data$merge_by=""  #input$merge_by
    p_data$max_trans = input$maxtrans
    p_data$conf.int=input$conf.int
    p_data$method="logit";
    p_data$toplot = c(isolate(input$toplot5))#,isolate(input$toplot6))#,isolate(input$toplot7),isolate(input$toplot8))
    p_data$toplot = p_data$toplot[p_data$toplot!="-"]
    p_data$usegrep=F
    if(length(p_data$toplot)==0){
      p_data$toplot=c(isolate(input$toplot7),isolate(input$toplot8))
      p_data$toplot = p_data$toplot[unlist(lapply(p_data$toplot,nchar))>2]
      p_data$usegrep=T
    }
    p_data$tojoin=isolate(input$tojoin)

    if(length(p_data$toplot)==0 || is.null(session$userData$datafile ) ){
      return(ggplot())
    }
    datafile = session$userData$datafile ;
    total_reads = session$userData$total_reads
    countsTotal = session$userData$countsTotal
    countsHostVirus= session$userData$countsHostVirus
    prev_params = session$userData$prev_params
    reuseData=F
    if(!is.null(prev_params)){
      
      reuseData = identical(p_data, prev_params$p_data )
      reusePlot = reuseData && identical(p_plot, prev_params$p_data )
      if(reusePlot){
        print("reusing tpm plot")
        ggp=session$userData$tpm_plot
        if(!is.null(ggp)) return(ggp)
      }
    }
    session$userData$prev_params = list(p_plot=p_plot, p_data =p_data)
    if(reuseData && !is.null(session$userData$results)){
      print("reusing tpm data")
      results_ = session$userData$results
      subs = results_$data
      countsHostVirus1 = results_$totals
    }else{
      subs = .extractTPM(datafile ,  total_reads,countsTotal, p_data)
      countsHostVirus1 = NULL
      if(!is.null(countsHostVirus)){
       # countsHostVirus = countsHostVirus[countsHostVirus$ID!="Total",]
        countsHostVirus1 = countsHostVirus[which(countsHostVirus$sample %in% subs$sample),,drop=F]
        names(countsHostVirus1)[2]="Type"
        names(countsHostVirus1)[3]="Reads"
        cell = unlist(lapply(as.character(countsHostVirus1$sample), function(x) strsplit(x,"_")[[1]][1]))
        molecule_type = unlist(lapply(as.character(countsHostVirus1$sample), function(x) strsplit(x,"_")[[1]][2]))
        time = unlist(lapply(as.character(countsHostVirus1$sample), function(x) strsplit(x,"_")[[1]][3]))
        time =  factor(time,level= paste(sort(as.numeric(unique(sub("hpi","",time)))),"hpi",sep=""))
        countsHostVirus1 = cbind(countsHostVirus1, cell, molecule_type, time)
      }
    }
    session$userData$results = list(data=subs, totals=countsHostVirus1);
    yname='TPM'
    if(!p_data$showTPM){
          yname="Counts"
    }
      
    ggp=.plotTPMData(subs,countsHostVirus1, p_data,p_plot,yname)
    session$userData$tpm_plot = ggp
    ggp
    
  }
  
  infectivityPlot=function(){
    if(!"showInfectivity" %in% input$options1) return(ggplot())
    currdir = session$userData$currdir
    barchart="barchart" %in% input$options1
    reverseOrder="reverseOrder" %in% input$options1
    showSecondAxis="showSecondAxis" %in% input$options1
    textsize=input$textsize
    conf.int=input$conf.int
    countsTotal=session$userData$countsTotal
    infilesAnnot = paste(currdir,"0.annot.txt.gz", sep="/")
    total_reads = session$userData$total_reads
    type_nme= names(total_reads)
    total_reads1 = total_reads
    if(!is.null(countsTotal)){
      total_reads1=countsTotal$count[match(type_nme,countsTotal$sample)]
      names(total_reads1) = type_nme
    }
    if(file.exists(infilesAnnot)){
      molecules=input$molecules; cells=input$cells; times = input$times;
      
      
      levels1 = .getlevels(type_nme,molecules, cells, times,reverseOrder)
      orfs=input$orfs
      norm=F
      ratio1 = .readAnnotFile(infilesAnnot,total_reads1,norm=norm,levels= levels1, conf.level=conf.int, orfs=orfs)
      #.plotAnnotFile<-function(ratio1,molecule_type, cell, time, levels=NULL,barchart=F,showEB = F,showSecondAxis=F){
      y_text="Ratio"
      #   y_text="spliced"
      annots1=.plotAnnotFile(ratio1,barchart=barchart,showSecondAxis=showSecondAxis,showEB=T, levels=levels1, y_text=y_text,
                             diff=0, coeff=5, textsize=textsize)
     
#      resall = 
     # names(resall) =   session$userData$dirname
      session$userData$resultsInf = list(data=annots1$data)
      annots1$ggp
    }else{
      ggplot()
    }
  }
  

	
# DE_countdata <- reactive( {
	
	# DE_countdata = loadDE()
	# print(names(DE_countdata))
	# updateSelectInput(session, "DE_cell",  choices = names(DE_countdata))
	# updateSelectInput(session, "DE_time1",  choices = info$times)
	# updateSelectInput(session, "DE_time2",  choices = info$times)
	# return(DE_coundata)
		
	# })
	
		

  
  observeEvent(input$plottype,{
    dd= .process1(input$plottype,  session$userData$info)
    updateSelectInput(session,"toplot5", label = dd$label, choices=dd$ch, selected="-")
  })
  observeEvent(input$dir, readDir() )		  
	#THIS JUST EXAMPLE FOR RENDERING A PLOT
	#REACTIVE ON PLOT BUTTON



	output$infPlot<-renderPlot({
	  input$plotButton
	   try(validate(need(input$dir, 'Please select a directory to begin')))
	  infectivityPlot()
	})

	output$distPlot <- renderPlot({
		try(validate(need(input$dir, '')))
		#validate(need(input$toplot5 != "-", 'Select a transcript to plot'))
		try(validate(need(length(input$molecules) > 0 & length(input$cells) > 0 & length(input$times) > 0, 'At least one molecule, cell and time point must be supplied') ))
	    input$plotButton
  	    transcriptPlot()
  	  })

	output$depthPlot <- renderPlot({
	    input$plotButton
	  depthPlot("depth")
	})
	output$depthStartPlot <- renderPlot({
	  input$plotButton
	  depthPlot("depthStart")
	})
	output$depthEndPlot <- renderPlot({
	  input$plotButton
	  depthPlot("depthEnd")
		
	})


DE <- reactiveValues(counts = NULL, main_out = NULL)


observeEvent( input$LoadDE,	
{
  if(input$LoadDE){
	DE$counts = loadDE()
  }else{
    print("disactivating DE")
    DE$counts = list()
  }
	print(paste('loadDE output'))
	updateSelectInput(session, "DE_cell1",  choices = names(DE$counts))
	updateSelectInput(session, "DE_cell2", choices = names(DE$counts))
	updateSelectInput(session, "DE_time1",  choices = session$userData$info$times)
	updateSelectInput(session, "DE_time2",  choices = session$userData$info$times)



toggle("plotDE")
toggle("DE_cell1")
toggle("DE_cell2")
toggle("DE_time1")
toggle("DE_time2")

} )


#DE Plots
 observeEvent(input$plotDE, {
 head(DE$counts)
 if (is.null(DE$counts)) {print('DE_countdata is null') } else {
	#head(DE$counts)
   plot_params = list(toplot5=input$toplot5, toplot2=c(input$toplot7, input$toplot8), 
                      tojoin=input$tojoin, group_by=input$group_by, merge_by=input$merge_by)
	DE$main_out <- try(runDE(count_list = DE$counts, cell1 = input$DE_cell1 ,
	                     cell2 = input$DE_cell2, time1 = input$DE_time1, 
	                     time2 = input$DE_time2,  thresh=input$mean_count_thresh,
	                     plot_params=plot_params))
	if(!inherits(DE$main_out,"try-error")) {
	print(paste('DEout done', names(DE$main_out)))
	
	
	output$DEPlot_volcano <- renderPlot( {
	DE$main_out[['volcano_params']][['remove_spurious']] <- input$remove_spurious 
	do.call(volcanoplot, 
	DE$main_out[['volcano_params']] )
	})
	
	output$DEPlot_PCA <- renderPlot( {
	do.call(rld_pca, 
	DE$main_out[['rld_pca_params']]
	)
	})
	}
	
		}
	 } )
	 
#Downloads
output$downloadInf <- downloadHandler(filename = function() {'plotInfectivity.pdf'}, content = function(file) ggsave(file, infectivityPlot(), device='pdf', height = 20, width = 40, units='cm' ) )
output$downloadDepth <- downloadHandler(filename = function() {'plotDepth.pdf'}, content = function(file) ggsave(file, depthPlot("depth", T), device = 'pdf', height = 20, width = 40, units='cm') )
output$downloadDepthStart <- downloadHandler(filename = function() {'plotDepthStart.pdf'}, content = function(file) ggsave(file, depthPlot("depthStart", T), device = 'pdf', height = 20, width = 40, units='cm') )
output$downloadDepthEnd <- downloadHandler(filename = function() {'plotDepthEnd.pdf'}, content = function(file) ggsave(file, depthPlot("depthEnd",T), device = 'pdf', height = 20, width = 40, units='cm') )
output$downloadDist <- downloadHandler(filename = function() {'plotDist.pdf'}, content = function(file) ggsave(file, transcriptPlot(), device = 'pdf' , height = 20, width = 40, units='cm') )
output$downloadResults<-downloadHandler(filename = function() {'results.xlsx'}, content = function(file) write_xlsx( session$userData$results,file ) )
output$downloadResultsInf<-downloadHandler(filename = function() {'resultsInf.xlsx'}, content = function(file) write_xlsx( session$userData$resultsInf,file ) )
output$downloadResultsDepth<-downloadHandler(filename = function() {'resultsDepth.xlsx'}, content = function(file) write_xlsx( session$userData$dataDepth,file ) )
output$downloadSequence<-downloadHandler(filename = function() {'sequence.fa'}, content = function(file) write.fasta(extract_sequence(),extract_sequence_name(), file ) )



output$downloadPCA <- downloadHandler(filename = function() {'plotPCA.pdf'}, content = function(file) {
	pdf(file, 18, 18, pointsize=50)
	do.call(rld_pca, DE$main_out[['rld_pca_params']])
	dev.off()
		})
output$downloadVOLCANO <- downloadHandler(filename = function() {'plotVOLCANO.pdf'}, content = function(file) {
	pdf(file, 18, 18, pointsize=20)
	do.call(volcanoplot, DE$main_out[['volcano_params']] )
	dev.off()
		})
output$downloadDEdata <- downloadHandler(filename = function() {'DE_data.xlsx'}, content = function(file) {
	write_xlsx(DE$main_out[['data']], file)
	})
})
