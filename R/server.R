library(shiny)
library(reshape2)
library(ggplot2)
library(tidyr)
library(rhdf5)
library(RColorBrewer)
library(binom)
library(writexl)


source( "transcript_functions.R")

basedir="../data"
#toreplace=list(virion="RNA_virion_0hpi", whole_genome_mapped="RNA_vero_24hpi")
decodeFile = paste(basedir,"decode.txt",sep='/')
replace=read.table(decodeFile,sep="\t",head=F)
toreplace = replace[,2]
names(toreplace) = replace[,1]



reorder=T

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
    names(l) = c("5_3", "5_no3","no5_3","no5_no3") 
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
  l[[2]] = grep(group_by, x1,inv=T,v=t )
  names(l) = c(group_by,paste("!",group_by))
  }
  l[ unlist(lapply(l, length))>0]
}
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
.getlevels<-function(type_nme, molecules, cells, times){
  types_=data.frame(t(data.frame(strsplit(type_nme,"_"))))
  names(types_) = c("molecules","cell","time")
  inds1 =  which(types_$molecules %in% molecules & types_$cell %in% cells & types_$time %in% times)
  types1_ = types_[inds1,,drop=F]
  ord = order(as.numeric(factor(types1_$time, levels=c("0hpi", "2hpi","24hpi","48hpi"))),types1_$cell,types1_$molecules)
  levels1=type_nme[inds1][ord]
  attr(levels1,"inds1") = inds1
  levels1
}
#run_depth(h5file, total_reads=total_reads)
run_depth<-function(h5file, total_reads=NULL,  toplot=c("leader_leader,N_end", "N_end"),combinedID="combined", gapthresh=100, mergeGroups=NULL,molecules="RNA",cells="vero",times=c('2hpi','24hpi','48hpi'), 
                    span = 0.01, sumAll=F, xlim=null, fimo=NULL, alpha=1.0,t= NULL,logy=T, showMotifs=F,showORFs = F, path="depth"){

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
 #print(inds1)
 id_cols = c("molecule","cell","time")
 tot_reads=NULL
 if(!is.null(total_reads)){
   
   tot_reads =  total_reads[inds1]/rep(1e6,length(inds1))
 }
 #print(tot_reads)
   	clusters_ = readH5(h5file,tot_reads, c("pos",header[inds1+1]),toAdd = toAdd, mergeGroups=mergeGroups,sumID=sumID, path=path,toplot,id_cols=id_cols, gapthresh=gapthresh, dinds = dinds[inds1], pos =NULL, span = span, cumul=F, sumAll=sumAll)
#print(clusters_)

  # if(!is.null(xlim)){
  #   minx = min(clusters_$pos)
  #   maxx = max(clusters_$pos)
  #   xlim[1] = max(minx,xlim[1])
  #   xlim[2] = min(maxx, xlim[2])
  # }
   	if(is.null(clusters_)){
  print(paste("could not read ",toplot))
 return (ggplot())
   	}
   	
   	if(sumAll) levs=names(clusters_)[3]
   	
   	tpm_df = melt(clusters_,id.vars=c("clusterID","pos"), measure.vars=names(clusters_)[-(1:2)], variable.name="sampleID",value.name='count') %>%
   	transform(sampleID=factor(sampleID,levels=levs))
   	  #  separate(variable, c('molecule_type', 'cell', 'time'), sep='_', remove = T)  %>%
   	 # transform( count=as.numeric(count), molecule_type = factor(molecule_type), cell = factor(cell), time = factor(time, levels = time_vec)) 
   	
   		if(sumAll) type_nme = "combined"
	#mat=t(h5read(h5file,paste("depth", toplot[1],sep="/")))
rawdepth = T
leg_size1=10
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
 plotClusters(tpm_df, 4,  1, 
              t,
             fimo,
               rawdepth = rawdepth, linetype=linetype, colour=colour, alpha=alpha, xlim = xlim,ylab=ylab , title =path, logy=logy, leg_size =leg_size1, show=show, fill =fill)
}

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
  cis = matrix(0,ncol=3,nrow =nrows )
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
	readDir <- function() {
    print(input$dir)
    print(" updating input dir")
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
    coords_file = paste(currdir, "Coordinates.csv",sep="/")
    if(file.exists(coords_file)){
      t = readCoords(coords_file)
      session$userData$t=t
      orfs=paste(t$gene,collapse=",")
   
      fimo_file = paste(currdir,"fimo.tsv",sep="/")
      fimo=read.table(fimo_file, sep="\t", head=T)
      session$userData$fimo = fimo
    }else{
      orfs = c()
    }
    session$userData$dirname = gsub("/","_",input$dirname)
    session$userData$currdir=currdir
    session$userData$datafile=datafile
    session$userData$h5file=h5file
    session$userData$results = list()

    updateSelectInput(session,"plottype", label = "Category 1", choices=ch, selected=input$plottype)
   # updateSelectInput(session,"plottype1", label = "Category 2", choices=ch, selected=input$plottype1)
    updateSelectInput(session, "toplot5",label = paste("Transcript",names(info$choices1)[1]),choices=c("-",info$choices1[[1]]),selected='-')
  #  updateSelectInput(session, "toplot6",label = paste("Transcript",names(info$choices1)[2]),choices=c("-",info$choices1[[2]]),selected=input$toplot6)
    updateCheckboxGroupInput(session,"molecules", label = "Molecule type",  choices =info$molecules, selected = info$molecules)
    updateCheckboxGroupInput(session,"cells", label = "Cell type",  choices = info$cells, selected = info$cells)
    updateCheckboxGroupInput(session,"times", label = "Time points",  choices = info$times, selected = info$times)
    updateTextInput(session,"orfs", label="ORFs to include", value = orfs)
   # updateSelectInput(session, "depth_plot_type", label ="What to plot", choices=plot_type_ch, selected="depth")
    
  }
	
	
  depthPlot= function(plot_type) {
    #result = loadData();
    if(!file.exists(session$userData$h5file)) return(ggplot())
    showDepth  = "show_depth" %in% input$options3
    logy = "logy" %in% input$options3
    group_by=input$group_by
  #  plot_type=input$depth_plot_type
    
    showORFs="showORFs" %in% input$options3
    showMotifs="showMotifs" %in% input$options3
    mergeCounts='mergeCounts' %in% input$options3
    
    tpm = "TPM" %in% input$options3
    h5file=session$userData$h5file
    total_reads = NULL
    fimo = NULL
    t = NULL
   # print(showMotifs)
    if(showMotifs){
    fimo = session$userData$fimo
    }
    #print(fimo)
    if(showORFs){
    t = session$userData$t
    }
    if(tpm){
      total_reads = session$userData$total_reads
    }
    #  print(total_reads)
    # print(h5file)
    ggp=ggplot()
    if(length(grep(plot_type,h5ls(h5file)$group))<=0) return (ggp)
    span = 0
    if(plot_type=="depth") span=input$loess
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
          if(nchar(group_by)>0){
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
        #print("xlim")
        #print(xlim)
        #xlim= NULL
        if(xlim[2]<=xlim[1]) xlim = NULL
          ggplot=run_depth(h5file,total_reads,toplot, span = span, mergeGroups=mergeGroups,molecules=molecules, combinedID=combinedID, cells=cells, times = times,logy=logy, sumAll = sumAll,
                    showORFs = showORFs, fimo=fimo,xlim =xlim, t=t,path=plot_type, showMotifs =showMotifs, alpha=alpha) 
        }
        #run_depth(h5file,toplot=c("leader_leader,N_end")) 
      }
    }
  }
    
  
  
  transcriptPlot=function(){
    if(!file.exists(session$userData$datafile)) return(ggplot())
	print(paste('testinput', input$molecules))
    molecules <-  input$molecules 
    cells <- input$cells 
    times<-input$times
    splitby=input$splitby
    splitby_vec=NULL
    xy=FALSE
    if(splitby=="molecules"){
      xy=T
      splitby_vec=molecules
    }else if(splitby=="cells"){
      xy=T
      splitby_vec=cells
    }else if(splitby=="times"){
      xy=T
      splitby_vec=times
    }
    
    
    logy = "logy" %in% input$options2
    showCI = "showCI" %in% input$options2
    ribbon="ribbonCI" %in% input$options2
    showTPM="TPM" %in% input$options2
    showCI = "showCI" %in% input$options2
    merge='mergeCounts' %in% input$options2
    barchart="barchart" %in% input$options2
    stack = "stacked" %in% input$options2
    group_by=input$group_by
    max_trans = input$maxtrans
    conf.int=input$conf.int
    method="logit";
    toplot = c(isolate(input$toplot5))#,isolate(input$toplot6))#,isolate(input$toplot7),isolate(input$toplot8))
    toplot = toplot[toplot!="-"]
    usegrep=F
    if(length(toplot)==0){
      toplot=c(isolate(input$toplot7),isolate(input$toplot8))
      toplot = toplot[unlist(lapply(toplot,nchar))>2]
      
      usegrep=T
    }
    datafile=session$userData$datafile
   
    #  header = session$userData$header
    total_reads = session$userData$total_reads
    tojoin=isolate(input$tojoin)
    header = names(total_reads)
    if(length(toplot)>0 && !is.null(datafile) ){
      if(usegrep){
        x1 = .findEntries(toplot,datafile,"/trans",tojoin);
        #  merge=TRUE
        #    mat = .readIsoGrep(toplot,datafile,header, "/trans") 
      }else{
        x1 =toplot
        
      }
      mat = t(data.frame( lapply(x1, .readIso, datafile, header, "/trans")))
      if(is.null(dim(mat))) mat = matrix(mat,nrow=1,ncol=length(header))
      if(merge){
        mat = matrix(apply(mat,2,sum),nrow=1,ncol=dim(mat)[2])
      }else if(nchar(group_by)>0){
        groups = .getGroups(x1,group_by)
        toplot=names(groups)
        mat1 = matrix(NA, nrow = length(toplot), ncol  =dim(mat)[2])
        for(j in 1:length(groups)){
          indsj = which(x1 %in% groups[[j]])
          mat1[j,]=apply(mat[indsj,,drop=F],2,sum)
        }
        mat = mat1
      } else{
        
        if(length(x1)>max_trans){
          ord=order(apply(mat,1,sum),decreasing=T)
          mat = mat[ord[1:max_trans],,drop=F]
          toplot = x1[ord[1:max_trans]]
        }else{
          toplot=x1
        }
      }
      #print(toplot)
      
      if(xy){
        subs = list()
        for(i in 1:length(splitby_vec)){
          if(splitby=="molecules"){
            levs1=.getlevels(header,molecules[i], cells, times)
          } else if(splitby=="cells"){
            levs1=.getlevels(header,molecules, cells[i], times)
          }else{
            levs1=.getlevels(header,molecules, cells, times[i])
            
          }
          subs[[i]] =.processTPM(mat, header, toplot, levels=levs1,split=T)
          print(head(subs[[i]]))
          # subs[[i]]$sample=sub(molecules[i],"",subs[[i]]$sample)
        }
        by=c("ID","cell","time")
        if(splitby=="cells"){
          by=c("ID","molecule_type","time") 
        }else if(splitby=="times"){
          by=c("ID","molecule_type","cell")
        }
        subs = merge(subs[[1]],subs[[2]],by=by)
        sample=apply(subs[,2:3,drop=F],1,paste,collapse="_")
      } else if(barchart){
        levs1=.getlevels(header,molecules, cells, times)
        subs =.processTPM(mat, header, toplot, levels=levs1,split=F)
        sample = subs$sample
      }else{
        tpm_df= .processTPM(mat, header, toplot,split=T)
        subs=subset(tpm_df, molecule_type %in% molecules & cell %in% cells & time %in% times)
        sample=apply(subs[,2:4,drop=F],1,paste,collapse="_")
      }
      yname='TPM'
      if(!showTPM){
        yname="Counts"
      }
      if(xy){
        if(splitby=="times") after =FALSE else after=TRUE
        cis.x = .getCIs(subs,sample,total_reads[ grep(splitby_vec[1],names(total_reads))],method, showTPM, prefix=splitby_vec[1],suffix=".x", after=after)
        cis.y = .getCIs(subs,sample,total_reads[ grep(splitby_vec[2],names(total_reads))],method, showTPM, prefix=splitby_vec[2],suffix=".y", after=after)
        cis = cbind(cis.x,cis.y)
        
      }else{
        cis = .getCIs(subs,sample,total_reads,method, showTPM)
        
      }
      subs = cbind(subs,cis)
     
       # resall = 
      #names(resall) =   session$userData$dirname
        session$userData$results = list(data=subs);
      if(xy){
        colorby=names(subs)[1]
        
        shapeby=names(subs)[2] 
        fillby = names(subs)[3]
        #   shapeby=names(subs)[2]
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
        
        if(showCI){
          ggp<-ggp+geom_errorbar(colour="black")
        } #
        trans="identity"
        if(logy){
          trans="log10"
        }
        ggp<-ggp+ scale_y_continuous(trans=trans,name=splitby_vec[2], limits=ylim)+ scale_x_continuous(limits = ylim,trans=trans,name=splitby_vec[1])
        
      }else if(barchart){
        ORF="ID"
        y_text="TPM"
        if(!showTPM) y_text = "Counts";
        #  ord="Start"
        # x1 =  paste("reorder(", ORF, ",", ord,")", sep="") 
        if(stack){
          ggp<-ggplot(subs, aes_string(x="sample",y=y_text,fill=ORF, colour=ORF,ymin="lower" ,ymax="upper"))
          ggp<-ggp+ geom_bar(position="stack", aes_string(y="TPM"),stat="identity")
          
        }else{
          ggp<-ggplot(subs, aes_string(x=ORF,y=y_text,fill="sample", colour='sample',ymin="lower" ,ymax="upper"))
          ggp<-ggp+ geom_bar(position=position_dodge(), aes_string(y="TPM"),stat="identity")
          if(showCI){
            ggp<-ggp+geom_errorbar(position=position_dodge(width=0.9),colour="black")
          } #ggp<-ggp+geom_errorbar(aes_string(x=x1,ymin="lower", ymax="upper"), width=.2)#, position="dodge")
        }
        ggp<-ggp+theme(text = element_text(size=18), axis.text.x = element_text(size = rel(0.7), angle = 25, hjust=0.75))
        
        #geom_bar(aes_string(x=x1, y="Ratio", fill = "type", colour = "type"),stat="identity", position = "dodge")
       
        
        if(logy){
          ggp<-ggp+ scale_y_continuous(trans="log10",name=yname)
        }
        
        
        ggp<-ggp+ggtitle(yname)
        ggp<-ggp+xlab("ORF")
      }else if(!xy){
        if(showCI){
          ggp<-ggplot(subs, aes(x=time, y=TPM ,ymin=lower ,ymax=upper,group=interaction(molecule_type, cell, ID), color = cell, linetype=ID))
          if(ribbon){
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
        ggp<-ggp+theme_bw();#+ylim(c(min(subs$TPM, na.rm=T), max(subs$TPM, na.rm=T)))
        if(!showTPM)ggp<-ggp+ylab("Counts")
        # ggp<-ggp+ geom_errorbar(aes(linetype=molecule_type))
        if(logy){
          ggp<-ggp+ scale_y_log10()
        }
      }
    }else{
      ggp = ggplot()
    }
    ggp
  }
  infectivityPlot=function(){
    currdir = session$userData$currdir
    barchart="barchart" %in% input$options1
    showSecondAxis="showSecondAxis" %in% input$options1
    conf.int=input$conf.int
    
    infilesAnnot = paste(currdir,"0.annot.txt.gz", sep="/")
    total_reads = session$userData$total_reads
    if(file.exists(infilesAnnot)){
      molecules=input$molecules; cells=input$cells; times = input$times;
      type_nme= names(total_reads)
      
      
      levels1 = .getlevels(type_nme,molecules, cells, times)
      orfs=input$orfs
      norm=F
      ratio1 = .readAnnotFile(infilesAnnot,total_reads,norm=norm,levels= levels1, conf.level=conf.int, orfs=orfs)
      #.plotAnnotFile<-function(ratio1,molecule_type, cell, time, levels=NULL,barchart=F,showEB = F,showSecondAxis=F){
      y_text="Ratio"
      #   y_text="spliced"
      annots1=.plotAnnotFile(ratio1,barchart=barchart,showSecondAxis=showSecondAxis,showEB=T, levels=levels1, y_text=y_text)
     
#      resall = 
     # names(resall) =   session$userData$dirname
      session$userData$resultsInf = list(data=annots1$data)
      annots1$ggp
    }else{
      ggplot()
    }
  }
  
  observeEvent(input$plottype,{
    dd= .process1(input$plottype,  session$userData$info)
    updateSelectInput(session,"toplot5", label = dd$label, choices=dd$ch, selected="-")
  })
  observeEvent(input$dir, readDir() )		  
	#THIS JUST EXAMPLE FOR RENDERING A PLOT
	#REACTIVE ON PLOT BUTTON

	output$infPlot<-renderPlot({
	  input$plotButton
	   validate(need(input$dir, 'Please select a directory to begin'))
	  infectivityPlot()
	})

	output$distPlot <- renderPlot({
		validate(need(input$dir, ''))
		#validate(need(input$toplot5 != "-", 'Select a transcript to plot'))
		validate(need(length(input$molecules) > 0 & length(input$cells) > 0 & length(input$times) > 0, 'At least one molecule, cell and time point must be supplied') )
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
#Downloads
output$downloadInf <- downloadHandler(filename = function() {'plotInfectivity.pdf'}, content = function(file) ggsave(file, infectivityPlot(), device='pdf', height = 20, width = 40, units='cm' ) )
output$downloadDepth <- downloadHandler(filename = function() {'plotDepth.pdf'}, content = function(file) ggsave(file, depthPlot("depth"), device = 'pdf', height = 20, width = 40, units='cm') )
output$downloadDepthStart <- downloadHandler(filename = function() {'plotDepthStart.pdf'}, content = function(file) ggsave(file, depthPlot("depthStart"), device = 'pdf', height = 20, width = 40, units='cm') )
output$downloadDepthEnd <- downloadHandler(filename = function() {'plotDepthEnd.pdf'}, content = function(file) ggsave(file, depthPlot("depthEnd"), device = 'pdf', height = 20, width = 40, units='cm') )
output$downloadDist <- downloadHandler(filename = function() {'plotDist.pdf'}, content = function(file) ggsave(file, transcriptPlot(), device = 'pdf' , height = 20, width = 40, units='cm') )
output$downloadResults<-downloadHandler(filename = function() {'results.xlsx'}, content = function(file) write_xlsx( session$userData$results,file ) )
output$downloadResultsInf<-downloadHandler(filename = function() {'resultsInf.xlsx'}, content = function(file) write_xlsx( session$userData$resultsInf,file ) )

	 })

