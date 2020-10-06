library(shiny)
library(reshape2)
library(ggplot2)
library(tidyr)
library(rhdf5)
library(RColorBrewer)
library(binom)



source( "transcript_functions.R")

basedir="../data"
toreplace=list(virion="RNA_virion_0hpi", whole_genome_mapped="RNA_vero_24hpi")

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

#run_depth(h5file, total_reads=total_reads)
run_depth<-function(h5file, total_reads=NULL,  toplot=c("leader_leader,N_end", "N_end"),molecules="RNA",cells="vero",times=c('2hpi','24hpi','48hpi'), 
                    span = 0.01, sumAll=F, fimo=NULL, t= NULL,logy=T, showMotifs=F,showORFs = F){

  	header =.getHeaderH5(h5file,toreplace)
	dinds  = 2*(2:(length(header))-2)+2
	type_nme = header[-1]
  types_=data.frame(t(data.frame(strsplit(type_nme,"_"))))
  names(types_) = c("molecules","cell","time")
 inds1 =  which(types_$molecules %in% molecules & types_$cell %in% cells & types_$time %in% times)
 #print(inds1)
 id_cols = c("molecule","cell","time")
 tot_reads=NULL
 if(!is.null(total_reads)){
   
   tot_reads =  total_reads[inds1]/rep(1e6,length(inds1))
 }
 #print(tot_reads)
   	clusters_ = readH5(h5file,tot_reads, c("pos",header[inds1+1]), toplot,id_cols=id_cols, dinds = dinds[inds1], pos =NULL, span = span, cumul=F, sumAll=sumAll)
#print(clusters_)
   	if(is.null(clusters_)){
  print(paste("could not read ",toplot))
 return (ggplot())
}
   		if(sumAll) type_nme = "all"
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
 plotClusters(clusters_, 2,  1, 
              if(showORFs)t else NULL, 
              if(showMotifs)fimo else NULL,
               rawdepth = rawdepth, linetype=linetype, colour=colour,  xlim = NULL,ylab=ylab , title ="depth", logy=logy, leg_size =leg_size1, show=show, fill =fill)
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
  observeEvent(input$plottype,{
    dd= .process1(input$plottype,  session$userData$info)
    updateSelectInput(session,"toplot5", label = dd$label, choices=dd$ch, selected="-")
  })
  observeEvent(input$plottype1,{
    dd= .process1(input$plottype1, session$userData$info)
  updateSelectInput(session,"toplot6", label = dd$label, choices=dd$ch, selected="-")
  })
  observeEvent(input$dir, {
    print(input$dir)
    print(" updating input dir")
    currdir = paste(basedir,input$dir,sep="/")
    datafile = paste(currdir,"0.isoforms.h5",sep="/")
    h5file=paste(currdir,"0.clusters.h5",sep="/")   
    isoInfo = .getIsoInfo(datafile, h5file,toreplace)
    total_reads = isoInfo$total_reads
    
    ##this gets order by counts
    order = .readTotalIso(datafile, group="/trans", trans=as.character(isoInfo$orfs$ORFs))
    isoInfo$orfs = isoInfo$orfs[order,,drop=F]
    
    info=.processInfo(isoInfo)
    print(paste("set", datafile))
    t = readCoords(paste(currdir, "Coordinates.csv",sep="/"))
    session$userData$t=t
    orfs=paste(t$gene,collapse=",")
    
    fimo_file = paste(currdir,"fimo.tsv",sep="/")
    ch=c(names(info$choices1), names(info$choices))
    
    
    session$userData$currdir=currdir
    session$userData$datafile=datafile
    session$userData$h5file=h5file
    session$userData$fimo = read.table(fimo_file, sep="\t", head=T)
    session$userData$isoInfo=isoInfo
    session$userData$info=info
    session$userData$total_reads = total_reads
    session$userData$header=names(total_reads)

    updateSelectInput(session,"plottype", label = "Category 1", choices=ch, selected=input$plottype)
    updateSelectInput(session,"plottype1", label = "Category 2", choices=ch, selected=input$plottype1)
    updateSelectInput(session, "toplot5",label = paste("Transcript",names(info$choices1)[1]),choices=c("-",info$choices1[[1]]),selected=input$toplot5)
    updateSelectInput(session, "toplot6",label = paste("Transcript",names(info$choices1)[2]),choices=c("-",info$choices1[[2]]),selected=input$toplot6)
    updateCheckboxGroupInput(session,"molecules", label = "Molecule type",  choices =info$molecules, selected = info$molecules)
    updateCheckboxGroupInput(session,"cells", label = "Cell type",  choices = info$cells, selected = info$cells)
    updateCheckboxGroupInput(session,"times", label = "Time points",  choices = info$times, selected = info$times)
    updateTextInput(session,"orfs", label="ORFs to include", value = orfs)
    
  })  
	#THIS JUST EXAMPLE FOR RENDERING A PLOT
	#REACTIVE ON PLOT BUTTON
	output$distPlot <- renderPlot({
	    input$plotButton
  	     molecules <-  input$molecules 
  	     cells <- input$cells 
  	     times<-input$times
  	     logy = "logy" %in% input$options
  	     showCI = "showCI" %in% input$options
  	     showTPM="TPM" %in% input$options
  	     conf.int=0.95; method="prop.test";
  	     toplot = c(isolate(input$toplot5),isolate(input$toplot6))#,isolate(input$toplot7),isolate(input$toplot8))
  	     toplot = toplot[toplot!="-"]
  	     datafile=session$userData$datafile
  	     header = session$userData$header
  	     total_reads = session$userData$total_reads
  	     if(length(toplot)>0 && !is.null(datafile)){
  	       mat = t(data.frame( lapply(toplot, .readIso, datafile, header, "trans")))
  	     #  print(mat)
  	     #  print(header)
  	       tpm_df =.processTPM(mat, header, toplot)
  	   
  	       subs=subset(tpm_df, molecule_type %in% molecules & cell %in% cells & time %in% times)
  	     
  	       sample=apply(subs[,2:4,drop=F],1,paste,collapse="_")
  	       inds1 =  match(sample, names(total_reads))
  	       subs1 = cbind(subs$count,total_reads[inds1])
  	       nrows = dim(subs)[1]
  	       cis = matrix(0,ncol=3,nrow =nrows )
  	       dimnames(cis)[[2]] = c("TPM","lower","upper")
  	      
  	       for(i in 1:nrows){
  	         cis[i,1:3] = .calcPropCI(subs1[i,], method=method, conf.int = conf.int)
  	       
  	       }
  	       if(!showTPM){
  	         cis = cis *subs1[,2]/1e6
  	       }
  	       subs = cbind(subs,cis)
  	     
  	     if(showCI){
  	      ggp<-ggplot(subs, aes(x=time, y=TPM ,ymin=lower ,ymax=upper,group=interaction(molecule_type, cell, ID), color = cell, linetype=ID))
  	      ggp<-ggp+ geom_line(position=position_dodge(width=0.1))  + geom_point(position=position_dodge(width=0.1),inherit.aes=T,aes(shape = molecule_type,size=10))
  	        ggp<-ggp+geom_errorbar(position=position_dodge(width=0.1)) #,colour="black")
  	     }else{
  	      ggp<-ggplot(subs, aes(x=time, y=TPM ,group=interaction(molecule_type, cell, ID), color = cell, linetype=ID))
  	      ggp<-ggp+ geom_line()  + geom_point(inherit.aes=T,aes(shape = molecule_type,size=10))
  	     }
  	     
  	     if(!showTPM)ggp<-ggp+ylab("Counts")
  	     # ggp<-ggp+ geom_errorbar(aes(linetype=molecule_type))
  	      if(logy){
  	        ggp<-ggp+ scale_y_log10()
  	      }
  	      
  	     }else{
  	       ggp = ggplot()
  	     }
  	     ggp
  	  })
output$infPlot<-renderPlot({
  input$plotButton
  currdir = session$userData$currdir
  barchart="barchart" %in% input$options
  showSecondAxis="showSecondAxis" %in% input$options
  infilesAnnot = paste(currdir,"0.annot.txt.gz", sep="/")
  total_reads = session$userData$total_reads
  if(file.exists(infilesAnnot)){
    molecules=input$molecules; cells=input$cells; times = input$times;
  type_nme= session$userData$header
  types_=data.frame(t(data.frame(strsplit(type_nme,"_"))))
  names(types_) = c("molecules","cell","time")
  inds1 =  which(types_$molecules %in% molecules & types_$cell %in% cells & types_$time %in% times)
  types1_ = types_[inds1,,drop=F]
  ord = order(as.numeric(factor(types1_$time, levels=c("0hpi", "2hpi","24hpi","48hpi"))),types1_$cell,types1_$molecules)
##we should use total reads from combined human viral.  For now we set norm=F
  print(input$orfs)
  annots1 = .readAnnotFile(infilesAnnot,total_reads, plot=T,barchart=barchart,norm=F,levels= type_nme[inds1][ord],
                           showSecondAxis=showSecondAxis, annot0 = NULL,conf.level=0.95,showEB=T, orfs=input$orfs)
  annots1
  }else{
    ggplot()
  }
})
	output$depthPlot <- renderPlot({
	    input$plotButton
	      #result = loadData();
	  showDepth  = "show_depth" %in% input$options
	  logy = "logy" %in% input$options
	  showORFs="showORFs" %in% input$options
	  showMotifs="showMotifs" %in% input$options
	  tpm = "TPM" %in% input$options
	  h5file=session$userData$h5file
	  total_reads = NULL
	  fimo = session$userData$fimo
	  t = session$userData$t
  if(tpm){
    total_reads = session$userData$total_reads
  }
	  	#  print(total_reads)
	 # print(h5file)
	    if(showDepth  && !is.null(h5file)){
	    
	      if(file.exists(h5file)){
	       # toplot = c(isolate(input$toplot1),isolate(input$toplot2),isolate(input$toplot3),isolate(input$toplot4))
	        toplot = c(isolate(input$toplot5),isolate(input$toplot6))#,isolate(input$toplot7),isolate(input$toplot8))
	        
	        toplot = toplot[toplot!="-"]
	        if(length(toplot)>0){
	          print(toplot)
  		run_depth(h5file,total_reads,toplot, molecules=input$molecules, cells=input$cells, times = input$times,logy=logy, sumAll = length(toplot)>1,
  		          showORFs = showORFs, fimo=fimo, t=t, showMotifs =showMotifs) 
	        }
  		#run_depth(h5file,toplot=c("leader_leader,N_end")) 
	      }
	    }
	   
	 })

})
