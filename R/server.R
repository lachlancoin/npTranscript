library(shiny)
library(reshape2)
library(ggplot2)
library(tidyr)
library(rhdf5)
library(RColorBrewer)
library(binom)


source( "transcript_functions.R")
datafile="../data/shiny/0.transcripts.txt.gz"
 data_src ="../data/SARS-Cov2/VIC01"

tpm_df = .readTPM(datafile)
h5file="../data/shiny/0.clusters.h5"      

total_reads = attr(tpm_df,"total_reads")

t = readCoords(paste(data_src, "Coordinates.csv",sep="/"))
subt_inds = t$gene!="none" & t$gene!="leader"
t1 = t[subt_inds,]
dimnames(t1)[[1]] = t[subt_inds,which(names(t)=='gene')]
fimo_file = paste(data_src,"fimo.tsv",sep="/")
  fimo = read.table(fimo_file, sep="\t", head=T)



run_depth<-function(h5file,   toplot=c("leader_leader,N_end"),molecules="RNA",cells="vero",times=c('2hpi','24hpi','48hpi'), span = 0.01, sumAll=F, logy=T, showMotifs=F,showORFs = F){

  	header = h5read(h5file,"header")
	dinds  = 2*(2:(length(header))-2)+2
	type_nme = header[-1]
  types_=data.frame(t(data.frame(strsplit(type_nme,"_"))))
  names(types_) = c("molecules","cell","time")
 inds1 =  which(types_$molecules %in% molecules & types_$cell %in% cells & types_$time %in% times)
 #print(inds1)
 id_cols = c("molecule","cell","time")
   	clusters_ = readH5(h5file, c("pos",header[inds1+1]), toplot,id_cols=id_cols, dinds = dinds[inds1], pos =NULL, span = span, cumul=F, sumAll=sumAll)
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

 plotClusters(clusters_, 2,  1, 
              if(showORFs)t else NULL, 
              if(showMotifs)fimo else NULL,
               rawdepth = rawdepth, linetype=linetype, colour=colour,  xlim = NULL,  title ="depth", logy=logy, leg_size =leg_size1, show=show, fill =fill)
}



# HERE IS THE SERVER PART
##OUTPUT IS PASSED TO THE UI
##INPUT PASSES IN INFORMATION
shinyServer(function(input, output) {
	#output$instructions <- renderPrint({
	#	print("Upload 0.transcripts.txt.gz file produced by npTranscript");
	#})

	#THIS JUST EXAMPLE FOR RENDERING A PLOT
	#REACTIVE ON PLOT BUTTON
	output$distPlot <- renderPlot({
	    input$plotButton
	      #result = loadData();
  	      toplot = c(isolate(input$toplot1),isolate(input$toplot2),isolate(input$toplot3),isolate(input$toplot4))
  	     molecules <-  input$molecules 
  	     cells <- input$cells 
  	     times<-input$times
  	     logy = "logy" %in% input$options
  	     showCI = "showCI" %in% input$options
  	     conf.int=0.95; method="prop.test";
  	     toplot = toplot[toplot!="-"]
  	     if(length(toplot)>0){
  	       subs=subset(tpm_df, ID %in% toplot & molecule_type %in% molecules & cell %in% cells & time %in% times)
  	       sample=apply(subs[,2:4],1,paste,collapse="_")
  	       inds1 =  match(sample, names(total_reads))
  	       subs1 = cbind(subs$count,total_reads[inds1])
  	       nrows = dim(subs)[1]
  	     
  	       cis = matrix(0,ncol=3,nrow =nrows )
  	       dimnames(cis)[[2]] = c("TPM","lower","upper")
  	       
  	       for(i in 1:nrows){
  	         cis[i,1:3] = .calcPropCI(subs1[i,], method=method, conf.int = conf.int)
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
  	     # ggp<-ggp+ geom_errorbar(aes(linetype=molecule_type))
  	      if(logy){
  	        ggp<-ggp+ scale_y_log10()
  	      }
  	      ggp
  	     }
  	  })

	output$depthPlot <- renderPlot({
	    input$plotButton
	      #result = loadData();
	  showDepth  = "show_depth" %in% input$options
	  logy = "logy" %in% input$options
	  showORFs="showORFs" %in% input$options
	  showMotifs="showMotifs" %in% input$options
	  
	    if(showDepth){
	    
	      if(file.exists(h5file)){
	        toplot = c(isolate(input$toplot1),isolate(input$toplot2),isolate(input$toplot3),isolate(input$toplot4))
	        toplot = toplot[toplot!="-"]
	        if(length(toplot)>0){
  		run_depth(h5file,toplot, molecules=input$molecules, cells=input$cells, times = input$times,logy=logy, sumAll = length(toplot)>1,
  		          showORFs = showORFs, showMotifs =showMotifs) 
	        }
  		#run_depth(h5file,toplot=c("leader_leader,N_end")) 
	      }
	    }
	   
	 })

})
