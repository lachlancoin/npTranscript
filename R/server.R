library(shiny)
library(reshape2)
library(ggplot2)
library(tidyr)
library(rhdf5)



source( "transcript_functions.R")
datafile="../data/shiny/0.transcripts.txt.gz"

tpm_df = .readTPM(datafile)
h5file="../data/shiny/0.clusters.h5"      






run_depth<-function(h5file,   toplot=c("leader_leader,N_end"),molecules="RNA",cells="vero",times=c('2hpi','24hpi','48hpi'), span = 0.01, sumAll=F, logy=T){
	header = h5read(h5file,"header")
	dinds  = 2*(2:length(header)-2)+2
	t = NULL; fimo=NULL;
	type_nme = header[-1]
  types_=data.frame(t(data.frame(strsplit(type_nme,"_"))))
  names(types_) = c("molecules","cell","time")
 inds1 =  which(types_$molecules %in% molecules & types_$cell %in% cells & types_$time %in% times)
   	clusters_ = readH5(h5file, c("pos",header[inds1]), toplot, dinds = dinds[inds1], pos =NULL, span = span, cumul=F, sumAll=sumAll)
	if(sumAll) type_nme = "all"
	#mat=t(h5read(h5file,paste("depth", toplot[1],sep="/")))
rawdepth = T
leg_size1=10
show=T
fill=F
k = 1
linetype="clusterID"
colour="sampID"
if(sumAll){
colour="clusterID"
linetype="sampID"
}
 plotClusters(clusters_, 2,  1, t, fimo, rawdepth = rawdepth, linetype=linetype, colour=colour,  xlim = NULL,  title ="depth", logy=logy, leg_size =leg_size1, show=show, fill =fill)
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
  	      toplot = c(isolate(input$toplot),isolate(input$toplot1),isolate(input$toplot2))
  	     molecules <-  input$molecules 
  	     cells <- input$cells 
  	    
  	      ggplot(subset(tpm_df, ID %in% toplot & molecule_type %in% molecules & cell %in% cells & time %in% input$times), aes(x=time, y=TPM, group=interaction(molecule_type, cell, ID), color = ID, linetype=molecule_type)) + geom_line() + scale_y_log10() + geom_point(inherit.aes=T,aes(shape = cell))
	 })

	output$depthPlot <- renderPlot({
	    input$plotButton
	      #result = loadData();
	  
	      if(file.exists(h5file)){
	        toplot = c(isolate(input$toplot),isolate(input$toplot1),isolate(input$toplot2))
  		run_depth(h5file,toplot, molecules=input$molecules, cells=input$cells, times = input$times) 
  		#run_depth(h5file,toplot=c("leader_leader,N_end")) 
	      }
	   
	 })

})
