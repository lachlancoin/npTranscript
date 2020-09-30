library(shiny)
library(reshape2)
library(ggplot2)
library(tidyr)

options("np.install"="FALSE")
options("np.libs_to_install"="VGAM,ggplot2,writexl,ggrepel,grDevices,gridExtra,abind,seqinr,RColorBrewer,gplots,seqinr,rhdf5");
options("np.datasource"="~/github/npTranscript/data/SARS-Cov2/VIC01")
options("np.source"="~/github/npTranscript/R" );
options("np.datasource"="~/github/npTranscript/data/SARS-Cov2/VIC01" );
options("np.depth_thresh" = "100" );


options("np.source"="../../R")
options("np.datasource"="../../data/SARS-Cov2/VIC01")
options("np.libdir"="C:/Users/LCOIN/R-4.0.2/library")


source( "transcript_functions.R")


run_analy<-function(transcripts, toplot=c("leader_leader,N_end+", "leader_leader,ORF7a_end+")){
  experiments <- attr(transcripts, 'info')
  count_idx <- which(grepl(pattern = 'count[0-9]', x = colnames(transcripts)))
  props <- prop.table(x = as.matrix(transcripts[,count_idx]), margin = 2)
  colnames(props) <- experiments
  
  #for use in 'split_by', if I can get it to work
  time_vec <- c('2hpi','24hpi','48hpi')
  exp_vec <- c('vero','calu','caco')
  
  #calculate tpm
  tpm <- props*1e-6
  tpm <- cbind(ID=as.character(transcripts$ORFs),tpm)
  
  
  #prep tpm_df
  as.data.frame(tpm, stringsAsFactors=F) %>% melt(id.vars='ID', measure.vars=experiments, value.name = 'TPM') %>%
    separate(variable, c('experiment', 'time'), sep='_', remove = T) %>%
    transform( TPM = as.numeric(TPM), experiment = factor(experiment), time = factor(time, levels = time_vec)) -> tpm_df
  
  #input IDs in vector to plot
  ggplot(subset(tpm_df, ID %in% toplot), aes(x=time, y=TPM, group=interaction(experiment, ID), color = ID)) + geom_line() + scale_y_log10() + geom_point(inherit.aes=T, aes(shape = experiment))
}


#SAVE TRAINED RESULTS
saveData <- function(result) {
  savedResults <<-result
}

#LOAD TRAIN RESULTS
loadData <- function() {
  if (exists("savedResults")) {
    savedResults
  }
}


# HERE IS THE SERVER PART
##OUTPUT IS PASSED TO THE UI
##INPUT PASSES IN INFORMATION
shinyServer(function(input, output) {

	## THIS RUNS FSPLS TRAINING.   IT IS REACTIVE ON PRESSING THE TRAIN BUTTON
	allresults<-reactive({
		input$trainButton
		family = isolate(input$family)
		debug = isolate(input$debug)
		options("debug"= debug);
		text = family;
		 datafile = isolate(input$datafile$datapath[[1]])
		result = NULL
	    if(family!='NA' && !is.null(datafile)){
	   	t = read.table(datafile, sep=",", header=T)
		result = run_analy(t,  family)
	    }
	    result
  	 })


	

	output$instructions <- renderPrint({
		print("Some text");
	})



	#THIS JUST EXAMPLE FOR RENDERING A PLOT
	#REACTIVE ON PLOT BUTTON
	output$distPlot <- renderPlot({
	    input$plotButton
	      #result = loadData();
	      dir="../data/results_20200928213742/"
	      toplot = isolate(input$toplot)
	      transcripts <- .readTranscripts(paste(dir,'0.transcripts.txt.gz',sep="/"))
	  
	      run_analy(transcripts,  toplot)
	   
	 })

})
