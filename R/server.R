library(shiny)
library(reshape2)
library(ggplot2)
library(tidyr)



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
    separate(variable, c('molecule' ,'experiment', 'time'), sep='_', remove = T) %>%
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
	output$instructions <- renderPrint({
		print("Upload 0.transcripts.txt.gz file produced by npTranscript");
	})

	#THIS JUST EXAMPLE FOR RENDERING A PLOT
	#REACTIVE ON PLOT BUTTON
	output$distPlot <- renderPlot({
	    input$plotButton
	      #result = loadData();
	  if(length(input$datafile)>0){
	   datafile=input$datafile$datapath[[1]]
	  }else{
	   datafile="../data/shiny/0.transcripts.txt.gz"
	  }
	      if(file.exists(datafile)){
  	      toplot = isolate(input$toplot)
  	      toplot1 = isolate(input$toplot1)
  	      toplot2 = isolate(input$toplot2)
  	      
  	      transcripts <- .readTranscripts(datafile)
  	  
  	      run_analy(transcripts,  c(toplot, toplot1, toplot2))
	      }
	   
	 })

})
