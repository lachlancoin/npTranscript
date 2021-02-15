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
library(ggforce)




# HERE IS THE SERVER PART
##OUTPUT IS PASSED TO THE UI
##INPUT PASSES IN INFORMATION
shinyServer(function(input, output,session) {

	#init
	
	source( "transcript_functions.R")
	source("shiny-DE.R")
  basedir="../data"
  toreplace1 = c("ORF1ab,S_ORF1ab",    "leader,S_ORF1ab",    "ORF1ab,end_3UTR")
  toreplace2= c("ORF1ab,ORF1ab_ORF1ab","leader,ORF1ab_ORF1ab","ORF1ab,end_end")
  decodeFile = paste(basedir,"decode.txt",sep='/')
  reorder = FALSE
  update = TRUE
  debug = FALSE
  counter = list()
  counter$n=0
  
  track_transcripts <- reactiveValues(update=0)
  session$userData$DE=NULL
  shinyjs::disable('LoadDE') 
  shinyjs::disable("add_btn")
  
	#functions

readDir <- function(inputdir, update=T, debug=F) {
		counter$n=0
	  replace=read.table(decodeFile,sep="\t",head=F)
	  toreplace = replace[,2]
	  names(toreplace) = replace[,1]
	  if(debug ) {
	    session = list()
	    if(is.null(session$userData)) session$userData = list()
	  }
	  h5closeAll()
	  print(inputdir)
	  print(nchar(inputdir))
	  if(nchar(inputdir)==0) return(NULL);
	  currdir = paste(basedir,inputdir,sep="/")
	  datafile = paste(currdir,"0.isoforms.h5",sep="/")
	  print(datafile)
	  if(!file.exists(datafile)) return(NULL);
	
	  replace=read.table(decodeFile,sep="\t",head=F)
	  toreplace = replace[,2]
	  names(toreplace) = replace[,1]
	  
    print(" updating input dir")
    session$userData$dataDepth = list("depth"=data.frame(), "depthStart"=data.frame(), "depthEnd"=data.frame())
    session$userData$dataPlot = list("depth"=NULL, "depthStart"=NULL, "depthEnd"=NULL)
    session$userData$depthOpts = list("depth"=list(), "depthStart"=list(), "depthEnd"=list())
   
    h5file=paste(currdir,"0.clusters.h5",sep="/")
    
    if(file.exists(datafile)){
      head =   h5read(datafile,"header")

      head = head[which(!head %in% names(toreplace))]
      if(length(head)>0){
        nme = head
        head =  unlist(lapply(head,.fixName))
        names(head) = nme
        toreplace = c(toreplace,head)
      }
      session$userData$toreplace=toreplace
      isoInfo = .getIsoInfo(datafile, h5file,toreplace)
      total_reads = isoInfo$total_reads
	  shinyjs::enable("add_btn")
      #  #   toreplace1=names(isoInfo$total_reads)
      #   grep("SRR",toreplace1)
      #   names(toreplace) = unlist( lapply(toreplace,.fixName))
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
   annot_file = .findFile(currdir,"annotation.csv.gz") # paste(currdir, ,sep="/")
    if(!is.null(annot_file)){
      annots = read.table(annot_file,sep="\t", head=F)
      names(annots) = c("chr","ens","name","ens1","type")
      session$userData$annots=annots
    }
   
  
   defs = list(
     options1=c("showInfectivity" , "showCI" ,"barchart", "reverseOrder", "showLegend"),
     totick1 = c("showCI" ,"barchart","showInfectivity", "showLegend"),
     options2 = c("showTranscriptPlot","logy","showCI", "TPM_amongst_all" ,"TPM_amongst_viral","barchart","ribbonCI","mergeCounts", "stacked", "reverseOrder", "showLegend"),
     totick2 = c("showTranscriptPlot","ribbonCI","barchart","TPM_amongst_viral", "showLegend"),
     options3 = c("show_depth","logy", "TPM_amongst_viral", "showORFs", "sumDepth","mergeCounts", "showPeptides", "showSequence","showWaterfall", "plotCorr", "showErrors","downsample", "showCI","zoom", "showLegend"),
     totick3 = c("show_depth", "mergeCounts", "sumDepth","zoom", "showLegend")
   )
   
   
   counts_file = paste(currdir, "Counts_genome1.csv",sep="/")
   
   if(file.exists(counts_file)){
      defs$options2 = c("showTranscriptPlot","logy","showCI", "TPM_amongst_all" ,"TPM_amongst_viral","barchart","ribbonCI","mergeCounts", "stacked", "reverseOrder","useReadCount", "showLegend")
      session$userData$counts_file=counts_file
      countsHostVirus = .readCountsHostVirus(counts_file,F)
      session$userData$countsTotal = countsHostVirus[countsHostVirus$ID=="Total",]
   }
     
    coords_file = .findFile(currdir, "Coordinates.csv")
    fimo_file = .findFile(currdir,"fimo.tsv",sep="/")
    
    motifText = ""
    if(!is.null(coords_file)){
      t = readCoords(coords_file)
      session$userData$t=t
      orfs=paste(t$gene,collapse=",")
      fimo=read.table(fimo_file, sep="\t", head=T)
      print(names(fimo))
      motifText = paste(levels(factor(fimo$matched_sequence)),collapse="|")
      session$userData$fimo = fimo
    }else{
      orfs = c()
     
    }
    fastafile = .findFile(currdir,"fasta.gz",grep=T,excl="extra")
   
    if(!is.null(fastafile)){
      
      fastaseq = read.fasta(fastafile)
      session$userData$fasta =fastaseq
      session$userData$fastaseq = lapply(fastaseq,paste,collapse="")
    }else{
      session$userData$fastaseq=NULL
      session$userData$fasta=NULL
    }
    peptide_file =.findFile(currdir,"peptides.csv",sep="/")
    if(!is.null(peptide_file)){
      
      peptide=read.csv(peptide_file,  head=F, comment.char="#")
      trans_vals= as.numeric(sub("#","",read.csv(peptide_file, head=F,nrow=1)))
      peptide = peptide[!duplicated(peptide[,2]),-1,drop=F]
      peptide= apply(peptide,c(1,2),function(x) (x-1)*trans_vals[2]+trans_vals[1])
      names(peptide)=c("start","end")
      session$userData$peptide = peptide
    }
    session$userData$dirname = gsub("/","_",inputdir)
    session$userData$currdir=currdir
    session$userData$datafile=datafile
    session$userData$h5file=h5file
    session$userData$results = list()

    ##NOW UPDATE
    
    if(update){
    	updateCheckboxInput(session, "LoadDE", value = FALSE)
    	hide("plotDE")
    	hide("DE_cell1")
    	hide("DE_cell2")
    	hide("DE_time1")
    	hide("DE_time2")	
  	  session$userData$DE$counts = NULL
    	if (!dir.exists(file.path(currdir,'DE'))) {
    		shinyjs::disable('LoadDE') 
    		} else {
    		session$userData$DE = file.path(currdir, 'DE')
    		shinyjs::enable('LoadDE')
    	}
    
    	#toggleState('ActivateDE', condition = dir.exists(file.path(currdir,'DE')))
		session$userData$transcript_categories = ch
      #updateSelectInput(session,"plottype", label = "Transcript category", choices=ch, selected=input$plottype)
      # updateSelectInput(session, "toplot5",label = paste("Transcript",names(info$choices1)[1]),choices=c("-",info$choices1[[1]]),selected='-')
      updateCheckboxGroupInput(session,"molecules", label = "Molecule type",  choices =info$molecules, selected = info$molecules)
      updateCheckboxGroupInput(session,"cells", label = "Cell type",  choices = info$cells, selected = info$cells)
      updateCheckboxGroupInput(session,"times", label = "Time points",  choices = info$times, selected = info$times)
      updateTextInput(session,"orfs", label="ORFs to include", value = orfs)
      updateTextInput(session,"motif", label="Show motif", value = motifText)
      updateCheckboxGroupInput(session,"options1", label = "Transcript activity panel", choices = defs$options1, selected=defs$totick1) 
      updateCheckboxGroupInput(session,"options2", label = "Transcript abundance panel", choices = defs$options2, selected=defs$totick2) 
      updateCheckboxGroupInput(session,"options3", label = "Depth panel", choices = defs$options3, selected=defs$totick3) 
	  
    }else{
      session$input = list(facet="off",facet1="off",textsize=20,angle=20, min_x = 0, max_x = 29865,
                           group_by = "No grouping",maxtrans=10,tojoin="OR",
                           loess=0, waterfallKmer=3, waterfallOffset=0,maxKmers=20, alpha=0.5, depth_thresh =1000)
      session$input$options2 = c(defs$totick2,"stacked")
      session$input$options3 = defs$totick3
      session$input$options1 = defs$totick1
      session$input$motif=motifText
      session$input$molecules = info$molecules
      session$input$cells = info$cells
      session$input$times = info$times
      session$input$orfs = orfs
      # session$input$toplot5 = "-" 
      # session$input$toplot5="-" 
      session$input$toplot7="5_3" 
      session$input$test="fisher"
    }
    if(debug){
      invisible(session)
    }
  } 
  


loadDE <- function() {
  
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
	
	
	
	
	#plotting functions

infectivityPlot <- function(input){
    if(!"showInfectivity" %in% input$options1) return(ggplot())
    currdir = session$userData$currdir
    barchart="barchart" %in% input$options1
    reverseOrder="reverseOrder" %in% input$options1
    showSecondAxis="showSecondAxis" %in% input$options1
    textsize=input$textsize
    conf.int=input$conf.int
    facet = input$facet1
    showLegend = "showLegend" %in% input$options1
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
      print(type_nme)
      
      levels1 = .getlevels(type_nme,molecules, cells, times,reverseOrder)
      orfs=input$orfs
      norm=F
      ratio1 = .readAnnotFile(infilesAnnot,total_reads1,norm=norm,levels= levels1, conf.level=conf.int, orfs=orfs)
      #.plotAnnotFile<-function(ratio1,molecule_type, cell, time, levels=NULL,barchart=F,showEB = F,showSecondAxis=F){
      y_text="Ratio"
      #   y_text="spliced"
      annots1=.plotAnnotFile(ratio1,barchart=barchart,facet=facet,showSecondAxis=showSecondAxis,showEB=T, levels=levels1, y_text=y_text,
                             diff=0, coeff=5, textsize=textsize)
     if(!showLegend){
       annots1$ggp<-annots1$ggp+ theme(legend.position = "none")
     }
#      resall = 
     # names(resall) =   session$userData$dirname
      session$userData$resultsInf = list(data=annots1$data)
      annots1$ggp
    }else{
      ggplot()
    }
  }
  
 
transcriptPlot<-function(input, selected_transcripts, regex_list){
    if(is.null(session$userData$datafile)) return(ggplot())
    if(!file.exists(session$userData$datafile)) return(ggplot())
    if(!"showTranscriptPlot" %in% input$options2) return(ggplot())

#	p_data= list(molecules=c("RNA","cDNA"),cells=c("calu","vero"),times=c("24hpi","48hpi"),toplot="all",splitby="NA",xy=FALSE,
	            # showTPM=T,merge=F,barchart=T,reverseOrder=T,stack=T,calcTPMFromAll=T,group_by="type",
	           #  tojoin="OR",usegrep=T,merge_by="",method="logit",conf.int=0.95)
#	p_plot = list(textsize=20, logy=T, showCI=F,riboon=F,angle=25,facet="none")
  p_data = list()
  p_plot = list()
	
	p_plot$textsize=input$textsize
  p_data$molecules <-  input$molecules 
  p_data$useReadCount<- "useReadCount" %in% input$options2
  p_data$cells <- input$cells 
  p_data$times<-input$times
  p_data$splitby="off" #input$splitby
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
    p_plot$showLegend = "showLegend" %in% input$options2
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
    useReadCount  = p_data$useReadCount;
	
	#toplot
    p_data$toplot = c(isolate(unique(selected_transcripts)))#,isolate(input$toplot6))#,isolate(input$toplot7),isolate(input$toplot8))
    p_data$toplot = p_data$toplot[p_data$toplot!="-"]
    p_data$usegrep=F
    if(length(p_data$toplot)==0){
      p_data$toplot=c(isolate(regex_list$regex1),isolate(regex_list$regex2))
      p_data$toplot = p_data$toplot[unlist(lapply(p_data$toplot,nchar))>2]
      p_data$usegrep=T
    }
    p_data$tojoin=isolate(regex_list$regex_join)

    if(length(p_data$toplot)==0 || is.null(session$userData$datafile ) ){
      return(ggplot())
    }
	#endtoplot
    datafile = session$userData$datafile ;
    total_reads = session$userData$total_reads
    counts_file=session$userData$counts_file
 
   
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
    sec_axis_name = "Proportion (%)"
    if(reuseData && !is.null(session$userData$results)){
      print("reusing tpm data")
      results_ = session$userData$results
      subs = results_$data
      countsHostVirus1 = results_$totals
    }else{
      
      countsTotal = session$userData$countsTotal
      countsHostVirus1 = NULL
      subs = .extractTPM(datafile ,  total_reads,countsTotal, p_data, toreplace1 = toreplace1, toreplace2 = toreplace2)
      
      if(!is.null(counts_file)){
       
        countsHostVirus = .readCountsHostVirus(counts_file,useReadCount)
        
        countsHostVirus = countsHostVirus[countsHostVirus$ID!="Total",]
        countsHostVirus1 = countsHostVirus[which(countsHostVirus$sample %in% subs$sample),,drop=F]
        names(countsHostVirus1)[2]="Type"
        names(countsHostVirus1)[3]="Reads"
        molecule_type = factor(unlist(lapply(as.character(countsHostVirus1$sample), function(x) strsplit(x,"_")[[1]][1])))
        countsHostVirus1 = .expand(countsHostVirus1,"sample")
        
      }
      
    }
	print(list(data = subs, totals = countsHostVirus1))
    session$userData$results = list(data=subs, totals=countsHostVirus1)
    yname='TPM'
    if(!p_data$showTPM){
          yname="Counts"
    }
      print("plot tpm")
    #  print(head(subs))
      print(paste("sec_axis",sec_axis_name))
      if(useReadCount){
        sec_axis_name="Read count"
      }else{
        sec_axis_name = "Proportion (%) "
      }
    ggp=.plotTPMData(subs,countsHostVirus1, p_data,p_plot,yname, sec_axis_name=sec_axis_name, normaliseToVirus=useReadCount)
    session$userData$tpm_plot = ggp
    ggp
  }
 
 	

depthPlot= function(input, selected_transcripts, regex_list, plot_type, reuse=F) {

 reuseData=F
 if(is.null(session$userData$h5file)) return(ggplot())
   if(!file.exists(session$userData$h5file)) return(ggplot())
	if(is.null(selected_transcripts)) return(ggplot())
 
   
    showDepth  = "show_depth" %in% input$options3
    logy = "logy" %in% input$options3
    textsize=input$textsize
    ci=0
    if(  "showCI" %in% input$options3){
      ci=input$conf.int
    }
      zoom = "zoom" %in% input$options3
    
    depth_thresh = input$depth_thresh
    group_by=input$group_by
    showLegend = "showLegend" %in% input$options3
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
        toplot = unique(c(isolate(selected_transcripts)))#,isolate(input$toplot6))#,isolate(input$toplot7),isolate(input$toplot8))
        tojoin=isolate(regex_list$regex_join)
        merge=F
        toplot = toplot[toplot!="-"]
        #combinedID="combined"
        sumAll = length(toplot)>1
        if(length(toplot)==0){
          ##need to get list from this
          toplot=c(isolate(regex_list$regex1),isolate(regex_list$regex2))
          toplot = toplot[unlist(lapply(toplot,nchar))>2]
          if(length(toplot)>0){
            #combinedID=toplot[1]
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
       plotOpts = list(showORFs = showORFs, motifpos=motifpos,peptides=peptides,xlim =xlim, t=t,showLegend = showLegend,
                       alpha=alpha,linesize=linesize, zoom = zoom,textsize=textsize,logy=logy )
       depthOpts = 
         list(total_reads=total_reads, toplot=toplot, downsample = downsample, span = span, 
              mergeGroups=mergeGroups,molecules=molecules, combinedID='ignore', cells=cells, 
              times = times, sumAll = sumAll,showErrors=showErrors,
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
              toreplace= session$userData$toreplace
          tpm_df=run_depth(h5file,total_reads,toplot, seq_df=seq_df, downsample = downsample, span = span, 
                           mergeGroups=mergeGroups,molecules=molecules, combinedID='ignore', cells=cells, 
                           times = times,logy=logy, sumAll = sumAll,zoom=zoom,
                    showORFs = showORFs, motifpos=motifpos,peptides=peptides,xlim =xlim, t=t,path=plot_type,
                     alpha=alpha,plotCorr=plotCorr,linesize=linesize, reverseOrder=reverseOrder,
                    textsize=textsize, calcErrors=showErrors,fisher=fisher,
                    ci = ci, depth_thresh = depth_thresh,toreplace=toreplace,toreplace1 = toreplace1, toreplace2 = toreplace2,
                    showWaterfall=showWaterfall,waterfallKmer=waterfallKmer,waterfallOffset=waterfallOffset, top10=maxKmers
                    )
            # print(head(tpm_df))
          if(!is.null(tpm_df)){
          session$userData$dataDepth[[which(names(session$userData$dataDepth)==plot_type)]] = tpm_df
          }
              }
            if(!is.null(tpm_df)){
         ggp =plot_depth(tpm_df,total_reads,toplot,zoom=zoom, seq_df=seq_df, downsample = downsample, span = span, mergeGroups=mergeGroups,molecules=molecules, combinedID=combinedID, cells=cells, times = times,logy=logy, sumAll = sumAll,
                     showORFs = showORFs, motifpos=motifpos,peptides=peptides,xlim =xlim, t=t,path=plot_type,
                      alpha=alpha,plotCorr=plotCorr,linesize=linesize, reverseOrder=reverseOrder,
                     textsize=textsize, calcErrors=showErrors,fisher=fisher,
                     ci = ci, depth_thresh = depth_thresh,
                     showWaterfall=showWaterfall,waterfallKmer=waterfallKmer,waterfallOffset=waterfallOffset, top10=maxKmers
          )
       }
         if(!showLegend){
         ggp<-ggp+ theme(legend.position = "none")
         }
          session$userData$dataPlot[[which(names(session$userData$dataDepth)==plot_type)]] =ggp
       }
          return(ggp)
        }
        #run_depth(h5file,toplot=c("leader_leader,N_end")) 
      }
    }
  }



	#observers
observeEvent(input$print_txn, { print(transcript_list())} )


observe({
    if(input$showpanel == TRUE) {
      js$showSidebar()
    }
    else {
      js$hideSidebar()
    }
  })


textboxes <- reactive({

    n <- counter$n

    if (n > 0) {
      lapply(seq_len(n), function(i) {
	  fluidRow(
	 div(style="margin-left: 15px; display: inline-block; width: 40%;",selectInput(inputId = paste0("txncat", i),
                  label = paste0("Category", i), choices = session$userData$transcript_categories) ),
				  
     div(style="display: inline-block; width: 50%;",selectInput(inputId = paste0("txnin", i),
				label = paste0("Transcript", i), choices = "-") )
      )}
	  
    

  )}})

# transcript_list <- reactive( {
	# sapply(seq_len(counter$n), function(i) {
	# input[[paste0("txnin", i)]] } )
	# })


output$textbox_ui <- renderUI({ textboxes() })

	
observeEvent(counter$n, {

track_transcripts$update <- counter$n

lapply(seq_len(counter$n), FUN = function(i) {
	#print(i)
	observeEvent(input[[paste0('txncat',i)]], {
		current_box <- paste0('txncat',i)
		update_box <- paste0('txnin',i)
		print(paste(current_box, update_box))
		updateSelectInput(session, update_box, 
			label = paste('Transcript', input[[current_box]]), 
			choices= .process1(input[[current_box]],  session$userData$info)$ch
			)}
	)}
	)
		})
	
 transcript_list <- eventReactive( input$plotButton, {
	 sapply(seq_len(counter$n), function(i) {
	 input[[paste0("txnin", i)]] } )

		
		
	 })
	
regex_list <- eventReactive( input$plotButton, {
	list(regex1 = input$toplot7, regex2 = input$toplot8, regex_join = input$tojoin) } )




# observeEvent( transcript_list(), {
	# print(paste(transcript_list(), collapse=' '))
	
	# })
	
# observeEvent( transcript_list(), {
# print('cat changed')
# })
  
  
# observeEvent(input$plottype,{
    # dd= .process1(input$plottype,  session$userData$info)
    # updateSelectInput(session,"toplot5", label = dd$label, choices=dd$ch, selected="-")

  # })
  
observeEvent(input$dir, readDir(input$dir, update=T) )
  
  counter <- reactiveValues(n = 0)

  observeEvent(input$add_btn, {
	if(counter$n <10) counter$n <- counter$n + 1
	else warning('maximum transcripts is 10')
	
	})
  observeEvent(input$rm_btn, {
    if (counter$n > 0) counter$n <- counter$n - 1
  })

  output$counter <- renderPrint(cat(counter$n, ' transcripts'))
  
  
  
observeEvent( input$LoadDE,	
{
  if(input$LoadDE){
	print('loadingDE')
	
	
	shinyjs::show("plotDE")
	shinyjs::show("DE_cell1")
	shinyjs::show("DE_cell2")
	shinyjs::show("DE_time1")
	shinyjs::show("DE_time2")
	shinyjs::show("remove_spurious")
	shinyjs::show("mean_count_thresh")
	shinyjs::show("merge_by")
	
	
	session$userData$DE$counts = loadDE()
	#print(head(session$userData$DE$counts))
	
  }else{
    print("de-activating DE")
	print('reset DE')
	shinyjs::hide("plotDE")
	shinyjs::hide("DE_cell1")
	shinyjs::hide("DE_cell2")
	shinyjs::hide("DE_time1")
	shinyjs::hide("DE_time2")
	shinyjs::hide("remove_spurious")
	shinyjs::hide("mean_count_thresh")
	shinyjs::hide("merge_by")
	
	
    session$userData$DE$counts = list()

	
  }
	updateSelectInput(session, "DE_cell1",  choices = names(session$userData$DE$counts))
	updateSelectInput(session, "DE_cell2", choices = names(session$userData$DE$counts))
	updateSelectInput(session, "DE_time1",  choices = session$userData$info$times)
	updateSelectInput(session, "DE_time2",  choices = session$userData$info$times)
	
} )


#DE Plots
observeEvent(input$plotDE, {
print('running DE')
 if (is.null(session$userData$DE$counts)) {print('DE_countdata is null') } else {
	#head(session$userData$DE$counts)
	#not working maybe because of isolate???
   plot_params <- isolate(list(toplot5= '-', toplot2= '', tojoin='', group_by=isolate(input$group_by), merge_by=isolate(input$merge_by)))
	#print(plot_params)
	print('arrived here')
	session$userData$DE$main_out <- tryCatch( {runDE(count_list = session$userData$DE$counts, cell1 = input$DE_cell1 ,
	                     cell2 = input$DE_cell2, time1 = input$DE_time1, 
	                     time2 = input$DE_time2,  thresh=input$mean_count_thresh,
	                     plot_params=isolate(plot_params)) }, error = function(e) message(paste(e)))
						 
	if(!inherits(session$userData$DE$main_out,"try-error")) {
	print(paste('DEout done'))
	print(names(session$userData$DE$main_out))
	
	
	output$DEPlot_volcano <- renderPlot( {
	session$userData$DE$main_out[['volcano_params']][['remove_spurious']] <- input$remove_spurious 
	do.call(volcanoplot, 
	session$userData$DE$main_out[['volcano_params']] )
	})
	
	output$DEPlot_PCA <- renderPlot( {
	do.call(rld_pca, 
	session$userData$DE$main_out[['rld_pca_params']]
	)
	})
	} else {print('try-error caught')}
	
		}
	 } )


	#render plots

	output$infPlot<-renderPlot({
	   shiny::validate(need(input$dir, 'Please select a directory to begin'))
	  input$plotButton
	  infectivityPlot(input)
	})

	output$distPlot <- renderPlot({
		shiny::validate(need(input$dir, ''))
		shiny::validate(need(any(transcript_list() != "-") | regex_list()$regex1 != "", 'Select transcripts to plot'))
		shiny::validate(need(length(input$molecules) > 0 & length(input$cells) > 0 & length(input$times) > 0, 'At least one molecule, cell and time point must be supplied') )
	    input$plotButton
  	    transcriptPlot(isolate(input), isolate(transcript_list()), isolate(regex_list()))
  	  })

	output$depthPlot <- renderPlot({
	shiny::validate(need(any(transcript_list() != "-") | regex_list()$regex1 != "", ''))
	input$plotButton
	depthPlot(isolate(input), isolate(transcript_list()), isolate(regex_list()), "depth")
	})
	output$depthStartPlot <- renderPlot({
	shiny::validate(need(any(transcript_list() != "-") | regex_list()$regex1 != "", ''))
	  input$plotButton
	  depthPlot(isolate(input), isolate(transcript_list()), isolate(regex_list()), "depthStart")
	})
	output$depthEndPlot <- renderPlot({
	shiny::validate(need(any(transcript_list() != "-") | regex_list()$regex1 != "", ''))
	  input$plotButton
	  depthPlot(isolate(input), isolate(transcript_list()), isolate(regex_list()), "depthEnd")
		
	})

	 
#Downloads
output$downloadInf <- downloadHandler(filename = function() {'plotInfectivity.pdf'}, content = function(file) ggsave(file, infectivityPlot(input), device='pdf', height = 20, width = 40, units='cm' ) )
output$downloadDepth <- downloadHandler(filename = function() {'plotDepth.pdf'}, content = function(file) ggsave(file, depthPlot(input, transcript_list(),regex_list(), "depth", T), device = 'pdf', height = 20, width = 40, units='cm') )
output$downloadDepthStart <- downloadHandler(filename = function() {'plotDepthStart.pdf'}, content = function(file) ggsave(file, depthPlot(input, transcript_list(),regex_list(), "depthStart", T), device = 'pdf', height = 20, width = 40, units='cm') )
output$downloadDepthEnd <- downloadHandler(filename = function() {'plotDepthEnd.pdf'}, content = function(file) ggsave(file, depthPlot(input, transcript_list(), regex_list(),"depthEnd",T), device = 'pdf', height = 20, width = 40, units='cm') )
output$downloadDist <- downloadHandler(filename = function() {'plotDist.pdf'}, content = function(file) ggsave(file, transcriptPlot(input, transcript_list(), regex_list()), device = 'pdf' , height = 20, width = 40, units='cm') )
output$downloadResults<-downloadHandler(filename = function() {'results.xlsx'}, content = function(file) write_xlsx( session$userData$results,file ) )
output$downloadResultsInf<-downloadHandler(filename = function() {'resultsInf.xlsx'}, content = function(file) write_xlsx( session$userData$resultsInf,file ) )
output$downloadResultsDepth<-downloadHandler(filename = function() {'resultsDepth.xlsx'}, content = function(file) write_xlsx( session$userData$dataDepth,file ) )
output$downloadSequence<-downloadHandler(filename = function() {'sequence.fa'}, content = function(file) write.fasta(extract_sequence(),extract_sequence_name(), file ) )



output$downloadPCA <- downloadHandler(filename = function() {'plotPCA.pdf'}, content = function(file) {
	pdf(file, 18, 18, pointsize=50)
	do.call(rld_pca, session$userData$DE$main_out[['rld_pca_params']])
	dev.off()
		})
output$downloadVOLCANO <- downloadHandler(filename = function() {'plotVOLCANO.pdf'}, content = function(file) {
	pdf(file, 18, 18, pointsize=20)
	do.call(volcanoplot, session$userData$DE$main_out[['volcano_params']] )
	dev.off()
		})
output$downloadDEdata <- downloadHandler(filename = function() {'DE_data.xlsx'}, content = function(file) {
	write_xlsx(session$userData$DE$main_out[['data']], file)
	})






##THIS IS FOR DEBUGGING
  if(FALSE){
    inputdir = "SARS-Cov2/VIC01" #
    inputdir="229E_new"
    session=readDir(inputdir,update=F,debug=T)
   # session$input$options2
   # session$input$group_by="type"
    session$input$toplot7="leader_leader,M_3UTR"
    session$input$toplot8="leader_leader,N_3UTR"
    session$input$options3 = c("showErrors",grep("mergeCounts" ,session$input$options3,v=T,inv=T))
   # session$input$group_by=
    infectivityPlot(session$input)
   # transcriptPlot(session$input, selected_transcripts = c(session$input$toplot7, session$input$toplot8), regex_list())
 #   session$input$selected_transcripts = c("^leader_leader,M_3UTR$","^leader_leader,N_3UTR$")
    transcriptPlot(session$input)
  #  depthPlot(session$input, "depth")
  
    depthPlot(session$input, selected_transcripts = c(session$input$toplot7, session$input$toplot8), regex_list = list(regex1 = '', regex2='', to_join=''), "depth")
  }

})
