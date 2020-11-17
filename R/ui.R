library(shiny)
library(reshape2)
library(tidyr)
library(rhdf5)
library(RColorBrewer)
library(writexl)
library(shinycssloaders)
library(shinyjs)
library(abind)
library(ggrepel)
#library(GGally)

source( "transcript_functions.R")

basedir="../data"
dirs = list.dirs(basedir,full.names=F, rec=T)
required=c("0.isoforms.h5") #,"Coordinates.csv")
for(i in 1:length(required)){
dirs=dirs[which(unlist(lapply(dirs,function(x) file.exists(paste(basedir,x,required[i],sep="/")))))]
}
seldir=1
dirs = c("",dirs)
#print(seldir)
currdir = paste(basedir,dirs[seldir],sep="/")
#print(currdir)
datafile=paste(currdir,"0.isoforms.h5",sep="/")
h5file=paste(currdir,"0.clusters.h5",sep="/")      
decodeFile = paste(basedir,"decode.txt",sep='/')
replace=read.table(decodeFile,sep="\t",head=F)
toreplace = replace[,2]

names(toreplace) = replace[,1]


#datafile=paste(currdir,"0.transcripts.txt.gz",sep="/")

#timevec= c('2hpi','24hpi','48hpi')
if(file.exists(datafile)){
  isoInfo = .getIsoInfo(datafile,h5file, toreplace)
  info = .processInfo(isoInfo)
  total_reads = isoInfo$total_reads
  ch=c(names(info$choices1), names(info$choices))
  
}else{
  isoInfo=NULL
  info=NULL
  total_reads=NULL
  ch = c()
}
if(file.exists(h5file)){
plot_type_ch =   sub("/","",grep("depth",unique(h5ls(h5file)[,1]),v=T))
}else{
  plot_type_ch  = c("-");
}

options1=c("showInfectivity" , "showCI" ,"barchart", "reverseOrder") #"showSecondAxis", 
totick1 = c("showCI" ,"barchart","showInfectivity")
#options2 = c("logy","showCI", "TPM" ,"barchart","ribbonCI","mergeCounts")
 options2 = c("showTranscriptPlot","logy","showCI", "TPM_amongst_all" ,"TPM_amongst_viral","barchart","ribbonCI","mergeCounts", "stacked", "reverseOrder")
totick2 = c("showTranscriptPlot","ribbonCI","barchart","TPM_amongst_viral")
options3 = c("show_depth","logy", "TPM_amongst_viral", "showORFs", "sumDepth","mergeCounts", "showPeptides", "showSequence","showWaterfall", "plotCorr", "showErrors","downsample", "showCI")
totick3 = c("show_depth", "mergeCounts", "sumDepth")

coordsFile = paste(currdir, "Coordinates.csv",sep="/")
if(file.exists(coordsFile)){
  t=readCoords(coordsFile)
  orfs=paste(t$gene,collapse=",")
}else{
  orfs = c()
}
# Define UI for application that plots random distributions 
shinyUI(fluidPage(
	useShinyjs(),
   theme="https://d2h9b02ioca40d.cloudfront.net/v8.0.1/uom.css",
	tags$head(includeHTML(file.path(basedir, "shiny-common/unset_shiny.html"))),
	htmlTemplate(file.path(basedir, "shiny-common/uomheader.html"),
            title = "Coin Lab",
			apptitle = "SARS-COV-2 Transcriptome",
			subapptitle = ""
               	),

  # Application title
  #headerPanel("SARS-COV2 transcriptome"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    #fileInput("datafile", "Transcripts file", multiple = FALSE, accept = NULL),\    
    selectInput("dir", label = "Directory", choices=dirs, selected=currdir),
    selectInput("plottype", label = "Category 1", choices=ch, selected=ch[1]),
  #  selectInput("plottype1", label = "Category 2", choices=ch, selected=ch[1]),
    
    selectInput("toplot5", label = paste("Transcript",names(info$choices1)[1]), choices=c("-",info$choices1[[1]]), selected="-"),
   # selectInput("toplot6", label = paste("Transcript",names(info$choices1)[1]), choices=c("-",info$choices1[[1]]), selected="-"),
    
    textInput("toplot7", label="All transcripts matching", value = ""),
    selectInput("tojoin", label ="Join", choices=c("AND","OR","AND NOT"), selected="OR"),
    
    textInput("toplot8", label="All transcripts matching", value = ""),
  selectInput("group_by", label="Group transcripts by", choices = c('No grouping' ,'all', 'type', 'juncts','leader,',',ORF10','ORF1ab,','type:juncts'), selected = 'No grouping'),
       actionButton("plotButton", "Generate plots"),
   checkboxGroupInput("molecules", label = "Molecule type",  choices =info$molecules, selected = info$molecules),
   checkboxGroupInput("cells", label = "Cell type",  choices = info$cells, selected = info$cells),
   checkboxGroupInput("times", label = "Time points",  choices = info$times, selected = info$times),
  numericInput("textsize", label = "Text size", value = 20.0, min=3.0, max=100),
   checkboxGroupInput("options1", label = h3("Top panel"), choices = options1, selected=totick1) ,
   textInput("orfs", label="ORFs to include", value = orfs),
   checkboxGroupInput("options2", label = h3("Middle panel"), choices = options2, selected=totick2) ,
   numericInput("conf.int", label = "Confidence intervals", value = 0.95),
   numericInput("maxtrans", label = "Maximum number of transcripts", value = 10),
  selectInput("splitby", label ="Plot x vs y", choices=c("off","molecules","cells","times"), selected="off"),
  
   checkboxGroupInput("options3", label = h3("Bottom panel"), choices = options3, selected=totick3) ,
 # selectInput("depth_plot_type", label ="What to plot", choices=plot_type_ch, selected="OR"),
  numericInput("min_x", label = "Min position", value = 0),
  numericInput("max_x", label = "Max position", value = 30000),
  numericInput("loess", label = "Loess span", value = 0.0,max=1,min=0),
 numericInput("waterfallKmer", label = "Kmer for waterfall plot", value = 3,max=11,min=1,step=2),
 numericInput("waterfallOffset", label = "Offset for waterfall plot", value = 0,min=-50,max=50),
 numericInput("maxKmers", label = "Max kmers for waterfall plot", value = 20,max=100,min=2),
 selectInput("test", label ="test", choices=c("fisher","chisq"), selected="chisq"),
 
 numericInput("alpha", label = "Transparency", value = 0.5),
 numericInput("depth_thesh", label = "Depth threshold for errors", value = 1000),
 
 #numericInput("linesize", label = "Thickness", value = 0.1, min = 0.0,max=1),
 
 textInput("motif", label="Show motif", value = ""),
 #CTAAAC|TTAAAC
 #ACGAAC|ACGATC|ATGAAC
 
 
 
   checkboxInput('LoadDE', label = 'Activate DE'),
   selectInput("DE_cell1", label = "DE in cell 1", choices = c()),
   selectInput("DE_cell2", label = "DE in cell 2", choices = c()),
   selectInput("DE_time1", label = "Time 1", choices = c()),
   selectInput("DE_time2", label = "Time 2", choices = c()),
   checkboxInput("remove_spurious", label = "Remove spurious results"),
 numericInput("mean_count_thresh", label = "Mean count thresh", value = 0.0,max=500,min=0),
 textInput("merge_by", label="Transcript group to collapse", value = ""),
     actionButton('plotDE', 'Do DE')
  ),
 
  
  # Show a plot of the generated distribution
  mainPanel(
    # verbatimTextOutput("instructions"),
    # verbatimTextOutput("variables"),
    # verbatimTextOutput("validation"),
	h2("Infectivity Plots"),
    withSpinner(plotOutput("infPlot", height=400)),
        downloadButton('downloadInf', 'Download plot'),
        downloadButton("downloadResultsInf", 'Download data'),
	h2("TPM Plots"),
    withSpinner(plotOutput("distPlot", height=400)),
        downloadButton('downloadDist', 'Download plot'),
        downloadButton("downloadResults", 'Download data'),
	h2("Depth Plots"),
    plotOutput("depthPlot", height=400),
        downloadButton("downloadDepth", 'Download plot'),
    plotOutput("depthStartPlot", height=400),
        downloadButton("downloadDepthStart", 'Download plot'),
    plotOutput("depthEndPlot", height=400),
        downloadButton("downloadDepthEnd", 'Download plot'),
	h2("DE Plots"),
	plotOutput("DEPlot_PCA", height = 400),
	downloadButton('downloadPCA', 'Download PCA'),
	plotOutput("DEPlot_volcano", height = 400),
	downloadButton("downloadVOLCANO", 'Download Volcano'),
	downloadButton("downloadDEdata", 'Download DE data')
	

  )
  #,
  #htmlTemplate(file.path(basedir, "shiny-common/uomfooter.html"))
  
))



