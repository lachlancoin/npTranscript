library(shiny)
library(reshape2)
library(tidyr)
library(rhdf5)
library(RColorBrewer)
library(writexl)

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

options1=c("showCI" ,"barchart", "showSecondAxis")
totick1 = c("showCI" ,"barchart")
#options2 = c("logy","showCI", "TPM" ,"barchart","ribbonCI","mergeCounts")
 options2 = c("logy","showCI", "TPM" ,"barchart","ribbonCI","mergeCounts", "stacked")
totick2 = c("TPM","ribbonCI")
options3 = c("show_depth","logy", "TPM","showMotifs","showORFs", "sumDepth","mergeCounts")
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
    
    selectInput("toplot5", label = paste("Transcript",names(info$choices1)[1]), choices=c("-",info$choices1[[1]]), selected="leader_leader,N_end"),
   # selectInput("toplot6", label = paste("Transcript",names(info$choices1)[1]), choices=c("-",info$choices1[[1]]), selected="-"),
    
    textInput("toplot7", label="All transcripts matching", value = ""),
    selectInput("tojoin", label ="Join", choices=c("AND","OR"), selected="OR"),
    
    textInput("toplot8", label="All transcripts matching", value = ""),
  textInput("group_by", label="Group transcripts by", value = ""),
  
    actionButton("plotButton", "Generate plots"),
   checkboxGroupInput("molecules", label = "Molecule type",  choices =info$molecules, selected = info$molecules),
   checkboxGroupInput("cells", label = "Cell type",  choices = info$cells, selected = info$cells),
   checkboxGroupInput("times", label = "Time points",  choices = info$times, selected = info$times),
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
 numericInput("alpha", label = "Transparency", value = 1.0)
  
  ),
 
  
  # Show a plot of the generated distribution
  mainPanel(
    # verbatimTextOutput("instructions"),
    # verbatimTextOutput("variables"),
    # verbatimTextOutput("validation"),
    plotOutput("infPlot", height=400),
	downloadButton('downloadInf'),
	downloadButton("downloadResultsInf"),
     plotOutput("distPlot", height=400),
	downloadButton('downloadDist'),
	downloadButton("downloadResults"),
   plotOutput("depthPlot", height=400),
   downloadButton("downloadDepth"),
	plotOutput("depthStartPlot", height=400),
	downloadButton("downloadDepthStart"),
	plotOutput("depthEndPlot", height=400),
	downloadButton("downloadDepthEnd")
  )
  #,
  #htmlTemplate(file.path(basedir, "shiny-common/uomfooter.html"))
  
))



