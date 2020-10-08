library(shiny)
library(reshape2)
library(tidyr)
library(rhdf5)
library(RColorBrewer)

source( "transcript_functions.R")

basedir="../data"
dirs = list.dirs(basedir,full.names=F, rec=T)
required=c("0.isoforms.h5","Coordinates.csv")
for(i in 1:length(required)){
dirs=dirs[which(unlist(lapply(dirs,function(x) file.exists(paste(basedir,x,required[i],sep="/")))))]
}
seldir=1
dirs = rev(dirs)
#print(seldir)
currdir = paste(basedir,dirs[seldir],sep="/")
#print(currdir)
datafile=paste(currdir,"0.isoforms.h5",sep="/")
h5file=paste(currdir,"0.clusters.h5",sep="/")      


#datafile=paste(currdir,"0.transcripts.txt.gz",sep="/")
toreplace=list(virion="RNA_virion_0hpi", whole_genome_mapped="RNA_vero_24hpi")
#timevec= c('2hpi','24hpi','48hpi')

isoInfo = .getIsoInfo(datafile,h5file, toreplace)
info = .processInfo(isoInfo)
total_reads = isoInfo$total_reads




options=c("show_depth", "logy", "showCI", "TPM","showMotifs","showORFs","barchart", "showSecondAxis","ribbonCI", "sumDepth","mergeCounts")
totick=c("show_depth", "TPM","barchart","ribbonCI")
ch=c(names(info$choices1), names(info$choices))
t=readCoords(paste(currdir, "Coordinates.csv",sep="/"))
orfs=paste(t$gene,collapse=",")
# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(
 
  
  # Application title
  headerPanel("SARS-COV2 transcriptome"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    #fileInput("datafile", "Transcripts file", multiple = FALSE, accept = NULL),\    
    selectInput("dir", label = "Directory", choices=dirs, selected=currdir),
    selectInput("plottype", label = "Category 1", choices=ch, selected=ch[1]),
    selectInput("plottype1", label = "Category 2", choices=ch, selected=ch[1]),
    
    selectInput("toplot5", label = paste("Transcript",names(info$choices1)[1]), choices=c("-",info$choices1[[1]]), selected="-"),
    selectInput("toplot6", label = paste("Transcript",names(info$choices1)[1]), choices=c("-",info$choices1[[1]]), selected="-"),
    textInput("toplot7", label="All transcripts matching", value = ""),
    
    actionButton("plotButton", "Generate plots"),
   checkboxGroupInput("molecules", label = "Molecule type",  choices =info$molecules, selected = info$molecules),
   checkboxGroupInput("cells", label = "Cell type",  choices = info$cells, selected = info$cells),
   checkboxGroupInput("times", label = "Time points",  choices = info$times, selected = info$times),
   checkboxGroupInput("options", label = "Plotting options", choices = options, selected=totick) ,
   textInput("orfs", label="ORFs to include", value = orfs)
 #  numericInput("conf.int", label = h3("Confidence intervals"), value = 0.95),
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    # verbatimTextOutput("instructions"),
    # verbatimTextOutput("variables"),
    # verbatimTextOutput("validation"),
#    plotOutput("exprPlot", height=400),
    plotOutput("infPlot", height=400),
     plotOutput("distPlot", height=400),
   plotOutput("depthPlot", height=400)
  )
))



