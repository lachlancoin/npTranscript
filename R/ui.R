library(shiny)
library(reshape2)
library(tidyr)
library(rhdf5)
library(RColorBrewer)


source( "transcript_functions.R")

basedir="../data"
dirs = list.dirs(basedir,full.names=F, rec=T)
dirs=dirs[which(unlist(lapply(dirs,function(x) file.exists(paste(basedir,x,"0.isoforms.h5",sep="/")))))]

seldir=1
currdir = paste(basedir,dirs[seldir],sep="/")

datafile=paste(currdir,"0.isoforms.h5",sep="/")

#datafile=paste(currdir,"0.transcripts.txt.gz",sep="/")
toreplace=list(virion="RNA_virion_0hpi", whole_genome_mapped="RNA_vero_24hpi")
#timevec= c('2hpi','24hpi','48hpi')
isoInfo = .getIsoInfo(datafile, toreplace)
info = .processInfo(isoInfo)
options=c("show_depth", "logy", "showCI", "TPM","showMotifs","showORFs")
totick=c("show_depth", "TPM")

# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(
 
  
  # Application title
  headerPanel("SARS-COV2 transcriptome"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    #fileInput("datafile", "Transcripts file", multiple = FALSE, accept = NULL),\    
    selectInput("dir", label = "Directory", choices=dirs, selected=currdir),

    selectInput("toplot5", label = paste("Transcript",names(info$choices1)[1]), choices=c("-",info$choices1[[1]]), selected="-"),
    selectInput("toplot6", label = paste("Transcript",names(info$choices1)[2]), choices=c("-",info$choices1[[2]]), selected="-"),
    selectInput("toplot7", label = paste("Transcript",names(info$choices1)[3]), choices=c("-",info$choices1[[3]]), selected="-"),
    selectInput("toplot8", label = paste("Transcript",names(info$choices1)[4]),choices=c("-",info$choices1[[4]]), selected="-"),
    selectInput("toplot1", label = paste("Transcript",names(info$choices)[1]), choices=c("-",info$choices[[1]]), selected="-"),
    selectInput("toplot2", label = paste("Transcript",names(info$choices)[2]), choices=c("-",info$choices[[2]]), selected="-"),
    selectInput("toplot3", label = paste("Transcript",names(info$choices)[3]), choices=c("-",info$choices[[3]]), selected="-"),
    actionButton("plotButton", "Generate plots"),
   checkboxGroupInput("molecules", label = "Molecule type",  choices =info$molecules, selected = info$molecules),
   checkboxGroupInput("cells", label = "Cell type",  choices = info$cells, selected = info$cells),
   checkboxGroupInput("times", label = "Time points",  choices = info$times, selected = info$times),
   checkboxGroupInput("options", label = "Plotting options", choices = options, selected=totick) 
 #  numericInput("conf.int", label = h3("Confidence intervals"), value = 0.95),
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    # verbatimTextOutput("instructions"),
    # verbatimTextOutput("variables"),
    # verbatimTextOutput("validation"),
     plotOutput("distPlot", height=400),
   plotOutput("depthPlot", height=400)
  )
))



