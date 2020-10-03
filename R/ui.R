library(shiny)
library(reshape2)
library(tidyr)
library(rhdf5)
library(RColorBrewer)


source( "transcript_functions.R")

datafile="../data/shiny/0.transcripts.txt.gz"
ORFs = .getORFs(datafile)
exps = attr(ORFs,"experiments")
choices = lapply(levels(ORFs$num_breaks), function(x) as.character(ORFs$ORFs[ORFs$num_breaks==x]))
names(choices) = levels(ORFs$num_breaks)
choices[[4]] = as.character(unlist(choices[-(1:3)]))
choices = choices[1:4]
molecules=levels(exps$molecule_type)
times = levels(exps$time)

cells = levels(exps$cell)
options=c("show_depth", "logy", "showCI", "showMotifs","showORFs")
# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(
 
  
  # Application title
  headerPanel("SARS-COV2 transcriptome"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    #fileInput("datafile", "Transcripts file", multiple = FALSE, accept = NULL),
    
    selectInput("toplot1", label = paste("Transcript",names(choices)[1]), choices=c("-",choices[[1]]), selected="-"),
    selectInput("toplot2", label = paste("Transcript",names(choices)[2]), choices=c("-",choices[[2]]), selected="-"),
    selectInput("toplot3", label = paste("Transcript",names(choices)[3]), choices=c("-",choices[[3]]), selected="-"),
    selectInput("toplot4", label = paste("Transcript",names(choices)[4]),choices=c("-",choices[[4]]), selected="-"),
    
   checkboxGroupInput("molecules", label = "Molecule type",  choices =molecules, selected = molecules),
   checkboxGroupInput("cells", label = "Cell type",  choices = cells, selected = cells),
   checkboxGroupInput("times", label = "Time points",  choices = times, selected = times),
   checkboxGroupInput("options", label = "Plotting options", choices = options, selected=options[options!="logy"]) ,
 #  numericInput("conf.int", label = h3("Confidence intervals"), value = 0.95),
	actionButton("plotButton", "Generate plots")
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    # verbatimTextOutput("instructions"),
    # verbatimTextOutput("variables"),
    # verbatimTextOutput("validation"),
     plotOutput("distPlot", height=300),
   plotOutput("depthPlot", height=300)
  )
))



