library(shiny)
datafile="../data/shiny/0.transcripts.txt.gz"
transcripts <- .readTranscripts(datafile)
choices=c("-",as.character(transcripts$ORFs))

# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(
 
  
  # Application title
  headerPanel("SARS-COV2 transcriptome"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    fileInput("datafile", "Transcripts file", multiple = FALSE, accept = NULL),
   selectInput("toplot", label = "Transcript 1", choices=choices, selected=choices[2]),
   selectInput("toplot1", label = "Transcript 2", choices=choices, selected=choices[1]),
   selectInput("toplot2", label = "Transcript 3", choices=choices, selected=choices[1]),
	actionButton("plotButton", "Generate plots")
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    # verbatimTextOutput("instructions"),
    # verbatimTextOutput("variables"),
    # verbatimTextOutput("validation"),
     plotOutput("distPlot", height=500)
  )
))



