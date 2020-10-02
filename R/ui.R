library(shiny)
source( "transcript_functions.R")

datafile="../data/shiny/0.transcripts.txt.gz"
tpm_df = .readTPM(datafile)

choices=c("-",as.character(attr(tpm_df,"ORFs")))
molecules=levels(tpm_df$molecule_type)

cells = levels(tpm_df$cell)

# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(
 
  
  # Application title
  headerPanel("SARS-COV2 transcriptome"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    #fileInput("datafile", "Transcripts file", multiple = FALSE, accept = NULL),
   selectInput("toplot", label = "Transcript 1", choices=choices, selected=choices[2]),
   selectInput("toplot1", label = "Transcript 2", choices=choices, selected=choices[1]),
   selectInput("toplot2", label = "Transcript 3", choices=choices, selected=choices[1]),
   checkboxGroupInput("molecules", label = h3("Molecules"),  choices =molecules, selected = molecules),
   checkboxGroupInput("cells", label = h3("Molecules"),  choices = cells, selected = cells),
   
                     
	actionButton("plotButton", "Generate plots")
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    # verbatimTextOutput("instructions"),
    # verbatimTextOutput("variables"),
    # verbatimTextOutput("validation"),
     plotOutput("distPlot", height=500)
    # plotOutput("depthPlot", height=500)
  )
))



