library(shiny)

# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("SARS-COV2 transcriptome"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    #fileInput("datafile", "Transcripts file", multiple = FALSE, accept = NULL),
    textInput("toplot", "Transcripts to plot", value = "leader_leader,N_end+"),#, width = NULL, placeholder = NULL)
    #selectInput("family", label = "Distribution of Y", 
     #   choices = list("--Select--"="NA", "Binomial" = "binomial", "Gaussian" = "gaussian",
      #                 "Multinomial" = "multinomial"), selected = 1),
    #selectInput("debug", label = "Debug mode", 
    #choices = list( "Off" = "off", "On" = "on"  ), selected = 1),
    #actionButton("run", "Plot"),
 #actionButton("testButton", "Test model"),
	actionButton("plotButton", "Generate plots")
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
     verbatimTextOutput("instructions"),
     verbatimTextOutput("variables"),
     verbatimTextOutput("validation"),
     plotOutput("distPlot", height=500)
  )
))



