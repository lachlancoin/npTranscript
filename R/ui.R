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
library(ggforce)
#library(jsonlite);
#library(GGally)

source( "transcript_functions.R")

basedir="../data"
dirs = list.dirs(basedir,full.names=F, rec=T)
required=c("0.isoforms.h5") #,"Coordinates.csv")
for(i in 1:length(required)){
dirs=dirs[which(unlist(lapply(dirs,function(x) file.exists(paste(basedir,x,required[i],sep="/")))))]
}

dirs = c("",dirs)


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

  sidebarPanel(
    selectInput("dir", label = "Directory", choices=dirs),
    selectInput("plottype", label = "Category 1", choices=c("-" )),
    selectInput("toplot5", label = "", choices=c("-")),
    textInput("toplot7", label="All transcripts matching", value = ""),
    selectInput("tojoin", label ="Join", choices=c("AND","OR","AND NOT"), selected="OR"),
    textInput("toplot8", label="All transcripts matching", value = ""),
    selectInput("group_by", label="Group transcripts by", choices = c('No grouping' ,'all', 'type', 'juncts','leader,',',ORF10','ORF1ab,','type:juncts'), selected = 'No grouping'),
       actionButton("plotButton", "Generate plots"),
    checkboxGroupInput("molecules", label = "Molecule type",  choices =c()),
   checkboxGroupInput("cells", label = "Cell type"),
   checkboxGroupInput("times", label = "Time points"),
  numericInput("textsize", label = "Text size", value = 20.0, min=3.0, max=100),
  numericInput("angle", label = "Text angle", value = 25.0, min=0.0, max=90),
   checkboxGroupInput("options1", label = h3("Top panel")) ,
  selectInput("facet1", label ="Grouping", choices=c("off","molecules","cells","times","ORF"), selected="off"),
 
   textInput("orfs", label="ORFs to include"),
   checkboxGroupInput("options2", label = h3("Middle panel")) ,
   numericInput("conf.int", label = "Confidence intervals", value = 0.95),
   numericInput("maxtrans", label = "Maximum number of transcripts", value = 10),
  #selectInput("splitby", label ="Plot x vs y", choices=c("off","molecules","cells","times"), selected="off"),
  selectInput("facet", label ="Grouping", choices=c("off","molecules","cells","times","molecules_and_cells","molecules_and_times","times_and_cells"), selected="off"),
  
   checkboxGroupInput("options3", label = h3("Bottom panel")) ,
 # selectInput("depth_plot_type", label ="What to plot", choices=plot_type_ch, selected="OR"),
  numericInput("min_x", label = "Min position", value = 25000),
  numericInput("max_x", label = "Max position", value = 29865),
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
	downloadButton("downloadResultsDepth", 'Download data'),
	downloadButton("downloadSequence", 'Download sequence'),
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



