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
library(shinydashboard)
library(prompter)
library(shinyBS)
#library(jsonlite);
#library(GGally)

source( "transcript_functions.R")
#.libPaths("C:/Users/LCOIN/R-4.0.2/library")
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
	use_prompt(),
   theme = "https://d2h9b02ioca40d.cloudfront.net/v8.0.1/uom.css",
	tags$head(includeHTML(file.path(basedir, "shiny-common/unset_shiny.html"))),
	tags$head(includeHTML(file.path(basedir, "shiny-common/google-analytics.html"))),
	tags$style(HTML("
	

	 .skin-blue .main-sidebar {
                              background-color: #f5f5f5;
							  padding-top:0px;
							  border:1px solid #e3e3e3;
                              }	
	.skin-blue .main-sidebar .sidebar {
									color:black;
									}
	
	.
	
			"
			
			)),
	htmlTemplate(file.path(basedir, "shiny-common/uomheader.html"),
            title = "Coin Lab",
			apptitle = "SARS-COV-2 Transcriptome",
			subapptitle = "From The Coin Group at University of Melbourne",
			papercitation = 'Chang, J. J.-Y., et al. (2020). "Transcriptional and epi-transcriptional dynamics of SARS-CoV-2 during cellular infection." BioRxiv: 2020.2012.2022.423893.'
               	),
	dashboardPage(
	
header = dashboardHeader(disable = TRUE),
	sidebar = dashboardSidebar(
	 tags$style(type='text/css',
                   ".selectize-dropdown-content{
                 height: 1000px;
                 width: 1000px;
                 background-color: #b0c4de;
                }"),
	
	width=300,
					  
    add_prompt(selectInput("dir", label = "Directory", choices=dirs), position = 'right', message = 'Select directory for desired CoV data'),
    #add_prompt(selectInput("plottype", label = "Transcript category", choices=c("-" )), position = 'right', message = 'Narrow the below field to transcripts of this category'),
	add_prompt(fluidRow( 
		column(2,div(style="display: inline-block; width: 20%;",actionButton("rm_btn", "-"))),
		  column(4,div(style="position-left: 200; display: inline-block; width: 100px; text-align: center; vertical-align: center;", textOutput("counter"))),
		  column(2,div(style="display: linline-block; width: 20%;", actionButton("add_btn", "+")))), position='right', message = 'Choose number of single transcripts to be input'),
	actionButton('print_txn', '?'),
	uiOutput("textbox_ui"),
      
    #add_prompt(selectInput("toplot5", label = "", choices=c("-")), position = 'right', message = 'Select a single transcript to plot'),
    add_prompt(textInput("toplot7", label="All transcripts matching", value = ""), position = 'right', message = 'Select transcripts by R-style regex. Above field must be set to 0 transcripts for this option'),
    add_prompt(selectInput("tojoin", label ="Join", choices=c("AND","OR","AND NOT"), selected="OR"), position = 'right', message = 'How to combine first and second regex fields'),
    add_prompt(textInput("toplot8", label="All transcripts matching", value = ""), position = 'right', message = 'Enter second R-style regex in conjunction with first'),
    add_prompt(selectInput("group_by", label="Group transcripts by", choices = c('No grouping' ,'all', 'type', 'juncts','leader,',',ORF10','ORF1ab,','type:juncts'), selected = 'No grouping'), position = 'right', message = 'Collapse selected transcripts around this grouping'),
    add_prompt(actionButton("plotButton", "Generate plots"), position = 'right', message = 'Run/refresh plots'),
	checkboxGroupInput("molecules", label = "Molecule type",  choices =c()),
   checkboxGroupInput("cells", label = "Cell type"),
   checkboxGroupInput("times", label = "Time points"),
  numericInput("textsize", label = "Text size", value = 20.0, min=3.0, max=100),
  numericInput("angle", label = "Text angle", value = 25.0, min=0.0, max=90),
   checkboxGroupInput("options1", label = h3("Transcript activity panel")) ,
  selectInput("facet1", label ="Grouping", choices=c("off","molecules","cells","times","ORF"), selected="off"),
 
   textInput("orfs", label="ORFs to include"),
   checkboxGroupInput("options2", label = h3("Transcript abundance panel")) ,
   numericInput("conf.int", label = "Confidence intervals", value = 0.95),
   numericInput("maxtrans", label = "Maximum number of transcripts", value = 10),
  #selectInput("splitby", label ="Plot x vs y", choices=c("off","molecules","cells","times"), selected="off"),
  selectInput("facet", label ="Grouping", choices=c("off","molecules","cells","times","molecules_and_cells","molecules_and_times","times_and_cells"), selected="off"),
  
   checkboxGroupInput("options3", label = h3("Depth panel")) ,
 # selectInput("depth_plot_type", label ="What to plot", choices=plot_type_ch, selected="OR"),
  numericInput("min_x", label = "Zoom Min", value = 25000),
  numericInput("max_x", label = "Zoom Max", value = 29865),
  numericInput("loess", label = "Loess span", value = 0.0,max=1,min=0),
 numericInput("waterfallKmer", label = "Kmer for waterfall plot", value = 3,max=11,min=1,step=2),
 numericInput("waterfallOffset", label = "Offset for waterfall plot", value = 0,min=-50,max=50),
 numericInput("maxKmers", label = "Max kmers for waterfall plot", value = 20,max=100,min=2),
 selectInput("test", label ="Test", choices=c("fisher","chisq"), selected="chisq"),
 
 
 numericInput("alpha", label = "Transparency", value = 0.5),
 numericInput("depth_thresh", label = "Depth threshold for errors", value = 100),
 
 #numericInput("linesize", label = "Thickness", value = 0.1, min = 0.0,max=1),
 
 textInput("motif", label="Show motif", value = ""),
 #CTAAAC|TTAAAC
 #ACGAAC|ACGATC|ATGAAC
 
 
 
   add_prompt(checkboxInput('LoadDE', label = 'Activate DE'), position = 'right', message = 'Load DE data from this directory. Button inactive if no data available'),
   selectInput("DE_cell1", label = "DE in cell 1", choices = c()),
   selectInput("DE_cell2", label = "DE in cell 2", choices = c()),
   selectInput("DE_time1", label = "Time 1", choices = c()),
   selectInput("DE_time2", label = "Time 2", choices = c()),
	   checkboxInput("remove_spurious", label = "Remove spurious results"),
	 numericInput("mean_count_thresh", label = "Mean count thresh", value = 0.0,max=500,min=0),
	 textInput("merge_by", label="Transcript group to collapse", value = ""),
     actionButton('plotDE', 'Do DE')
	 ),
	 
	 
  
body = dashboardBody(
useShinyjs(),
 extendShinyjs(text = 'shinyjs.hideSidebar = function(params) { $("body").addClass("sidebar-collapse"); 
              $(window).trigger("resize"); }', functions = 'hideSidebar'),
    extendShinyjs(text='shinyjs.showSidebar = function(params) { $("body").removeClass("sidebar-collapse"); 
                  $(window).trigger("resize"); }', functions = 'showSidebar'),
    bsButton("showpanel", "Show/Hide sidebar",icon = icon("toggle-off"), type = "toggle",style = "info", value = TRUE),

	h2("Transcriptional activity"),
    withSpinner(plotOutput("infPlot", height=400)),
        downloadButton('downloadInf', 'Download plot'),
        downloadButton("downloadResultsInf", 'Download data'),
	h2("Transcript abundance"),
    withSpinner(plotOutput("distPlot", height=400)),
        downloadButton('downloadDist', 'Download plot'),
        downloadButton("downloadResults", 'Download data'),
	h2("Depth of Coverage"),
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
  
)
)
)




