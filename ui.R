require(shiny)
require(plotly)
require(tidyverse)
require(DT)
require(shinycssloaders)
require(mixOmics)
require(shinydashboard)

ui <- dashboardPage(
  
#logo position
  dashboardHeader(
    title = img(src = "RZlab.jpg", height = 100, align = "left"),
    titleWidth = 350
  ),
  
  dashboardSidebar(
    width = 350,
      radioButtons("radio", label = h3("Select your proteomics data set"),
                   choices = list("Original ProTargetMiner (n=55)" = 4, "A549 (n=9)" = 1,
                                  "MCF7 (n=9)" = 2, "RKO (n=9)" = 3, "Own data set" = 5), 
                   selected = 4),
      
      conditionalPanel(condition = "input.radio == 5",
                        fileInput("file1", "Please choose csv file",
                                  accept = c(
                                    "text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv"))),
    
    conditionalPanel(condition = "input.radio == 5",
                     checkboxGroupInput("checkGroup", label = h5("Please select data set(s) to merge with"), 
                                        choices = c("Original ProTargetMiner (n=55)" = 4,
                                                    "A549 (n=9)" = 1, "MCF7 (n=9)" = 2, "RKO (n=9)" = 3),
                                        selected = NULL)),

    uiOutput("choose_columns"),
    
    uiOutput("download"),
    
    h1(""),
  
    p(strong(uiOutput("citation")))
    
    ),
 
  #colour and design settings
    dashboardBody(
      tags$head(tags$style(HTML('
                                .skin-blue .main-header .navbar {
                                background-color: #870052;
                                }
                                /* logo */
                                .skin-blue .main-header .logo {
                                background-color: #ffffff ;
                                }
                                /* logo when hovered */
                                .skin-blue .main-header .logo:hover {
                                background-color: #ffffff ;
                                }
                                /* main sidebar */
                                .skin-blue .main-sidebar {
                                background-color: #808080;
                                }                         
                                /* active selected tab in the sidebarmenu */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                                background-color: #808080;
                                }                          
                                /* other links in the sidebarmenu */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                                background-color: #808080;
                                }
                                /* other links in the sidebarmenu when hovered */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                                background-color: #808080;
                                }
                                /* toggle button when hovered  */                    
                                .skin-blue .main-header .navbar .sidebar-toggle:hover{
                                background-color: #870052;
                                }
                                .main-header { max-height: 100px; 
                                font-size:24px; 
                                font-weight:bold; 
                                line-height:24px;
                                }
                                .main-header .logo {
                                height: 100px;
                                font-size:24px; 
                                font-weight:bold; 
                                line-height:24px;align:
                                }
                                .skin-blue .sidebar a {
                                 color: #444;
                                 }
                                .main-sidebar {
                                float:top; margin-top:40px; padding-left:15px; padding-right:15px
                                }
                                '))),
#result section      
      tabsetPanel(type = "tabs",
                  tabPanel("ProTargetMiner", 
                           h1(""),
                           h1(""),
                           h1(textOutput('dynamicTitle1')),
                           dataTableOutput("contents") %>% withSpinner(),
                           
                           h1(""),  
                           h1(""),
                           h1(textOutput('dynamicTitle2')),
                           plotlyOutput("PLSDA") %>% withSpinner(),
                           
                           verbatimTextOutput("click"),
                           
                           h1(""),
                           h1(""),
                           h1(textOutput('dynamicTitle3')),
                           dataTableOutput("top_plsda") %>% withSpinner(),
                           
                           h1(""),
                           h1(""),
                           plotlyOutput("POI")
                           )
      )
    )
)
