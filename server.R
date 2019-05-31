library(shiny)
library(plotly)
library(tidyverse)
library(DT)
library(shinycssloaders)
library(mixOmics)
library(shinydashboard)

server <- function(input, output){
  options(shiny.maxRequestSize = 50*1024^2)
  
  #defining variables
  df <- NULL
  selected <- NULL
  plsda.df <- NULL
  ds <- NULL
  drug_of_interest <- NULL
  df.merge <- NULL
  inp <- NULL
  
  #code, name and location of all available data sets
  dt <- data.frame(code = c(1, 2, 3, 4),
                   data.set = c('A549', 'MCF7', 'RKO', 'ProTargetMiner'),
                   file.location = c('data/A549_deep_proteome.csv', 'data/MCF7_deep_proteome.csv', 'data/RKO_deep_proteome.csv', 'data/ProTargetMiner.csv'))
  
  #assign location (name) to the select input
  data.set <- reactive({
    if(is.character(input$radio)){
      if(input$radio < 5){
        ds <- as.character(dt$file.location[which(dt$code == input$radio)])
      }
      if(input$radio == 5){
        inFile <- input$file1
        if(!is.null(inFile)){
          ds <- inFile$datapath
        }
      }
    }
    return(ds)
  })
  
  # read csv file(s) either as part of the privided data bases or from user
  df <- reactive({
    ds <- data.set()
    if(is.character(ds)){
      data <- suppressMessages(read_csv(ds))
      
      if(input$radio == 5 & length(input$checkGroup) > 0){
        for(i in 1:length(input$checkGroup)){
          ds <- as.character(dt$file.location[which(dt$code == input$checkGroup[i])])
          
          mer <- suppressMessages(read_csv(ds)) %>%
            dplyr::select(-c('Majority protein IDs', 'Protein names', 'Peptides', 'Sequence coverage [%]'))
          colnames(mer) <- paste(as.character(dt$data.set[which(dt$code == input$checkGroup[i])]), colnames(mer), sep = '_')
          colnames(mer)[1] <- 'Gene names'
          
          suppressWarnings(data <- data %>%
                             left_join(mer, by = 'Gene names') %>%
                             na.omit()
          )
          
          data <- data[!duplicated(data$`Majority protein IDs`),]
        }
      }
    }
    else{
      data <- NULL
    }
    return(data)
  })
  
  #table containing the user select file(s)
  output$contents <- DT::renderDataTable({
    d <- df()
    
    if(!is.null(d)){
      good <- which(!grepl('MCF7', colnames(d)) &
                      !grepl('A549', colnames(d)) &
                      !grepl('RKO', colnames(d)))
      datatable(d[,good], options = list(scrollX = T))
    }
    else{
      datatable(data = NULL)
    }
  })
  
  #generation of all allowed drugs for plsda calculation
  output$choose_columns <- renderUI({
    data <- df()
    if(is.data.frame(data)){
      good <- which(!grepl('MCF7', colnames(data)) &
                      !grepl('A549', colnames(data)) &
                      !grepl('RKO', colnames(data)))
      data <- data[,good]
      data <- data[,-c(1:5)]
      colnames <- str_sub(colnames(data), 1, str_length(colnames(data))-2)
      
      selectInput('columns', 'Please choose your compound of interest', 
                  choices  = colnames,
                  selected = colnames)
    }
  })
  
  #generation of plsda model
  #all drug contrasting the user defined drug
  plsda.df <- reactive({
    req(input$columns)
    drug_of_interest = input$columns
    data <- df()
    
    if(is.data.frame(data) & !is.null(drug_of_interest)){
      names <- data$`Majority protein IDs`
      data <- data[,-c(2:5)]
      data <- as.data.frame(t(data[,-1]))
      colnames(data) <- names
      data$Sample <- row.names(data)
      #rownames(data) <- NULL
      
      data$Drug <- str_sub(data$Sample, 1, str_length(data$Sample)-2)
      
      drugs <- data$Drug
      dmatrix <- data %>%
        dplyr::select(c(-Drug, -Sample))
      
      res <- data.frame(id = names(dmatrix))
      resrank <- data.frame(id = names(dmatrix))
      X = as.matrix(dmatrix)
      
      Y = as.factor(drugs == drug_of_interest)
      colnames(X) <- names(dmatrix)
      plsda.res <- mixOmics::plsda(scale(X), Y, ncomp = 2)
      wyloadings <- plsda.res$loadings$Y[, 1]
      
      load.plsda <- as.data.frame(plsda.res$loadings$X)
      load.plsda$`Majority protein IDs` <- rownames(load.plsda)
      load.plsda$label <- 'Protein'
      
      orientation_min <- load.plsda[which(load.plsda$`comp 1` == min(load.plsda$`comp 1`, na.rm = T)),]
      orientation_max <- load.plsda[which(load.plsda$`comp 1` == max(load.plsda$`comp 1`, na.rm = T)),]
      
      d <- df() %>%
        dplyr::select(c(`Majority protein IDs`, contains(drug_of_interest)))
      
      d_min <- d %>%
        dplyr::filter(`Majority protein IDs` == orientation_min$`Majority protein IDs`) %>%
        gather(Treatment, value, 2:4) %>%
        summarise(mean.value = mean(value, na.rm = T)) 
      
      d_max <- d %>%
        dplyr::filter(`Majority protein IDs` == orientation_max$`Majority protein IDs`) %>%
        gather(Treatment, value, 2:4) %>%
        summarise(mean.value = mean(value, na.rm = T))
      
      #define orientation of plsda model and plot
      #positive x-axis corresponds to specific upregulation in user definied drug
      if(d_min > d_max){
        load.plsda$`comp 1` <- - load.plsda$`comp 1`
      }
      
      cont <- data.frame('comp 1' = c(min(load.plsda$`comp 1`) * 1.5, max(load.plsda$`comp 1`) * 1.5),
                         'comp 2' = 0)
      colnames(cont) = c('comp 1', 'comp 2')
      
      cont$`Majority protein IDs` = c('all other drugs', drug_of_interest)
      cont$label <- cont$`Majority protein IDs`
      
      load.plsda <- rbind(load.plsda, cont)
      
      load.plsda$`comp 1` <- load.plsda$`comp 1`
      
      info <- df()
      info <- info[,c(1:5)]
      
      suppressWarnings(load.plsda <- info %>%
                         right_join(load.plsda, by = 'Majority protein IDs'))
      
      return(load.plsda)
    }
    else{
      return(NULL)
    }
  })
  
  #intarctive plot from plsda model
  output$PLSDA <- renderPlotly({
    req(plsda.df())
    load.plsda <- plsda.df()
    if(is.data.frame(load.plsda)){
      load.plsda$pointsize <- ifelse(load.plsda$label != 'Protein', 1.5, 1)
      
      load.plsda$ID <- paste(paste('Majority protein IDs', load.plsda$`Majority protein IDs`, sep = ' = '),
                             paste('Gene names', load.plsda$`Gene names`, sep = ' = '),
                             paste('Protein names', load.plsda$`Protein names`, sep = ' = '),
                             paste('Peptides', load.plsda$Peptides, sep = ' = '),
                             paste('Sequence coverage', load.plsda$`Sequence coverage [%]`, sep = ' = '), sep = '\n')
      
      g <- ggplot(load.plsda) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0) +
        geom_point(aes(x = `comp 1`, y = `comp 2`, colour = `label`, group = `ID`, size = `pointsize`), alpha = 0.5) +
        theme_classic()+
        theme(legend.position = 'none')
      
      ggplotly(g, source = 'PLSDA.plot',
               tooltip = c('group')) %>%
        layout(dragmode = 'select')
    }
    else{
      ggplotly(ggplot())
    }
  })
  
  #expression proteomics for the protein selected in the plsda plot
  output$POI <- renderPlotly({
    req(plsda.df())
    req(input$columns)
    req(df())
    load.plsda <- plsda.df()
    drug_of_interest <- input$columns
    data <- df()[,-c(2:5)]
    if(is.data.frame(load.plsda)){
      poi <- event_data('plotly_click', source = 'PLSDA.plot')
      
      if(!is.null(poi)){
        poi.df <- load.plsda$`Majority protein IDs`[(poi$pointNumber + 1)]
        
        load.plsda <- load.plsda %>%
          filter(`Majority protein IDs` == poi.df)
        
        data <- data %>%
          dplyr::filter(`Majority protein IDs` == poi.df) %>%
          gather(Condition, value, -`Majority protein IDs`)
        
        data$Treatment <- str_sub(data$Condition, 1, str_length(data$Condition)-2)
        
        marked.protein <- data %>%
          filter(Treatment == drug_of_interest) %>%
          mutate(log.value = log2(value))
        
        load.plsda <- load.plsda %>%
        mutate(ID_main =  paste(paste("Majority protein IDs", load.plsda$`Majority protein IDs`, sep = " = "),
                               paste("Gene names", load.plsda$`Gene names`, sep = " = "),
                               paste("Protein names", load.plsda$`Protein names`, sep = " = "),
                               paste("Peptides", load.plsda$Peptides, sep = " = "),
                               paste("Sequence coverage = ", load.plsda$`Sequence coverage [%]`, "%", sep = ""),
                               paste("p.value vs. control", round(t.test(marked.protein$log.value, mu = 0)$p.value, 4), sep = " = "),
                               sep = "\n"))
        s <- data %>%
          group_by(Treatment) %>%
          mutate(log.value = log2(value)) %>%
          summarise(mean.value = mean(log.value, na.rm = T),
                    sd.value = sd(log.value, na.rm = T))
        
        s$mark <- ifelse(s$Treatment == drug_of_interest, 1, 0)
        s$mark <- as.factor(s$mark)
        
        g <- ggplot(s, aes(x = Treatment, y = mean.value, fill = mark)) +
          geom_bar(stat = 'identity', position = 'dodge') +
          geom_errorbar(aes(ymin = mean.value - sd.value, ymax = mean.value + sd.value), width = 0.25, position = position_dodge(0.1)) +
          geom_hline(yintercept = 0, col = 'black', linetype = 'solid', alpha = 0.5) +
          labs(title = load.plsda$ID_main) +
          ylab('Log2 Mean Fold Change +/- SD') +
          scale_fill_manual(values = c('grey', 'red')) +
          ylim(-Inf, max((s$mean.value + s$sd.value)) * 1.4) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = 'none',
                plot.title = element_text(size = 8, hjust = 0.5))
        
        ggplotly(g, tooltip = c('Treatment', 'mean.value'))
      }  
      else{
        ggplotly(ggplot())
      }
    }
    else{
      ggplotly(ggplot())
    }
  })
  
  #interactive table for the plsda results
  output$top_plsda <- DT::renderDataTable({
    req(plsda.df())
    load.plsda <- plsda.df()
    if(is.data.frame(load.plsda)){
      d <- load.plsda %>%
        dplyr::arrange(`comp 1`) %>%
        dplyr::filter(label == 'Protein') %>%
        dplyr::select(-label)
      
      datatable(d,
                options = list(scrollX = TRUE))
    }
    else{
      return(NULL)
    }
  })
  
  #user interaction with instruction for next step
  output$dynamicTitle1 <- renderText({
    if(nchar(input$radio) > 0){
      sprintf('Data Table of file %s', data.set())
    }
    else{
      sprintf('Please select your data file')
    }
    
  })
  
  output$dynamicTitle2 <- renderText({
    req(input$columns)
    sprintf('PLSDA model for %s', input$columns)
  })
  
  output$dynamicTitle3 <- renderText({
    req(input$columns)
    sprintf('PLSDA model for %s', input$columns)
  })
  
  #user interaction with instruction for next step
  output$dynamicTitle4 <- renderText({
    if(nchar(input$columns) > 0){
      sprintf('PLSDA model ranking for %s', input$columns)
    }
    else{
      sprintf('Please select your protein of interest')
    }
    
  })
  
  #download of plsda results with automatic naming
  output$download <- downloadHandler(
    filename = function(){
      file = gsub('\\..*','', data.set())
      paste(file, '_', input$columns, '.tsv', sep = '')
    }, 
    content = function(file.name){
      write_tsv(x = plsda.df() %>%
                  filter(label == 'Protein'),
                path = file.name)
    }
  )
  
  #link to paper
  url <- a('ProTargetMiner: A proteome signature library of anticancer molecules for functional discovery',
           href = 'https://www.biorxiv.org/content/10.1101/421115v1')
  output$citation <- renderUI({
    tagList('Please cite:', url)
  })
}
