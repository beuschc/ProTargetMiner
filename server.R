require(shiny)
require(plotly)
require(tidyverse)
require(DT)
require(shinycssloaders)
require(mixOmics)
require(shinydashboard)

server <- function(input, output){
   
  df <- NULL
  selected <- NULL
  opls.df <- NULL
  ds <- NULL
  drug_of_interest <- NULL
  df.merge <- NULL
  inp <- NULL
  
  dt <- data.frame(code = c(1, 2, 3),
                   data.set = c("A459", "MCF7", "RKO"),
                      file.location = c("data/MCf7_dummy.csv", "data/MCf7_dummy.csv", "data/MCf7_dummy.csv"))

  data.set <- reactive({
    if(is.character(input$radio)){
      if(input$radio < 4){
        ds <- as.character(dt$file.location[which(dt$code == input$radio)])
      }
      if(input$radio == 4){
        inFile <- input$file1
        if(!is.null(inFile)){
          ds <- inFile$datapath
        }
      }
    }
    return(ds)
  })
  
  df <- reactive({
    ds <- data.set()
    if(is.character(ds)){
      data <- read.csv2(ds, header = T)
      
      if(input$radio == 4 & length(input$checkGroup) > 0){
        for(i in 1:length(input$checkGroup)){
          ds <- as.character(dt$file.location[which(dt$code == input$checkGroup[i])])
          
          mer <- read.csv2(ds, header = T) %>%
            dplyr::select(-c("Majority.protein.IDs", "Protein.names", "Peptides", "Sequence.coverage...."))
          colnames(mer) <- paste(dt$data.set[which(dt$code == input$checkGroup[i])], colnames(mer), sep = "_")
          colnames(mer)[1] <- "Gene.names"

          suppressWarnings(data <- data %>%
            right_join(mer, by = "Gene.names"))
        }
      }
      
      for(i in 1:ncol(data)){
        w <- which(is.na(data[,i]))
        if(length(w) > 0){
          data <- data[-w,]
        }
        w <- which(data[,i] == Inf) 
        if(length(w) > 0){
          data <- data[-w,]
        }
        w <- which(data[,i] == -Inf) 
        if(length(w) > 0){
          data <- data[-w,]
        }
      }

      if(class(data[, 6]) == "factor"){
        data[, -c(1:3)] <- sapply(data[, -c(1:3)], FUN = as.character)
      }
      if(class(data[, 6]) == "character"){
        data[, -c(1:3)] <- sapply(data[, -c(1:3)], FUN = as.numeric)
      }

    }
    else{
      data <- NULL
    }
    return(data)
  })
  
  output$contents <- DT::renderDataTable({
    d <- df()
    
    if(!is.null(d)){
      good <- which(!grepl("MCF7", colnames(d)) &
                      !grepl("A459", colnames(d)) &
                    !grepl("RKO", colnames(d)))
      datatable(d[,good], options = list(scrollX = T))
    }
    else{
      datatable(data = NULL)
    }
  })
  
  output$choose_columns <- renderUI({
    data <- df()
    if(is.data.frame(data)){
      good <- which(!grepl("MCF7", colnames(data)) &
                      !grepl("A459", colnames(data)) &
                      !grepl("RKO", colnames(data)))
      data <- data[,good]
      data <- data[,-c(1:5)]
      colnames <- str_sub(colnames(data), 1, str_length(colnames(data))-2)
      
      w <- which(colnames == "Control")
      if(length(w) > 0){
        colnames <- colnames[-w]
      }else{
        colnames <- colnames
      }

      selectInput("columns", "Please choose your compound of interest", 
                  choices  = colnames,
                  selected = colnames)
    }
  })
  
  opls.df <- reactive({
    req(input$columns)
    drug_of_interest = input$columns
    data <- df()

    if(is.data.frame(data) & !is.null(drug_of_interest)){
      names <- data$Majority.protein.IDs
      data <- data[,-c(2:5)]
      data <- as.data.frame(t(data[,-1]))
      colnames(data) <- names
      data$Sample <- row.names(data)
      rownames(data) <- NULL
      
      data$Drug <- str_sub(data$Sample, 1, str_length(data$Sample)-2)
      #data <- data %>%
      #  filter(Drug != "Control")
      
      drugs <- data$Drug
      dmatrix <- data %>%
        dplyr::select(-Drug, -Sample)
     
      res <- data.frame(id = names(dmatrix))
      resrank <- data.frame(id = names(dmatrix))
      X = as.matrix(dmatrix)
      
      Y = as.factor(drugs == drug_of_interest)
      colnames(X) <- names(dmatrix)
      plsda.res <- mixOmics::plsda(scale(X), Y, ncomp = 2)
      wyloadings <- plsda.res$loadings$Y[, 1]
      
      load.opls <- as.data.frame(plsda.res$loadings$X)
      load.opls$Majority.protein.IDs <- rownames(load.opls)
      load.opls$label <- "Protein"

      cont <- data.frame("comp 1" = c(min(load.opls$`comp 1`) * 1.5, max(load.opls$`comp 1`) * 1.5),
                         "comp 2" = 0)
      colnames(cont) = c("comp 1", "comp 2")
      
      cont$Majority.protein.IDs = c(drug_of_interest, "all other drugs")
      cont$label <- cont$Majority.protein.IDs

      load.opls <- rbind(load.opls, cont)
      
      load.opls$`comp 1` <- load.opls$`comp 1`
      
      info <- df()
      info <- info[,c(1:5)]
      
      suppressWarnings(load.opls <- info %>%
                         right_join(load.opls, by = "Majority.protein.IDs"))
      
      return(load.opls)
    }
    else{
      return(NULL)
    }
  })
  
  output$OPLS <- renderPlotly({
    req(opls.df())
    load.opls <- opls.df()
    if(is.data.frame(load.opls)){
      load.opls$pointsize <- ifelse(load.opls$label != "Protein", 1.5, 1)
      
      load.opls$ID <- paste(paste("Majority protein IDs", load.opls$Majority.protein.IDs, sep = " = "),
                            paste("Gene names", load.opls$Gene.names, sep = " = "),
                            paste("Protein names", load.opls$Protein.names, sep = " = "),
                            paste("Peptides", load.opls$Peptides, sep = " = "),
                            paste("Sequence coverage", load.opls$Sequence.coverage, sep = " = "), sep = "\n")
      
      g <- ggplot(load.opls) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0) +
        geom_point(aes(x = `comp 1`, y = `comp 2`, colour = `label`, group = `ID`, size = `pointsize`), alpha = 0.5) +
        theme_classic()+
        theme(legend.position = "none")
      
      ggplotly(g, source = "OPLS.plot",
               tooltip = c("group")) %>%
        layout(dragmode = "select")
    }
    else{
      ggplotly(ggplot())
    }
  })
  
  output$POI <- renderPlotly({
    req(opls.df())
    req(input$columns)
    req(df())
    load.opls <- opls.df()
    drug_of_interest <- input$columns
    data <- df()[,-c(2:5)]
    if(is.data.frame(load.opls)){
      poi <- event_data("plotly_click", source = "OPLS.plot")
      
      if(!is.null(poi)){
        poi.df <- load.opls$Majority.protein.IDs[(poi$pointNumber + 1)]
        
        load.opls$ID <- paste(paste("Majority protein IDs", load.opls$Majority.protein.IDs, sep = " = "),
                              paste("Gene names", load.opls$Gene.names, sep = " = "),
                              paste("Protein names", load.opls$Protein.names, sep = " = "),
                              paste("Peptides", load.opls$Peptides, sep = " = "),
                              paste("Sequence coverage", load.opls$Sequence.coverage, sep = " = "), sep = "\n")
        
        load.opls <- load.opls %>%
          filter(Majority.protein.IDs == poi.df)
        
        data <- data %>%
          dplyr::filter(Majority.protein.IDs == poi.df) %>%
          gather(Condition, value, -Majority.protein.IDs)
        
        data$Treatment <- str_sub(data$Condition, 1, str_length(data$Condition)-2)
        
        s <- data %>%
          group_by(Treatment) %>%
          summarise(mean.value = mean(value, na.rm = T),
                    sd.value = sd(value, na.rm = T))
        
        s$mark <- ifelse(s$Treatment == drug_of_interest, 1, 0)
        s$mark <- as.factor(s$mark)
        
        g <- ggplot(s, aes(x = Treatment, y = mean.value, fill = mark)) +
          geom_bar(stat = "identity", position = "dodge") +
          geom_errorbar(aes(ymin = mean.value - sd.value, ymax = mean.value + sd.value), width = 0.25, position = position_dodge(0.1)) +
          ggtitle(load.opls$ID) +
          ylab("mean expression +/- SD") +
          scale_fill_manual(values = c("grey", "red")) +
          ylim(0, max((s$mean.value + s$sd.value)) * 1.2) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "none",
                plot.title = element_text(size = 10))
        
        ggplotly(g, tooltip = c("Treatment", "mean.value"))
      }  
      else{
        ggplotly(ggplot())
      }
    }
    else{
      ggplotly(ggplot())
    }
  })
  
  output$top_opls <- DT::renderDataTable({
    req(opls.df())
    load.opls <- opls.df()
    if(is.data.frame(load.opls)){
      d <- load.opls %>%
        dplyr::arrange(`comp 1`) %>%
        dplyr::filter(label == "Protein") %>%
        dplyr::select(-label)
      
      datatable(d,
                options = list(scrollX = TRUE))
    }
    else{
      return(NULL)
    }
  })
  
  output$dynamicTitle1 <- renderText({
    if(nchar(input$radio) > 0){
      sprintf("Data Table of file %s", data.set())
    }
    else{
      sprintf("Please select your data file")
    }
    
  })
  
  output$dynamicTitle2 <- renderText({
    req(input$columns)
    sprintf("OPLS model for %s", input$columns)
  })
  
  output$dynamicTitle3 <- renderText({
    req(input$columns)
    sprintf("OPLS model for %s", input$columns)
  })
  
  output$dynamicTitle4 <- renderText({
    if(nchar(input$columns) > 0){
      sprintf("OPLS model ranking for %s", input$columns)
    }
    else{
      sprintf("Please select your protein of interest")
    }
    
  })
  
  output$download <- downloadHandler(
    filename = function(){
      file = gsub("\\..*","", data.set())
      paste(file, "_", input$columns, ".tsv", sep = "")
      }, 
    content = function(file.name){
      write_tsv(x = opls.df(), path = file.name)
    }
  )
}