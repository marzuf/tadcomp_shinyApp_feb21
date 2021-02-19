#######################################################################################
########################## Genome Partition web application ###########################
########################## assessment of TAD calling methods ##########################
# Marie Zufferey - UNIL - December 2019 ###############################################
# License: Open Source AGPL v3 ########################################################
#######################################################################################

###########################################
##### Server implementation - web app #####
###########################################

# server implementation: we call shinyServer and pass it
# a function that accepts 2 parameters: input and output

# library(BiocManager)
# options(repos = BiocManager::repositories())

library(shiny)
library(shinyBS)
library(shinyjs)
library(graphics)
library(ggplot2)
library(ggpubr)

heightGG <- 7
widthGG <- 7

demo_files <- list.files("data", pattern="_domains.txt", full.names = TRUE)
# demo_files <- demo_files[1:2]

hicDemoResKb <- 25

possible_metrics <- list("Measure of Concordance (MoC)" = "get_MoC", 
                         "Jaccard Index (bins)" = "get_bin_JaccardIndex",
                         "Jaccard Index (boundaries)" = "get_boundaries_JaccardIndex", 
                         "Ratio matching TADs" = "get_ratioMatchingTADs",
                         "Variation of information" = "get_variationInformation")


# plot_TADlist_comparison(TAD_list= list(TopDom=topdom_dt, CaTCH=catch_dt, Arrowhead=arrowhead_dt, HiCseg=hicseg_dt), "get_MoC", nCpu=1)
# plot_TADlist_comparison(TAD_list= list(TopDom=topdom_dt, CaTCH=catch_dt, Arrowhead=arrowhead_dt, HiCseg=hicseg_dt), "get_bin_JaccardIndex", nCpu=1, binSize = bin_size)
# plot_TADlist_comparison(TAD_list= list(TopDom=topdom_dt, CaTCH=catch_dt, Arrowhead=arrowhead_dt, HiCseg=hicseg_dt),"get_boundaries_JaccardIndex", nCpu=1, tolRad = bin_size*2, matchFor="all")
# plot_TADlist_comparison(TAD_list= list(TopDom=topdom_dt, CaTCH=catch_dt, Arrowhead=arrowhead_dt, HiCseg=hicseg_dt),"get_ratioMatchingTADs", nCpu=1, coverMatchRatioThresh = 0.8, matchFor="all")
# plot_TADlist_comparison(TAD_list= list(TopDom=topdom_dt, CaTCH=catch_dt, Arrowhead=arrowhead_dt, HiCseg=hicseg_dt),"get_variationInformation", nCpu=1)

source("scripts/other_metrics.R")
source("scripts/plot_TADlist_similarities.R")

heatmapLowCol <- "blue"
heatmapHighCol <- "red"

shinyServer(function(input, output){
  
  getSimMetric <- function() {
    return(input$simMetric)
  }
  getFieldSep <- function() {
    return(input$fieldSep)
  }
  
  ###############################-----------------------------------------------------------------------------------
  ### Tab with help html file ###-----------------------------------------------------------------------------------
  ###############################-----------------------------------------------------------------------------------
  
  # function to get the file with help information
  getHelp <- function(){
    return(includeHTML("help.html"))
  }
  # display the help in the tab
  output$help <- renderUI({
    getHelp()
  })

  ############################################-----------------------------------------------------------------------------------
  ### Tab with headers of loaded data      ###-----------------------------------------------------------------------------------
  ############################################-----------------------------------------------------------------------------------
  
  output$loaded_data <- renderTable({
    if(input$demo) {
      domain_files <- demo_files
      fieldSep <- "\t"
      fileHeader <- FALSE
      filenames <- basename(domain_files)
    } else{
      domain_files <- input$allDomainFiles$datapath
      fieldSep <- input$fieldSep
      fileHeader <- input$fileHeader
      filenames <- input$allDomainFiles$name
      nameCols <- paste0("dataset", 1:length(domain_files))
    }
    file_heads <- lapply(domain_files, function(x) 
      head(read.table(x, header=fileHeader, sep=fieldSep, stringsAsFactors = FALSE, col.names=c("chromo", "start", "end"))))
    names(file_heads) <- gsub("_final_domains.txt", "", filenames)
    file_heads
  })
  
  #############################################-----------------------------------------------------------------------------------
  ### Tab with HTML with current options    ###-----------------------------------------------------------------------------------
  #############################################-----------------------------------------------------------------------------------
  
  output$showCurrentOptions <- renderUI({

  sepcharact <- input$fieldSep
  if(sepcharact == "\t") {
  	  sepcharact <- "tab"
  } else if(sepcharact == " ") {
	  sepcharact <- "whitespace"
  }

    if(input$demo) {
        filenames <- basename(demo_files)
      } else {
        filenames <- input$allDomainFiles$name
      }
    HTML(paste(
      "<p><u><h5>Demo data ?</u><br>",
      as.character(input$demo),"</p>",
      "<p><u>Similarity metric: </u><br>", 
      names(possible_metrics)[possible_metrics == input$simMetric], "</p>",
      "<p><u>File with header ? </u><br>",
      as.character(input$fileHeader),"</p>",
      "<p><u>Field separator </u><br>",
      as.character(sepcharact),"</p>",
      "<p><u>Datasets: </u><br>",
      paste0(paste0(filenames, collapse="<br>\n")),"</p>"
    ))
  })
  
  
  
  
  #####################################################-----------------------------------------------------------------------------------
  ### Tab with plot descriptive features table      ###-----------------------------------------------------------------------------------
  #####################################################-----------------------------------------------------------------------------------
  descTable <- reactive({
    
    if(input$demo) {
      domain_files <- demo_files
      fieldSep <- "\t"
      fileHeader <- FALSE
      filenames <- gsub("_final_domains.txt", "", basename(demo_files))
    } else{
      domain_files <- input$allDomainFiles$datapath
      fieldSep <- input$fieldSep
      fileHeader <- input$fileHeader
      filenames <- input$allDomainFiles$name
    }
    if(length(domain_files) < 1) {
      return(HTML("! At least one file should be provided ! "))
    }
    
    all_domains_dt <- do.call(rbind, lapply(1:length(domain_files), function(i) {
      curr_dt <- read.table(domain_files[i], header=fileHeader, stringsAsFactors = FALSE, col.names=c("chromo", "start", "end"), sep=fieldSep)
      head(curr_dt)
      curr_dt$dataset <- filenames[i]
      curr_dt
    }))
    
    head(all_domains_dt)
    
    all_domains_dt$size <- all_domains_dt$end - all_domains_dt$start + 1
    all_domains_dt$size_log10 <- log10(all_domains_dt$size)
    
    nbr_dt <- aggregate(size~dataset, FUN=length, data=all_domains_dt)
    
    return(nbr_dt)
    
  })
  descPlot <- reactive({
    
    if(input$demo) {
      domain_files <- demo_files
      fieldSep <- "\t"
      fileHeader <- FALSE
      filenames <- gsub("_final_domains.txt", "", basename(demo_files))
    } else{
      domain_files <- input$allDomainFiles$datapath
      fieldSep <- input$fieldSep
      fileHeader <- input$fileHeader
      filenames <- input$allDomainFiles$name
    }
    if(length(domain_files) < 1) {
      return(HTML("! At least one file should be provided ! "))
    }
    
    all_domains_dt <- do.call(rbind, lapply(1:length(domain_files), function(i) {
      curr_dt <- read.table(domain_files[i], header=fileHeader, stringsAsFactors = FALSE, col.names=c("chromo", "start", "end"), sep=fieldSep)
      head(curr_dt)
      curr_dt$dataset <- filenames[i]
      curr_dt
    }))
    
    head(all_domains_dt)
    
    all_domains_dt$size <- all_domains_dt$end - all_domains_dt$start + 1
    all_domains_dt$size_log10 <- log10(all_domains_dt$size)
    
    nbr_dt <- aggregate(size~dataset, FUN=length, data=all_domains_dt)
    
    p1 <- ggboxplot(data = all_domains_dt, 
                    x = "dataset",
                    y = "size_log10",
                    xlab="",
                    ylab="TAD size [log10")
    p2 <- ggbarplot(data = nbr_dt,
                    x = "dataset",
                    y= "size",
                    ylab="# of TADs")
    
    p_all <- ggarrange(
      p1, p2,
      ncol=2)
    return(p_all)
    
  })
  
  output$descFeatures <- renderPlot({
    
    if(input$demo) {
      domain_files <- demo_files
      fieldSep <- "\t"
      fileHeader <- FALSE
      filenames <- gsub("_final_domains.txt", "", basename(demo_files))
    } else{
      domain_files <- input$allDomainFiles$datapath
      fieldSep <- input$fieldSep
      fileHeader <- input$fileHeader
      filenames <- input$allDomainFiles$name
    }
    if(length(domain_files) < 1) {
      return(HTML("! At least one file should be provided ! "))
    }
    
    all_domains_dt <- do.call(rbind, lapply(1:length(domain_files), function(i) {
      curr_dt <- read.table(domain_files[i], header=fileHeader, stringsAsFactors = FALSE, col.names=c("chromo", "start", "end"), sep=fieldSep)
      head(curr_dt)
      curr_dt$dataset <- filenames[i]
      curr_dt
      }))
    
    head(all_domains_dt)
    
    all_domains_dt$size <- all_domains_dt$end - all_domains_dt$start + 1
    all_domains_dt$size_log10 <- log10(all_domains_dt$size)
    
    nbr_dt <- aggregate(size~dataset, FUN=length, data=all_domains_dt)

    p1 <- ggboxplot(data = all_domains_dt, 
                   x = "dataset",
                   y = "size_log10",
                   xlab="",
                   ylab="TAD size [log10")
    p2 <- ggbarplot(data = nbr_dt,
                    x = "dataset",
                    y= "size",
                    ylab="# of TADs")
    
    p_all <- ggarrange(
      p1, p2,
      ncol=2)
    return(p_all)
    
  })
  
  ###
  ##### Functions to download the desc plot PNG #####
  ###
  output$downloadDescPlot <- downloadHandler(
    filename=paste0("inputData_descPlot.png"),
    content = function(file){
      ggsave(filename = file, plot = descPlot(), height=heightGG, width=widthGG*1.5)
    },
    contentType="image"
  )
  ###
  ##### Functions to download the desc # table TXT #####
  ###
  output$downloadDescTable <- downloadHandler(
    filename=paste0("inputData_descNbr.txt"),
    content = function(file){
      write.table(descTable(), file, col.names=TRUE, row.names=FALSE, quote=FALSE, sep=getFieldSep())
    },
    contentType="text"
  )
  
  
  ############################################-----------------------------------------------------------------------------------
  ### Tab with similarity table      ###-----------------------------------------------------------------------------------
  ############################################-----------------------------------------------------------------------------------
  
  simTable <- reactive({
    
    getMetricArgs <- function(curr_metric) {
      if(curr_metric == "get_MoC") return(NULL)
      if(curr_metric == "get_bin_JaccardIndex") return(list("binSize" = input$bin_size))
      if(curr_metric == "get_boundaries_JaccardIndex") return(list("tolRad" = as.numeric(input$tol_rad), "matchFor"= input$simMatch))
      if(curr_metric == "get_ratioMatchingTADs") return(list("coverMatchRatioThresh" = input$ratioCovMatch, "matchFor"=input$simMatch))
      if(curr_metric == "get_variationInformation") return(NULL)
    } 
    if(input$demo) {
      domain_files <- demo_files
      fieldSep <- "\t"
      fileHeader <- FALSE
      filenames <- basename(demo_files)
    } else{
      domain_files <- input$allDomainFiles$datapath
      fieldSep <- input$fieldSep
      fileHeader <- input$fileHeader
      filenames <- input$allDomainFiles$name
    }
    if(length(domain_files) < 2) {
      return(HTML("! At least two files should be provided ! "))
    }
    all_domains_dt <- lapply(domain_files, function(x) read.table(x, header=fileHeader, stringsAsFactors = FALSE, col.names=c("chromo", "start", "end"), sep=fieldSep))
    names(all_domains_dt) <- gsub("_final_domains.txt", "", filenames)
    
    metricOut <- do.call(plot_TADlist_comparison,
                         c(list(all_domains_dt, metric=input$simMetric, lowColor=heatmapLowCol, highColor=heatmapHighCol, nCpu=input$nCpu),
                           getMetricArgs(input$simMetric)))
    outDt <- metricOut[[2]]
    colnames(outDt)[3] <- names(possible_metrics)[possible_metrics == input$simMetric]
    return(outDt)
  })
  
  output$simTable <- renderTable({
    getMetricArgs <- function(curr_metric) {
      if(curr_metric == "get_MoC") return(NULL)
      if(curr_metric == "get_bin_JaccardIndex") return(list("binSize" = input$bin_size))
      if(curr_metric == "get_boundaries_JaccardIndex") return(list("tolRad" = as.numeric(input$tol_rad), "matchFor"= input$simMatch))
      if(curr_metric == "get_ratioMatchingTADs") return(list("coverMatchRatioThresh" = input$ratioCovMatch, "matchFor"=input$simMatch))
      if(curr_metric == "get_variationInformation") return(NULL)
    } 
    if(input$demo) {
      domain_files <- demo_files
      fieldSep <- "\t"
      fileHeader <- FALSE
      filenames <- basename(demo_files)
    } else{
      domain_files <- input$allDomainFiles$datapath
      fieldSep <- input$fieldSep
      fileHeader <- input$fileHeader
      filenames <- input$allDomainFiles$name
    }
    if(length(domain_files) < 2) {
      return(HTML("! At least two files should be provided ! "))
    }
    all_domains_dt <- lapply(domain_files, function(x) read.table(x, header=fileHeader, stringsAsFactors = FALSE, col.names=c("chromo", "start", "end"), sep=fieldSep))
    names(all_domains_dt) <- gsub("_final_domains.txt", "", filenames)
    
    metricOut <- do.call(plot_TADlist_comparison,
                         c(list(all_domains_dt, metric=input$simMetric, lowColor=heatmapLowCol, highColor=heatmapHighCol, nCpu=input$nCpu),
                           getMetricArgs(input$simMetric)))
    outDt <- metricOut[[2]]
    colnames(outDt)[3] <- names(possible_metrics)[possible_metrics == input$simMetric]
    return(outDt)
  })
  
  ###
  ##### Functions to download the output DT #####
  ###

  output$downloadSimTable <- downloadHandler(
    filename=paste0(getSimMetric(), "_simTable.txt"), # does not work without function call
    content = function(file){
      write.table(simTable(), file, col.names=TRUE, row.names=FALSE, quote=FALSE, sep=getFieldSep())
    },
    contentType = "text"
  )
  
  
  ############################################-----------------------------------------------------------------------------------
  ### Tab with similarity heatmap      ###-----------------------------------------------------------------------------------
  ############################################-----------------------------------------------------------------------------------
  
  simHeatmap <- reactive({
    getMetricArgs <- function(curr_metric) {
      if(curr_metric == "get_MoC") return(NULL)
      if(curr_metric == "get_bin_JaccardIndex") return(list("binSize" = input$bin_size))
      if(curr_metric == "get_boundaries_JaccardIndex") return(list("tolRad" = as.numeric(input$tol_rad), "matchFor"= input$simMatch))
      if(curr_metric == "get_ratioMatchingTADs") return(list("coverMatchRatioThresh" = input$ratioCovMatch, "matchFor"=input$simMatch))
      if(curr_metric == "get_variationInformation") return(NULL)
    } 
    if(input$demo) {
      domain_files <- demo_files
      fieldSep <- "\t"
      fileHeader <- FALSE
      filenames <- basename(demo_files)
    } else{
      domain_files <- input$allDomainFiles$datapath
      fieldSep <- input$fieldSep
      fileHeader <- input$fileHeader
      filenames <- input$allDomainFiles$name
    }
    if(length(domain_files) < 2) {
      return(HTML("! At least two files should be provided ! "))
    }
    all_domains_dt <- lapply(domain_files, function(x) read.table(x, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "start", "end")))
    names(all_domains_dt) <- gsub("_final_domains.txt", "", filenames)
    
    metricOut <- do.call(plot_TADlist_comparison,
                         c(list(all_domains_dt, metric=input$simMetric, lowColor=heatmapLowCol, highColor=heatmapHighCol, nCpu=input$nCpu),
                           getMetricArgs(input$simMetric)))
    simPlot <- metricOut[[1]]
    
    return(simPlot)
  })
  
  
  output$simHeatmap <- renderPlot({
    getMetricArgs <- function(curr_metric) {
      if(curr_metric == "get_MoC") return(NULL)
      if(curr_metric == "get_bin_JaccardIndex") return(list("binSize" = input$bin_size))
      if(curr_metric == "get_boundaries_JaccardIndex") return(list("tolRad" = as.numeric(input$tol_rad), "matchFor"= input$simMatch))
      if(curr_metric == "get_ratioMatchingTADs") return(list("coverMatchRatioThresh" = input$ratioCovMatch, "matchFor"=input$simMatch))
      if(curr_metric == "get_variationInformation") return(NULL)
    } 
    if(input$demo) {
      domain_files <- demo_files
      fieldSep <- "\t"
      fileHeader <- FALSE
      filenames <- basename(demo_files)
    } else{
      domain_files <- input$allDomainFiles$datapath
      fieldSep <- input$fieldSep
      fileHeader <- input$fileHeader
      filenames <- input$allDomainFiles$name
    }
    if(length(domain_files) < 2) {
      return(HTML("! At least two files should be provided ! "))
    }
    all_domains_dt <- lapply(domain_files, function(x) read.table(x, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "start", "end")))
    names(all_domains_dt) <- gsub("_final_domains.txt", "", filenames)
    
    metricOut <- do.call(plot_TADlist_comparison,
                         c(list(all_domains_dt, metric=input$simMetric, lowColor=heatmapLowCol, highColor=heatmapHighCol, nCpu=input$nCpu),
                           getMetricArgs(input$simMetric)))
    simPlot <- metricOut[[1]]
    
    return(simPlot)
  })
  
  
  ###
  ##### Functions to download the output PNG #####
  ###
  output$downloadSimHeatmap <- downloadHandler(
    filename=paste0(getSimMetric(), "_simHeatmap.png"),
    content = function(file){
      ggsave(filename = file, plot = simHeatmap(), width=widthGG, height=heightGG)
    },
    contentType="image"
  )
  
  
  
  ##################################-----------------------------------------------------------------------------------
  ### Enrichment struct. prot. ###-----------------------------------------------------------------------------------
  ##################################-----------------------------------------------------------------------------------
  
  
  structPlot <- function(){
    
    
    if(input$demo) {
      filenames <- basename(demo_files)
      domain_file <- demo_files[which( basename(demo_files) == input$datasetForStructProt)]
    } else {
      filenames <- input$allDomainFiles$name
      domain_file <- input$allDomainFiles$datapath[which( input$allDomainFiles$name == input$datasetForStructProt)]
    }
    
    
    print("hello2\n")
    
    source("scripts/StructProt_EnrichBoundaries_script.R")
    # plot(x=c(1:10), y=c(1:10))
    
    print(paste0("input$datasetForStructProt = ", input$datasetForStructProt))
    print(paste0("input$allDomainFiles$name = ", input$allDomainFiles$name))
    
    
    print(paste0("domain_file = ", domain_file))
    print("hello3")
    
    fieldSep <- input$fieldSep
    fileHeader <- input$fileHeader
    
    if(length(domain_file) < 1) {
      return(HTML("! At least one file should be provided ! "))
    }
    
    print(paste0("domain_file = ", domain_file))
  
    domain_dt <- read.table(domain_file,
                          header=fileHeader, 
                          stringsAsFactors = FALSE,  sep=fieldSep)#col.names=c("chromo", "start", "end"),
    head(domain_dt)
    
    plotTit <- ifelse(input$enrichmentPlotTit == "" | length(input$enrichmentPlotTit) == 0 |is.null(input$enrichmentPlotTit) ,
                      basename(domain_file), input$enrichmentPlotTit)

    print(paste0("plotTit = ", plotTit))
    print(paste0("is.null plotTit = ", is.null(plotTit)))
    print(paste0("empty plotTit = ", plotTit==""))
    
    if(input$demo) {
      hic_res <- hicDemoResKb
    } else {
      hic_res <- input$structProt_resolutionKb
    }
    
    stopifnot(is.numeric(hic_res))
    stopifnot(is.numeric(input$structProt_resStep))
    
     get_structEnrichmentPlot(tad_dt=domain_dt, proteins = input$structProtSelect,
                              resolution_kb=hic_res, 
                              res_step=input$structProt_resStep,
                              plot_tit = plotTit)
    # get_structEnrichmentPlot()
  }
  
  output$structEnrichment <- renderPlot({
    
    print(input$structProtSelect)
    print("hello\n")
    return(structPlot())
  })
  
  
  
  
  
  ###
  ##### Functions to download the output PNG #####
  ###
  output$downloadStructProtEnrichment <- downloadHandler(
    filename=paste0("structProtEnrichmentPlot.svg"),
    content = function(file){
      # ggsave(filename = file, plot = get_structEnrichmentPlot(), width=widthGG, height=heightGG)
      # png(filename = file, plot = get_structEnrichmentPlot(), width=widthGG, height=heightGG)
      svg(filename = file, width=7, height=7)
      # get_structEnrichmentPlot()
      structPlot()
      dev.off()
    },
    contentType="image"
  )
  
  
  ###
  ### Show current datasets for selection
  ###
  
  output$showCurrentDatasets <- renderUI({
    
    if(input$demo) {
      filenames <- basename(demo_files)
    } else {
      filenames <- input$allDomainFiles$name
    }
    if(length(filenames) > 0){
      filenames2 <- paste0(seq_along(filenames), "- ", filenames)      
    } else {
      filenames2 <- ""
    }
    HTML(paste(
      "<p><u>Datasets: </u><br>",
      paste0(paste0(filenames2, collapse="<br>\n")),"</p>"
    ))
  })
  
  output$showSelectedDatasets <- renderUI({
    
    HTML(paste(
      "<p><u>Selected datasets: </u><br>",
      paste0(paste0(input$datasetForStructProt, collapse="<br>\n")),"</p>"
    ))
  })
  
  
  ###
  ### Select current datasets for selection
  ###
  output$secondSelection <- renderUI({
    if(input$demo) {
      filenames <- basename(demo_files)
    } else {
      filenames <- input$allDomainFiles$name
    }
    selectInput("datasetForStructProt", "Dataset:", choices = filenames)
  })
  
  
  
  
  
  ##################################-----------------------------------------------------------------------------------
  ### Tab with contact html file ###-----------------------------------------------------------------------------------
  ##################################-----------------------------------------------------------------------------------
  
  # function to get the file with help information
  getContact <- function(){
    return(includeHTML("contact.html"))
  }
  # display the help in the tab
  output$contact <- renderUI({
    getContact()
  })
  
  ##################################-----------------------------------------------------------------------------------
  ### Tab with details about the metrics currently implemented ###----------------------------------------------------
  ##################################-----------------------------------------------------------------------------------
  
  # display the help in the tab
  output$metricDetails <- renderUI({
    includeHTML("metric_details.html")
  })
  
  ##############################################-----------------------------------------------------------------------------------
  ### Reset selections in the sidebar        ###-----------------------------------------------------------------------------------
  ##############################################-----------------------------------------------------------------------------------
  
  ##### Reset button
  observeEvent(input$reset_button, {
    reset("simMetric")
    reset("fieldSep")
    reset("fileHeader")
    reset("allDomainFiles")
  })
  
  #########################################-----------------------------------------------------------------------------------
  ### Warning if many Cpu to use        ###-----------------------------------------------------------------------------------
  #########################################-----------------------------------------------------------------------------------
  
  observeEvent(input$nCpu, {
    if(input$nCpu > 40)
    showNotification(paste0("Do you really have ", input$nCpu," available ?\n"), type="error")  #  position of the warning set in ui.R (css)
  })
  
})

