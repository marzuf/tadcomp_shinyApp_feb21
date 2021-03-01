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

# demo_files <- list.files("data", pattern="_domains.txt", full.names = TRUE)
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

proteinPeaks_Folder = "data/struct_prot/"

demoFile_ctcf <- file.path(proteinPeaks_Folder, "chr6_CTCF_peaks.bed")
demoFile_rad21 <- file.path(proteinPeaks_Folder, "chr6_RAD21_peaks.bed")
demoFile_smc3 <- file.path(proteinPeaks_Folder, "chr6_SMC3_peaks.bed")

source("scripts/other_metrics.R")
source("scripts/plot_TADlist_similarities.R")
source("scripts/StructProt_EnrichBoundaries_script.R")

heatmapLowCol <- "blue"
heatmapHighCol <- "red"

shinyServer(function(input, output){
  
  get_demoFiles <- function() {
    print(paste0("input$demoBinSize = ", input$demoBinSize))
    print(paste0("input$demoNorm = ", input$demoNorm))
    print(paste0("input$demoSub = ", input$demoSub))
    if(input$demoNorm=="ICE" & input$demoBinSize == "50kb"){
        all_files <- list.files(file.path("data","cp_data",input$demoBinSize, input$demoNorm, paste0("sub",input$demoSub)),
                          pattern="_domains.txt", full.names = TRUE)
    } else {
      all_files <- list.files(file.path("data","cp_data",input$demoBinSize, input$demoNorm),
                              pattern="_domains.txt", full.names = TRUE)
    }
    return(all_files)
  }
  
  # getSimMetric <- function() {
  #   return(input$simMetric)
  # }
  getSimMetric_table <- function() {
    return(input$table_simMetric)
  }
  getSimMetric_heatmap <- function() {
    return(input$heatmap_simMetric)
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
      demo_files <- get_demoFiles()
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
  output$loaded_data_chip <- renderTable({
    # if(input$demo) {
    #   demo_files <- get_demoFiles()
    #   domain_files <- demo_files
    #   fieldSep <- "\t"
    #   fileHeader <- FALSE
    #   filenames <- basename(domain_files)
    # } else{
    #   domain_files <- input$allDomainFiles$datapath
    #   fieldSep <- input$fieldSep
    #   fileHeader <- input$fileHeader
    #   filenames <- input$allDomainFiles$name
    #   nameCols <- paste0("dataset", 1:length(domain_files))
    # }
    chip_files <-  get_structProtFiles()
    
    print(paste0("chip_files = ", chip_files))
    
    
    chip_files <- Filter(function(x)!is.null(x), chip_files)
    
    if(length(chip_files) > 0) {
      filenames <- basename(unlist(chip_files))
      file_heads <- lapply(seq_along(chip_files), function(i) 
        # head(read.table(x, header=fileHeader, sep=fieldSep, stringsAsFactors = FALSE, col.names=c("chromo", "start", "end"))))
        head(read.table(chip_files[[i]], quote = '', stringsAsFactors=FALSE, sep = "\t", header = FALSE)))
      names(file_heads) <- gsub("_peaks.bed", "", filenames)
      file_heads
    }
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
        demo_files <- get_demoFiles()
        print(paste0("demo_files = ", demo_files))
        filenames <- basename(demo_files)
      } else {
        filenames <- input$allDomainFiles$name
      }
  
  filenames_chip <- get_structProtFiles()
  
    HTML(paste(
      "<p><u><h5>Demo data ?</u><br>",
      as.character(input$demo),"</p>",
      # "<p><u>Similarity metric: </u><br>", 
      # names(possible_metrics)[possible_metrics == input$simMetric], "</p>",
      "<p><u>Similarity metric - table: </u><br>", 
      names(possible_metrics)[possible_metrics == input$table_simMetric], "</p>",
      "<p><u>Similarity metric - heatmap: </u><br>", 
      names(possible_metrics)[possible_metrics == input$heatmap_simMetric], "</p>",
      "<p><u>File with header ? </u><br>",
      as.character(input$fileHeader),"</p>",
      "<p><u>Field separator </u><br>",
      as.character(sepcharact),"</p>",
      "<p><u>Hi-C datasets: </u><br>",
      paste0(paste0(filenames, collapse="<br>\n")),"</p>",
      "<p><u>ChIP-seq datasets: </u><br>",
      paste0(paste0(names(filenames_chip), ": ", basename(unlist(filenames_chip)), collapse="<br>\n")),"</p>"
    ))
  })
  
  
  
  
  #####################################################-----------------------------------------------------------------------------------
  ### Tab with plot descriptive features table      ###-----------------------------------------------------------------------------------
  #####################################################-----------------------------------------------------------------------------------
  descTable <- reactive({
    
    print("DESCTABLE")
    
    if(input$demo) {
      demo_files <- get_demoFiles()
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
    print(paste0("domain_files = ", domain_files))
    if(length(domain_files) < 1) {
      return(data.frame("no data provided"))
      # return(HTML("! At least one file should be provided ! "))
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
      demo_files <- get_demoFiles()
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
      # return(HTML("! At least one file should be provided ! "))
      return(data.frame("no data provided"))
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
      demo_files <- get_demoFiles()
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
  
  ##### to display message if no data selected:
  # https://stackoverflow.com/a/21535587/1100107
  getData <- reactive({
    if(is.null(input$allDomainFiles)) return(NULL)
    return(1)
  })
  output$domainFilesUploaded <- reactive({
    print(paste0("input$allDomainFiles = ", input$allDomainFiles))
    return(!is.null(getData()))
  })
  outputOptions(output, 'domainFilesUploaded', suspendWhenHidden=FALSE)
  
  ### for the chipseq
  getDataChip <- reactive({
    if(is.null(input$smc3File) & 
       is.null(input$ctcfFile) &  is.null(input$rad21File)) return(NULL)
    return(1)
  })
  output$chipFilesUploaded <- reactive({
    return(!is.null(getDataChip()))
  })
  outputOptions(output, 'chipFilesUploaded', suspendWhenHidden=FALSE)
  
  ############################################-----------------------------------------------------------------------------------
  ### Tab with similarity table      ###-----------------------------------------------------------------------------------
  ############################################-----------------------------------------------------------------------------------
  
  simTable <- reactive({
    
    getMetricArgs <- function(curr_metric) {
      if(curr_metric == "get_MoC") return(NULL)
      # if(curr_metric == "get_bin_JaccardIndex") return(list("binSize" = input$bin_size))
      # if(curr_metric == "get_boundaries_JaccardIndex") return(list("tolRad" = as.numeric(input$tol_rad), "matchFor"= input$simMatch))
      # if(curr_metric == "get_ratioMatchingTADs") return(list("coverMatchRatioThresh" = input$ratioCovMatch, "matchFor"=input$simMatch))
      if(curr_metric == "get_bin_JaccardIndex") return(list("binSize" = input$table_bin_size))
      if(curr_metric == "get_boundaries_JaccardIndex") return(list("tolRad" = as.numeric(input$table_tol_rad), "matchFor"= input$table_simMatch))
      if(curr_metric == "get_ratioMatchingTADs") return(list("coverMatchRatioThresh" = input$table_ratioCovMatch, "matchFor"=input$table_simMatch))
      
      if(curr_metric == "get_variationInformation") return(NULL)
    } 
    if(input$demo) {
      demo_files <- get_demoFiles()
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
      # if(curr_metric == "get_bin_JaccardIndex") return(list("binSize" = input$bin_size))
      # if(curr_metric == "get_boundaries_JaccardIndex") return(list("tolRad" = as.numeric(input$tol_rad), "matchFor"= input$simMatch))
      # if(curr_metric == "get_ratioMatchingTADs") return(list("coverMatchRatioThresh" = input$ratioCovMatch, "matchFor"=input$simMatch))
      if(curr_metric == "get_bin_JaccardIndex") return(list("binSize" = input$table_bin_size))
      if(curr_metric == "get_boundaries_JaccardIndex") return(list("tolRad" = as.numeric(input$table_tol_rad), "matchFor"= input$table_simMatch))
      if(curr_metric == "get_ratioMatchingTADs") return(list("coverMatchRatioThresh" = input$table_ratioCovMatch, "matchFor"=input$table_simMatch))
      
      if(curr_metric == "get_variationInformation") return(NULL)
    } 
    if(input$demo) {
      demo_files <- get_demoFiles()
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
                         # c(list(all_domains_dt, metric=input$simMetric, lowColor=heatmapLowCol, highColor=heatmapHighCol, nCpu=input$nCpu),
                         #   getMetricArgs(input$simMetric)))
                         c(list(all_domains_dt, metric=input$table_simMetric, lowColor=heatmapLowCol, highColor=heatmapHighCol, nCpu=input$nCpu),
                           getMetricArgs(input$table_simMetric)))
    
    outDt <- metricOut[[2]]
    # colnames(outDt)[3] <- names(possible_metrics)[possible_metrics == input$simMetric]
    colnames(outDt)[3] <- names(possible_metrics)[possible_metrics == input$table_simMetric]
    return(outDt)
  })
  
  ###
  ##### Functions to download the output DT #####
  ###
  
  output$downloadSimTable <- downloadHandler(
    # filename=paste0(getSimMetric(), "_simTable.txt"), # does not work without function call
    filename=paste0(getSimMetric_table(), "_simTable.txt"), # does not work without function call
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
      # if(curr_metric == "get_bin_JaccardIndex") return(list("binSize" = input$bin_size))
      # if(curr_metric == "get_boundaries_JaccardIndex") return(list("tolRad" = as.numeric(input$tol_rad), "matchFor"= input$simMatch))
      # if(curr_metric == "get_ratioMatchingTADs") return(list("coverMatchRatioThresh" = input$ratioCovMatch, "matchFor"=input$simMatch))
      if(curr_metric == "get_bin_JaccardIndex") return(list("binSize" = input$heatmap_bin_size))
      if(curr_metric == "get_boundaries_JaccardIndex") return(list("tolRad" = as.numeric(input$heatmap_tol_rad), "matchFor"= input$heatmap_simMatch))
      if(curr_metric == "get_ratioMatchingTADs") return(list("coverMatchRatioThresh" = input$heatmap_ratioCovMatch, "matchFor"=input$heatmap_simMatch))
      
      if(curr_metric == "get_variationInformation") return(NULL)
    } 
    if(input$demo) {
      demo_files <- get_demoFiles()
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
                         #                      c(list(all_domains_dt, metric=input$simMetric, lowColor=heatmapLowCol, highColor=heatmapHighCol, nCpu=input$nCpu),
                         # getMetricArgs(input$simMetric)))
                         c(list(all_domains_dt, metric=input$heatmap_simMetric, lowColor=heatmapLowCol, highColor=heatmapHighCol, nCpu=input$nCpu),
                           getMetricArgs(input$heatmap_simMetric)))
    
    simPlot <- metricOut[[1]]
    
    return(simPlot)
  })
  
  
  output$simHeatmap <- renderPlot({
    getMetricArgs <- function(curr_metric) {
      if(curr_metric == "get_MoC") return(NULL)
      # if(curr_metric == "get_bin_JaccardIndex") return(list("binSize" = input$bin_size))
      # if(curr_metric == "get_boundaries_JaccardIndex") return(list("tolRad" = as.numeric(input$tol_rad), "matchFor"= input$simMatch))
      # if(curr_metric == "get_ratioMatchingTADs") return(list("coverMatchRatioThresh" = input$ratioCovMatch, "matchFor"=input$simMatch))
      if(curr_metric == "get_bin_JaccardIndex") return(list("binSize" = input$heatmap_bin_size))
      if(curr_metric == "get_boundaries_JaccardIndex") return(list("tolRad" = as.numeric(input$heatmap_tol_rad), "matchFor"= input$heatmap_simMatch))
      if(curr_metric == "get_ratioMatchingTADs") return(list("coverMatchRatioThresh" = input$heatmap_ratioCovMatch, "matchFor"=input$heatmap_simMatch))
      
      if(curr_metric == "get_variationInformation") return(NULL)
    } 
    if(input$demo) {
      demo_files <- get_demoFiles()
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
                         # c(list(all_domains_dt, metric=input$simMetric, lowColor=heatmapLowCol, highColor=heatmapHighCol, nCpu=input$nCpu),
                         #   getMetricArgs(input$simMetric)))
                         c(list(all_domains_dt, metric=input$heatmap_simMetric, lowColor=heatmapLowCol, highColor=heatmapHighCol, nCpu=input$nCpu),
                           getMetricArgs(input$heatmap_simMetric)))
    simPlot <- metricOut[[1]]
    
    return(simPlot)
  })
  
  
  ###
  ##### Functions to download the output PNG #####
  ###
  output$downloadSimHeatmapPNG <- downloadHandler(
    # filename=paste0(getSimMetric(), "_simHeatmap.png"),
    filename=paste0(getSimMetric_heatmap(), "_simHeatmap.png"),
    content = function(file){
      ggsave(filename = file, plot = simHeatmap(), width=widthGG, height=heightGG)
    },
    contentType="image"
  )
  ###
  ##### Functions to download the output SVG #####
  ###
  output$downloadSimHeatmapSVG <- downloadHandler(
    # filename=paste0(getSimMetric(), "_simHeatmap.png"),
    filename=paste0(getSimMetric_heatmap(), "_simHeatmap.svg"),
    content = function(file){
      ggsave(filename = file, plot = simHeatmap(), width=widthGG, height=heightGG)
    },
    contentType="image"
  )
  
  
  ##################################-----------------------------------------------------------------------------------
  ### Enrichment struct. prot. -1 dataset ###-----------------------------------------------------------------------------------
  ##################################-----------------------------------------------------------------------------------
  
  get_structProtFiles <- function() {
    if(input$demo_CTCF) {
      ctcf_file <- demoFile_ctcf
    } else {
      ctcf_file <- input$ctcfFile$datapath
    }
    if(input$demo_RAD21) {
      rad21_file <- demoFile_rad21
    } else {
      rad21_file <- input$rad21File$datapath
    }
    if(input$demo_SMC3) {
      smc3_file <- demoFile_smc3
    } else {
      smc3_file <- input$smc3File$datapath
    }
    print(paste0("ctcf_file=", ctcf_file))
    print(paste0("rad21_file=", rad21_file))
    print(paste0("smc3_file=", smc3_file))
    return(
      list(CTCF=ctcf_file, RAD21=rad21_file, SMC3=smc3_file)
    )
  }
  
  
  structPlot <- function(){
    
    
    if(input$demo) {
      demo_files <- get_demoFiles()
      filenames <- basename(demo_files)
      domain_file <- demo_files[which( basename(demo_files) == input$datasetForStructProt)]
    } else {
      filenames <- input$allDomainFiles$name
      domain_file <- input$allDomainFiles$datapath[which( input$allDomainFiles$name == input$datasetForStructProt)]
    }
    
    
    print("hello2\n")
    
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
                             # ProteinPeaks_Folder=proteinPeaks_Folder,
                             ProteinPeaks_Files=get_structProtFiles(),
                             plot_tit = plotTit)
    # get_structEnrichmentPlot()
  }
  
  
  
  
  ###
  ##### Functions to call struct prot plot  - 1 dataset
  ###
  output$structEnrichment <- renderPlot({
    
    print(input$structProtSelect)
    print("hello\n")
    return(structPlot())
  })
  
  ###
  ##### Functions to call struct prot plot  - multi dataset
  ###
  output$multi_structEnrichment <- renderPlot({
    
    print(input$multi_structProtSelect)
    print("hello multi\n")
    return(multi_structPlot())
  })
  
  ###
  ##### Functions to download the output PNG #####- struct prot enrichment - 1 dataset
  ###
  output$downloadStructProtEnrichmentPNG <- downloadHandler(
    filename=paste0("structProtEnrichmentPlot.png"),
    content = function(file){
      # ggsave(filename = file, plot = get_structEnrichmentPlot(), width=widthGG, height=heightGG)
      # png(filename = file, plot = get_structEnrichmentPlot(), width=widthGG, height=heightGG)
      png(filename = file, width=400, height=400)
      # get_structEnrichmentPlot()
      structPlot()
      dev.off()
    },
    contentType="image"
  )
  ###
  ##### Functions to download the output PNG #####- struct prot enrichment - 1 dataset
  ###
  output$downloadStructProtEnrichmentSVG <- downloadHandler(
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
  ##### Functions to download the output SVG #####- struct prot enrichment - multi datasets
  ###
  output$downloadMultiStructProtEnrichmentSVG <- downloadHandler(
    filename=paste0("multi_structProtEnrichmentPlot.svg"),
    content = function(file){
      # ggsave(filename = file, plot = get_structEnrichmentPlot(), width=widthGG, height=heightGG)
      # png(filename = file, plot = get_structEnrichmentPlot(), width=widthGG, height=heightGG)
      svg(filename = file, width=7, height=7)
      # get_structEnrichmentPlot()
      multi_structPlot()
      dev.off()
    },
    contentType="image"
  )
  
  ###
  ##### Functions to download the output PNG #####- struct prot enrichment - multi datasets
  ###
  output$downloadMultiStructProtEnrichmentPNG <- downloadHandler(
    filename=paste0("multi_structProtEnrichmentPlot.png"),
    content = function(file){
      # ggsave(filename = file, plot = get_structEnrichmentPlot(), width=widthGG, height=heightGG)
      # png(filename = file, plot = get_structEnrichmentPlot(), width=widthGG, height=heightGG)
      png(filename = file, width=400, height=400)
      # get_structEnrichmentPlot()
      multi_structPlot()
      dev.off()
    },
    contentType="image"
  )
  
  ###
  ### Show current datasets for selection
  ###
  
  output$showCurrentDatasets <- renderUI({
    
    if(input$demo) {
      demo_files <- get_demoFiles()
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
  ### Select current datasets for selection - for struct. plot 1 dataset; among the current data -> select which to plot
  ###
  output$secondSelection <- renderUI({
    if(input$demo) {
      demo_files <- get_demoFiles()
      filenames <- basename(demo_files)
    } else {
      filenames <- input$allDomainFiles$name
    }
    selectInput("datasetForStructProt", "Dataset:", choices = filenames)
  })
  
  ###
  ### Select current datasets for selection - for struct. plot multi datasets; among the current data -> select which to plot
  ###
  output$multi_secondSelection <- renderUI({
    if(input$demo) {
      demo_files <- get_demoFiles()
      filenames <- basename(demo_files)
    } else {
      filenames <- input$allDomainFiles$name
    }
    checkboxGroupInput("multi_datasetForStructProt", "Dataset:", choices = filenames, selected="all")
  })
  
  
  
  ##################################-----------------------------------------------------------------------------------
  ### Enrichment struct. prot. - multiple datasets comparison ###-----------------------------------------------------------------------------------
  ##################################-----------------------------------------------------------------------------------
  
  
  multi_structPlot <- function(){
    
    
    if(input$demo) {
      demo_files <- get_demoFiles()
      filenames <- basename(demo_files)
      all_domain_files <- demo_files[which( basename(demo_files) %in% input$multi_datasetForStructProt)]
    } else {
      filenames <- input$allDomainFiles$name
      all_domain_files <- input$allDomainFiles$datapath[which( input$allDomainFiles$name %in% input$multi_datasetForStructProt)]
    }
    
    print("get_structProtFiles")
    print(get_structProtFiles())
    print("hello2\n")
    
    # plot(x=c(1:10), y=c(1:10))
    
    print(paste0("input$multi_datasetForStructProt = ", input$multi_datasetForStructProt))
    print(paste0("input$allDomainFiles$name = ", input$allDomainFiles$name))
    
    
    print(paste0("all_domain_files = ", all_domain_files))
    print("hello3")
    
    fieldSep <- input$fieldSep
    fileHeader <- input$fileHeader
    
    if(length(all_domain_files) < 1) {
      return(HTML("! At least one file should be provided ! "))
    }
    
    # print(paste0("domain_file = ", domain_file))
    
    all_domain_dt <- lapply(all_domain_files, function(xfile) {
      read.table(xfile,
                 header=fileHeader, 
                 stringsAsFactors = FALSE,  sep=fieldSep)#col.names=c("chromo", "start", "end"),
    })
    stopifnot(length(all_domain_dt) == length(all_domain_files))
    ds_names <- paste0(seq_along(all_domain_files), "_", basename(all_domain_files))
    names(all_domain_dt) <- ds_names
    
    defaultTit <- "FC and pval enrichment comparison"
    
    plotTit <- ifelse(input$multi_enrichmentPlotTit == "" | 
                        length(input$multi_enrichmentPlotTit) == 0 |is.null(input$multi_enrichmentPlotTit) ,
                      defaultTit, input$multi_enrichmentPlotTit)
    # 
    plot_regexsub <- ifelse(input$multi_enrichmentRegexSub == "" |
                              length(input$multi_enrichmentRegexSub) == 0 |is.null(input$multi_enrichmentRegexSub) ,
                            "", input$multi_enrichmentRegexSub)
    
    
    
    print(paste0("plotTit = ", plotTit))
    print(paste0("is.null plotTit = ", is.null(plotTit)))
    print(paste0("empty plotTit = ", plotTit==""))
    
    if(input$demo) {
      hic_res <- hicDemoResKb
    } else {
      hic_res <- input$multi_structProt_resolutionKb
    }
    
    stopifnot(is.numeric(hic_res))
    stopifnot(is.numeric(input$multi_structProt_resStep))
    
    get_multiStructEnrichmentPlot(all_tad_dt=all_domain_dt, 
                                  proteins = input$multi_structProtSelect,
                                  resolution_kb=hic_res, 
                                  res_step=input$structProt_resStep,
                                  # ProteinPeaks_Folder=proteinPeaks_Folder,
                                  ProteinPeaks_Files=get_structProtFiles(),
                                  plot_tit = plotTit,
                                  regexsub = plot_regexsub)
    # get_structEnrichmentPlot()
  }
  
  
  
  
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

