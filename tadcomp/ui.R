#######################################################################################
############################# CARPnTrees web application ############################## 
############# building consensus tree from ancient and recent phylogenies #############  
# Marie Zufferey - UNIL - December 2015 ###############################################
# License: Open Source AGPL v3 ########################################################
#######################################################################################

#####################################
##### User interface - web app  #####
#####################################
#minimal code 
#headerPanel(), sidebarPanel(), mainPanel()

library(shiny)
library(shinythemes)
library(shinyBS)
library(shinyjs)
library(graphics)

possible_metrics <- list("Measure of Concordance (MoC)" = "get_MoC", 
                         "Jaccard Index (bins)" = "get_bin_JaccardIndex",
                         "Jaccard Index (boundaries)" = "get_boundaries_JaccardIndex", 
                         "Ratio matching TADs" = "get_ratioMatchingTADs",
                         "Variation of information" = "get_variationInformation")

possible_match <- list("Symmetric" = "all",
                      "Set1 as ref." = "set1",
                      "Set2 as ref." = "set2"
                      )

possible_structProt <- list("CTCF" = "CTCF", "RAD21" = "RAD21", "SMC3" = "SMC3")

shinyUI(fluidPage(    
  
  # => to place the notification in the middle of the page
  tags$head(
    tags$style(
      HTML(".shiny-notification {
             position:fixed;
             top: calc(50%);
             left: calc(50%);
             }
             "
      )
    )
  ),
  
  theme = shinytheme("cosmo"),
  HTML('<title>Compare TAD lists</title>
       <h1 style="color:#2780E3;font-size:400%;" align="center"><b>Compare TAD partitions </b></h1>
       <h4 align="center"><i>comparison of TAD caller outputs </i></h4>
       <hr>'),
  sidebarLayout(  
    sidebarPanel(
      
      checkboxInput('demo', 'Example with demo data', TRUE),
      
      tags$hr(),
      
      selectInput("simMetric", 
                  label = "Choose a comparison metric",
                  choices = possible_metrics,
                  selected = "get_boundaries_JaccardIndex"),
      
      sliderInput(inputId='nCpu', label="# available Cpu",
                    min = 1, max = 100, value=1),
        
        conditionalPanel(
          condition="input.simMetric=='get_bin_JaccardIndex'",
          sliderInput(inputId='bin_size', label="Hi-C matrix bin size",
                      min = 5000, max = 500000, value=25000, step=5000)
        ),
        
        
        conditionalPanel(
          condition="input.simMetric=='get_boundaries_JaccardIndex'",
          sliderInput(inputId="tol_rad", label="Matching tolerance radius",
                      min = 5000, max = 500000, value=50000, step=5000)
        ),
        conditionalPanel(
          condition="input.simMetric=='get_boundaries_JaccardIndex' | input.simMetric == 'get_ratioMatchingTADs'",
          selectInput("simMatch", 
                      label = "Direction of matching",
                      choices = possible_match,
                      selected = "all")
        ),
        conditionalPanel(
          condition="input.simMetric=='get_ratioMatchingTADs'",
          sliderInput("ratioCovMatch", 
                      label = "Ratio cov. needed for matching",
                      min = 0, max = 1, value = 0.8, step=0.1)
        ),
      conditionalPanel(
        condition="input.demo == false",
        tags$hr(),
        fileInput("allDomainFiles", label=("Files with TAD coordinates"), multiple=TRUE, placeholder="--- TAD list in bed-like format ---"),
        tags$hr()
      ),
      conditionalPanel(
        condition="input.demo == false",
        HTML('<u>Input format:</u>'),
        checkboxInput('fileHeader', 'File with header', FALSE),
        radioButtons('fieldSep', 'Separator',
                     c(Comma=',',
                       Semicolon=';',
                       Tab='\t'),
                     '\t'),
        useShinyjs(),                                                  # Include shinyjs in the UI
        actionButton("reset_button", "Reset your choices")
      )
    ),
    mainPanel(
      tabsetPanel(
        ###### TAB WITH OVERVIEW OF LOADED DATA #####
        tabPanel("Loaded data",
                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                  tags$div("Please wait - Datasets loading...",id="loading_message")),
                 tableOutput("loaded_data")
        ),
        ###### TAB WITH CURRENT SELECTED PARAMETERS #####
        tabPanel("Current parameters",
                 htmlOutput("showCurrentOptions")
        ),
        
        ###### TAB WITH CURRENT SELECTED PARAMETERS #####
        tabPanel("Dataset overview",
                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                  tags$div("Please wait - Analysis running...",id="running_message")),
                 plotOutput("descFeatures"),
                 HTML(    '<table style="width:100%" align="center">
                          <tr>
                               <td>
                              <style>a#outputDwld{width:200px;height:30px;padding-top:5px;}</style>'),
                 
                 conditionalPanel(condition="! $('html').hasClass('shiny-busy')",
                                  downloadButton('downloadDescTable', 'Download # TADs table (.txt)')),
                 HTML(    '</td>
                               <td>
                              <style>a#outputDwld{width:200px;height:30px;padding-top:5px;}</style>'),
                 
                 conditionalPanel(condition="! $('html').hasClass('shiny-busy')",
                                  downloadButton('downloadDescPlot', 'Download # TADs and size plots (.png)')),
                 HTML(    '</td>
                          </tr>
                          </table>')
        ),
        
        ###### TAB WITH THE SIMILARITY TABLE #####
        tabPanel("Similarity table",
                 
                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                  tags$div("Please wait - Analysis running...",id="running_message")),
                 
                 tableOutput("simTable"),
                 
                 HTML(    '<table style="width:100%">
                          <tr>
                               <td>
                              <style>a#outputDwld{width:200px;height:30px;padding-top:5px;}</style>'),
                 
                 conditionalPanel(condition="! $('html').hasClass('shiny-busy')",
                                  downloadButton('downloadSimTable', 'Download similarity table (.txt)')),
                 HTML(    '
                          </tr>
                               </td>
                          </table>')
                 
        ),
        ###### TAB WITH THE SIMILARITY HEATMAP #####
        tabPanel("Similarity heatmap",
                 
                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                  tags$div("Please wait - Analysis running...",id="running_message")),
                 plotOutput("simHeatmap"),
                 HTML(    '<table style="width:100%">
                          <tr>
                               <td>
                              <style>a#outputDwld{width:200px;height:30px;padding-top:5px;}</style>'),
                 conditionalPanel(condition="! $('html').hasClass('shiny-busy')",
                                  downloadButton('downloadSimHeatmap', 'Download similarity heatmap (.png)')),
                 HTML(    '
                          </tr>
                               </td>
                          </table>')
        ),
        ###### TAB WITH DETAILS ABOUT METRICS #####
        tabPanel("Details about the metrics",
                 htmlOutput("metricDetails")
            ),
        
        
        ###### TAB WITH STRUCT PROT ENRICHMENT #####
        tabPanel("Struct. prot. enrichment",
                 
                 # conditionalPanel(
                   # condition="input.structProt=='' | input.structProt == ''",
                 checkboxGroupInput("structProtSelect", 
                               label = "Struct. prot.",
                               choices = possible_structProt,
                               selected = "all"),
                 # ),
                 
                 conditionalPanel(condition="input.demo == false",
                 sliderInput(inputId='structProt_resolutionKb', label="Hi-C data resolution (kb)",
                             min = 1, max = 100, value=25)),
                 
                 
                 sliderInput(inputId='structProt_resStep', label="Peak agg. resolution step (default: 5)",
                             min = 1, max = 10, value=5),
                 
                 
                 textInput("enrichmentPlotTit", label = h3(""), value = "Enter plot title..."),
                 
                 
                 htmlOutput("showCurrentDatasets"),

                 uiOutput("secondSelection"),
                 
                 htmlOutput("showSelectedDatasets"),
                 
                 conditionalPanel(condition="input.structProtSelect == ''",
                                  tags$div("No protein selected",id="running_message")),
                 
                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                  tags$div("Please wait - Analysis running...",id="running_message")),
                 
                 conditionalPanel(condition="input.structProtSelect != ''",
                                  plotOutput("structEnrichment")),
                 
                                  
                 
                 
   # conditionalPanel(condition="input.structProtSelect != ''",
   #                   plotOutput("structEnrichment")),
                 
                 
                 
                 HTML(    '<table style="width:100%">
                          <tr>
                               <td>
                              <style>a#outputDwld{width:200px;height:30px;padding-top:5px;}</style>'),
                ### can be busy because redo everything...
                 # conditionalPanel(condition="! $('html').hasClass('shiny-busy')", 
                 #                  downloadButton('downloadStructProtEnrichment', 'Download enrichment plot (.png)')),
   conditionalPanel(condition="input.structProtSelect != ''",
                    downloadButton('downloadStructProtEnrichment', 'Download enrichment plot (.png)')),
   
   
                 HTML(    '
                          </tr>
                               </td>
                          </table>')
        ),
        
        
        
       ###### TAB WITH HELP #####
       tabPanel("Help",
                htmlOutput("help")
       ),
       ###### TAB WITH CONTACT #####
        tabPanel("Contact",
                 htmlOutput("contact")
        )
      ) # close tabsetPanel
    ) # close mainPanel
  ), # close sidebar layout
                 

  # ADD THE FOOTER WITH PICTURES FOR THE LOGO AND THE LICENSE TEXT
HTML('<hr><footer>
            <a href="http://www.unil.ch" target="_blank">
              <img src="images/logo_UNIL.png" alt="UNIL logo" width="5%" 
                  style="vertical-align:bottom;">
              </a><div style="text-align:center;margin-top:-21px">
              <a align=center href="https://opensource.org/licenses/AGPL-3.0" style="font-size:70%;"  target="_blank"> 
              Open Source AGPL v3 License </a></div>

<div style="text-align:right;margin-top:-21px">
        <a href="http://ciriellolab.org" target="_blank">

              <img src="images/logo_CSO.png" alt="CSO logo" width="5%" 
                  style="vertical-align:bottom;"></a><div style="text-align:center;margin-top:-21px">
</div>

              </p>
              </footer>')
        )  # close fluidPage
  ) # close shinyUI
