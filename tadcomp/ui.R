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

possible_demoBinSize <- list(
  "10kb" = "10kb",
  "50kb" = "50kb",
  "100kb" = "100kb"
)

possible_demoNorm <- list(
  "ICE" = "ICE",
  "LGF" = "LGF"
)

possible_demoSub <- list(
  "1" = "1",
  "0.1" = "0.1",
  "0.01" = "0.01",
  "0.001" = "0.001"
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
      
      h3("Hi-C data"),
      
      checkboxInput('demo', 'Example with demo data', TRUE),
      
      conditionalPanel(
        condition="input.demo==true",
        selectInput("demoNorm", 
                         label = "Demo data normalization",
                         choices = possible_demoNorm,
                         selected = "50kb")
      ),
      conditionalPanel(
        condition="input.demo==true",
        selectInput("demoBinSize", 
                    label = "Demo data bin size",
                    choices = possible_demoBinSize,
                    selected = "10kb")
      ),
      conditionalPanel(
        condition="input.demo==true && input.demoNorm=='ICE' && input.demoBinSize == '50kb'",
        selectInput("demoSub", 
                    label = "Demo data subsampling (1=no sub.)",
                    choices = possible_demoSub,
                    selected = "1")
      ),
      
      conditionalPanel(
        condition="input.demo == false",
        tags$hr(),
        fileInput("allDomainFiles", 
                  label=("Files with TAD coordinates"), multiple=TRUE, placeholder="--- TAD list in bed-like format ---"),
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
      )  ,
      tags$hr(),
      
      h3("ChIP-seq data"),
      
      checkboxInput('demo_CTCF', HTML('CTCF demo data<br/>(chr6 - GM12878)'), TRUE),
      checkboxInput('demo_RAD21', HTML('RAD21 demo data<br/>(chr6 - GM12878)'), TRUE),
      checkboxInput('demo_SMC3', HTML('SMC3 demo data<br/>(chr6 - GM12878)'), TRUE),
      
      
      conditionalPanel(
        condition="input.demo_CTCF == false",
        tags$hr(),
        fileInput("ctcfFile", label=("Files with CTCF peaks"), multiple=FALSE, placeholder="--- CTCF peaks (see Help for format) ---"),
        tags$hr()
      ),
      
      conditionalPanel(
        condition="input.demo_RAD21 == false",
        tags$hr(),
        fileInput("rad21File", label=("Files with RAD21 peaks"), multiple=FALSE, placeholder="--- RAD21 peaks (see Help for format) ---"),
        tags$hr()
      ),
      
      conditionalPanel(
        condition="input.demo_SMC3 == false",
        tags$hr(),
        fileInput("smc3File", label=("Files with SMC3 peaks"), multiple=FALSE, placeholder="--- SMC3 peaks (see Help for format) ---"),
        tags$hr()
      ),
      
      
      tags$hr(),
      
      ### update 20.02.21 -> now I separate heatmap and table!
      
      # selectInput("simMetric", 
      #             label = "Choose a comparison metric",
      #             choices = possible_metrics,
      #             selected = "get_boundaries_JaccardIndex"),
      
      h3("Ressources"),
      
      # 
      sliderInput(inputId='nCpu', label="# available Cpu",
                    min = 1, max = 100, value=1),
      #   
      #   conditionalPanel(
      #     condition="input.simMetric=='get_bin_JaccardIndex'",
      #     sliderInput(inputId='bin_size', label="Hi-C matrix bin size",
      #                 min = 5000, max = 500000, value=25000, step=5000)
      #   ),
      #   
      #   
      #   conditionalPanel(
      #     condition="input.simMetric=='get_boundaries_JaccardIndex'",
      #     sliderInput(inputId="tol_rad", label="Matching tolerance radius",
      #                 min = 5000, max = 500000, value=50000, step=5000)
      #   ),
      #   conditionalPanel(
      #     condition="input.simMetric=='get_boundaries_JaccardIndex' | input.simMetric == 'get_ratioMatchingTADs'",
      #     selectInput("simMatch", 
      #                 label = "Direction of matching",
      #                 choices = possible_match,
      #                 selected = "all")
      #   ),
      #   conditionalPanel(
      #     condition="input.simMetric=='get_ratioMatchingTADs'",
      #     sliderInput("ratioCovMatch", 
      #                 label = "Ratio cov. needed for matching",
      #                 min = 0, max = 1, value = 0.8, step=0.1)
      #   )
      # 
      
      
      

    ),
    mainPanel(
      tabsetPanel(
        ###### TAB WITH OVERVIEW OF LOADED DATA #####
        tabPanel("Loaded data",
                 
                 h2("Hi-C data"),
                 
                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                  # tags$div("Please wait - Datasets loading...",id="loading_message")
                                  tags$span(style="color:blue", tags$div("Please wait - Datasets loading...",id="loading_message"))
                                  ),
                 
                 conditionalPanel(condition="output.domainFilesUploaded == false & input.demo == false",
                                  # tags$div("Please wait - Analysis running...",id="running_message")
                                  tags$span(style="color:orange",tags$div("No Hi-C data selected",id="nohicdata_message"))
                 ),
                 tableOutput("loaded_data"),
                 
                 
                 h2("ChIP-seq data"),
                 
                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                  # tags$div("Please wait - Datasets loading...",id="loading_message")
                                  tags$span(style="color:blue", tags$div("Please wait - Datasets loading...",id="loading_message"))
                 ),
                 
                 conditionalPanel(condition="output.chipFilesUploaded == false & input.demo_CTCF == false & input.demo_RAD21 == false & input.demo_SMC3 == false",
                                  # tags$div("Please wait - Analysis running...",id="running_message")
                                  tags$span(style="color:orange",tags$div("No ChIP-seq data selected",id="nochipdata_message"))
                 ),
                 tableOutput("loaded_data_chip")
        ),
        ###### TAB WITH CURRENT SELECTED PARAMETERS #####
        tabPanel("Current parameters",
                 htmlOutput("showCurrentOptions")
        ),
        
        ###### TAB WITH CURRENT SELECTED PARAMETERS #####
        tabPanel("Hi-C data overview",
                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                  # tags$div("Please wait - Analysis running...",id="running_message")
                                  tags$span(style="color:blue",tags$div("Please wait - Analysis running...",id="running_message"))
                                  ),
                 
                 
                 conditionalPanel(condition="output.domainFilesUploaded == false & input.demo == false",
                                  # tags$div("Please wait - Analysis running...",id="running_message")
                                  tags$span(style="color:orange",tags$div("No Hi-C data selected",id="nohicdata_message"))
                 ),
                 
                 
                 
                 plotOutput("descFeatures"),
                 HTML(    '<table style="width:100%" align="center">
                          <tr>
                               <td>
                              <style>a#outputDwld{width:200px;height:30px;padding-top:5px;}</style>'),
                 
                 # conditionalPanel(condition="! $('html').hasClass('shiny-busy')",
                 #                  downloadButton('downloadDescTable', 'Download # TADs table (.txt)')),
                 # not need conditional because redo all the analysis for downloading...
                 conditionalPanel(condition="output.domainFilesUploaded == true | input.demo == true",
                                  downloadButton('downloadDescTable', 'Download # TADs table (.txt)')),
                 HTML(    '</td>
                               <td>
                              <style>a#outputDwld{width:200px;height:30px;padding-top:5px;}</style>'),
                 
                 # conditionalPanel(condition="! $('html').hasClass('shiny-busy')",
                 #                  downloadButton('downloadDescPlot', 'Download # TADs and size plots (.png)')),
                 # not need conditional because redo all the analysis for downloading...
                 conditionalPanel(condition="output.domainFilesUploaded == true | input.demo == true",
                                  downloadButton('downloadDescPlot', 'Download # TADs and size plots (.png)')),
                 
                 HTML(    '</td>
                          </tr>
                          </table>')
        ),
        
        ###### TAB WITH THE SIMILARITY TABLE #####
        tabPanel("Similarity table",
                 
                 
                 fluidRow(
                   column(3,
                 
                 
                 selectInput("table_simMetric", 
                             label = "Choose a comparison metric",
                             choices = possible_metrics,
                             selected = "get_boundaries_JaccardIndex"),
                 
                 tags$hr(),
                 conditionalPanel(
                   condition="input.table_simMetric=='get_boundaries_JaccardIndex' | input.table_simMetric == 'get_ratioMatchingTADs'",
                   selectInput("table_simMatch", 
                               label = "Direction of matching",
                               choices = possible_match,
                               selected = "all")
                 )
                 ),
                 column(4,
                 # sliderInput(inputId='table_nCpu', label="# available Cpu",
                 #             min = 1, max = 100, value=1),
                 # tags$hr(),
                 conditionalPanel(
                   condition="input.table_simMetric=='get_bin_JaccardIndex'",
                   sliderInput(inputId='table_bin_size', label="Hi-C matrix bin size",
                               min = 5000, max = 500000, value=25000, step=5000)
                 ),
                 conditionalPanel(
                   condition="input.table_simMetric=='get_boundaries_JaccardIndex'",
                   sliderInput(inputId="table_tol_rad", label="Matching tolerance radius",
                               min = 5000, max = 500000, value=50000, step=5000)
                 ),
                 conditionalPanel(
                   condition="input.table_simMetric=='get_ratioMatchingTADs'",
                   sliderInput("table_ratioCovMatch", 
                               label = "Ratio cov. needed for matching",
                               min = 0, max = 1, value = 0.8, step=0.1)
                 )
                 ), align="left"),
                 
                 
                 
                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                  tags$span(style="color:blue",tags$div("Please wait - Analysis running...",id="running_message"))),
                                            # tags$div("Please wait - Analysis running...",id="running_message")),
                 
                 conditionalPanel(condition="output.domainFilesUploaded == false & input.demo == false",
                                  # tags$div("Please wait - Analysis running...",id="running_message")
                                  tags$span(style="color:orange",tags$div("No Hi-C data selected",id="nohicdata_message"))
                 ),
                 
                 tableOutput("simTable"),
                 
                 HTML(    '<table style="width:100%">
                          <tr>
                               <td>
                              <style>a#outputDwld{width:200px;height:30px;padding-top:5px;}</style>'),
                 
                 # conditionalPanel(condition="! $('html').hasClass('shiny-busy')",
                 #                  downloadButton('downloadSimTable', 'Download similarity table (.txt)')),
                 # not need conditional because redo all the analysis for downloading...
                 conditionalPanel(condition="output.domainFilesUploaded == true | input.demo == true ",
                                  downloadButton('downloadSimTable', 'Download similarity table (.txt)')
                 ),
                 
                 HTML(    '
                          </tr>
                               </td>
                          </table>')
                 
        ),
        ###### TAB WITH THE SIMILARITY HEATMAP #####
        tabPanel("Similarity heatmap",
                 
                 
                 fluidRow(
                   column(3,
                          
                          
                          selectInput("heatmap_simMetric", 
                                      label = "Choose a comparison metric",
                                      choices = possible_metrics,
                                      selected = "get_boundaries_JaccardIndex"),
                          
                          tags$hr(),
                          conditionalPanel(
                            condition="input.heatmap_simMetric=='get_boundaries_JaccardIndex' | input.heatmap_simMetric == 'get_ratioMatchingTADs'",
                            selectInput("heatmap_simMatch", 
                                        label = "Direction of matching",
                                        choices = possible_match,
                                        selected = "all")
                          )
                   ),
                   column(4,
                          # sliderInput(inputId='heatmap_nCpu', label="# available Cpu",
                          #             min = 1, max = 100, value=1),
                          # tags$hr(),
                          conditionalPanel(
                            condition="input.heatmap_simMetric=='get_bin_JaccardIndex'",
                            sliderInput(inputId='heatmap_bin_size', label="Hi-C matrix bin size",
                                        min = 5000, max = 500000, value=25000, step=5000)
                          ),
                          conditionalPanel(
                            condition="input.heatmap_simMetric=='get_boundaries_JaccardIndex'",
                            sliderInput(inputId="heatmap_tol_rad", label="Matching tolerance radius",
                                        min = 5000, max = 500000, value=50000, step=5000)
                          ),
                          conditionalPanel(
                            condition="input.heatmap_simMetric=='get_ratioMatchingTADs'",
                            sliderInput("heatmap_ratioCovMatch", 
                                        label = "Ratio cov. needed for matching",
                                        min = 0, max = 1, value = 0.8, step=0.1)
                          )
                   ), align="left"),
                 
                 
                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                  # tags$div("Please wait - Analysis running...",id="running_message")
                                  tags$span(style="color:blue",tags$div("Please wait - Analysis running...",id="running_message"))
                                  ),
                 
                 
                 conditionalPanel(condition="output.domainFilesUploaded == false & input.demo == false",
                                  # tags$div("Please wait - Analysis running...",id="running_message")
                                  tags$span(style="color:orange",tags$div("No Hi-C data selected",id="nohicdata_message"))
                 ),
                 
                 plotOutput("simHeatmap"),
                 HTML(    '<table style="width:100%">
                          <tr>
                               <td>
                              <style>a#outputDwld{width:200px;height:30px;padding-top:5px;}</style>'),
                 # conditionalPanel(condition="! $('html').hasClass('shiny-busy')",
                 #                  downloadButton('downloadSimHeatmapPNG', 'Download similarity heatmap (.png)')),
                 # # not needed conditional because redo
                 
                 
                 conditionalPanel(condition="output.domainFilesUploaded == true | input.demo == true",
                                  downloadButton('downloadSimHeatmapPNG', 'Download similarity heatmap (.png)')
                 ),
                 
        
                 HTML(    '
                          </td>
                          <td>
                          <style>a#outputDwld{width:200px;height:30px;padding-top:5px;}</style>'),
                 
                 
                 conditionalPanel(condition="output.domainFilesUploaded == true | input.demo == true",
                                  downloadButton('downloadSimHeatmapSVG', 'Download similarity heatmap (.svg)')
                 ),
                 
                          
                          HTML('    
                          </td>
                          </tr>
                          </table>')
        ),
        ###### TAB WITH DETAILS ABOUT METRICS #####
        tabPanel("Details about the metrics",
                 htmlOutput("metricDetails")
            ),
        
        
        ###### TAB WITH STRUCT PROT ENRICHMENT - 1 dataset #####
        tabPanel("Struct. prot. enrichment\n(1 dataset)",
                 
                 
                 fluidRow(
                       # conditionalPanel(
                         # condition="input.structProt=='' | input.structProt == ''",
                       column(width=3, checkboxGroupInput("structProtSelect", 
                                     label = "Struct. prot.",
                                     choices = possible_structProt,
                                     selected = "all")),
                       # ),
                       
                       column(width=3, sliderInput(inputId='structProt_resStep', label="Peak agg. resolution [kb] (default: 5)",
                                   min = 1, max = 10, value=5)),
                       #  (average number of peaks in 5-kb intervals)
                       
                       column(width=3,conditionalPanel(condition="input.demo == false",
                                        sliderInput(inputId='structProt_resolutionKb', label="Hi-C data resolution [kb]",
                                                    min = 1, max = 100, value=25))),
                       column(12), align="left"
                 ),
                 
                 
                 
                 
                 # htmlOutput("showCurrentDatasets"),# >> FOR DEBUG
                 
                 fluidRow(

                   column(width=4, uiOutput("secondSelection")),
                   # label = h3("Plot title")
                          column(width=4,textInput("enrichmentPlotTit", label = "Plot title", value = "Enter plot title...")),
                   column(12), align="left"
                 ),
                 
                 # htmlOutput("showSelectedDatasets"),  # >> FOR DEBUG
                 
                 conditionalPanel(condition="input.structProtSelect == ''",
                                  # tags$div("No protein selected",id="noprot_message")
                                  tags$span(style="color:orange", tags$div("No protein selected",id="noprot_message"))
                                  ),
                 
                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                  # tags$div("Please wait - Analysis running...",id="running_message")
                                  tags$span(style="color:blue", tags$div("Please wait - Analysis running...",id="running_message"))
                                  ),
                 
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
                    downloadButton('downloadStructProtEnrichmentPNG', 'Download enrichment plot (.png)')),
   
   HTML(    '</td><td>  <style>a#outputDwld{width:200px;height:30px;padding-top:5px;}</style>'),
   conditionalPanel(condition="input.structProtSelect != ''",
                    downloadButton('downloadStructProtEnrichmentSVG', 'Download enrichment plot (.svg)')),
   
                 HTML(    '
                          </td>
                               </tr>
                          </table>')
        ),
        
   
   
   
   ###### TAB WITH STRUCT PROT ENRICHMENT comparison FC multiple datasets #####
   tabPanel("Struct. prot. enrichment\n(multiple datasets)",
            # https://shiny.rstudio.com/articles/layout-guide.html
            
            fluidRow(
              # conditionalPanel(
              # condition="input.structProt=='' | input.structProt == ''",
              column(width=3, checkboxGroupInput("multi_structProtSelect", 
                                                 label = "Struct. prot.",
                                                 choices = possible_structProt,
                                                 selected = "all")),
              # ),
              
              column(width=3, sliderInput(inputId='multi_structProt_resStep', label="Peak agg. resolution [kb] (default: 5)",
                                          min = 1, max = 10, value=5)),
              #  (average number of peaks in 5-kb intervals)
              
              column(width=3,conditionalPanel(condition="input.demo == false",
                                              sliderInput(inputId='multi_structProt_resolutionKb', label="Hi-C data resolution [kb]",
                                                          min = 1, max = 100, value=25))),
              # column(12),
              align="left"
            ),
            
            
            
            
            # htmlOutput("showCurrentDatasets"),# >> FOR DEBUG
            
            fluidRow(
              
              column(width=4, uiOutput("multi_secondSelection")),
              # label = h3("Plot title")
              column(width=4,textInput("multi_enrichmentPlotTit", label = "Plot title", value = "Comparison struct. prot. enrichment")),
              
              column(width=4,textInput("multi_enrichmentRegexSub", label = "Sub. regex dataset names", value = "_ICE_final_domains.txt")),
              # column(12), 
              align="left"
            ),
            
            # htmlOutput("showSelectedDatasets"),  # >> FOR DEBUG
            
            conditionalPanel(condition="input.multi_structProtSelect == ''",
                             # tags$div("No protein selected",id="running_message")),
                             tags$span(style="color:orange", tags$div("No protein selected",id="noprot_message"))
                                       ),
            
            conditionalPanel(condition="input.multi_datasetForStructProt == ''",
                             # tags$div("No protein selected",id="running_message")),
                             tags$span(style="color:orange", tags$div("No Hi-C dataset selected",id="nods_message"))
            ),
            
            
            conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                             # tags$div("Please wait - Analysis running...",id="running_message")),
                             tags$span(style="color:blue",  tags$div("Please wait - Analysis running...",id="running_message"))
            ),
            
            conditionalPanel(condition="input.multi_structProtSelect != '' & input.multi_datasetForStructProt != ''",
                             plotOutput("multi_structEnrichment")),
            
            
            
            
            HTML(    '<table style="width:100%">
                          <tr>
                               <td>
                              <style>a#outputDwld{width:200px;height:30px;padding-top:5px;}</style>'),
            ### can be busy because redo everything...
            # conditionalPanel(condition="! $('html').hasClass('shiny-busy')", 
            #                  downloadButton('downloadStructProtEnrichment', 'Download enrichment plot (.png)')),
            conditionalPanel(condition="input.multi_structProtSelect != '' & input.multi_datasetForStructProt != ''",
                             downloadButton('downloadMultiStructProtEnrichmentPNG', 'Download enrichment plot (.png)')),
            HTML('</td><td>   <style>a#outputDwld{width:200px;height:30px;padding-top:5px;}</style>'),
            conditionalPanel(condition="input.multi_structProtSelect != '' & input.multi_datasetForStructProt != ''",
                             downloadButton('downloadMultiStructProtEnrichmentSVG', 'Download enrichment plot (.svg)')),
            HTML(    '
                          </td>
                               </tr>
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
