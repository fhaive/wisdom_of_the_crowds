suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(InteractiveComplexHeatmap))
library(data.table)
library(shiny)
library(shinyBS)
library(shinyjs)
library(shinyWidgets)
library(DT)
library(ComplexHeatmap)
library(circlize)
library(InteractiveComplexHeatmap)
library(shinydashboard)
library(kableExtra)
library(plotly) 
library(ggplot2)
library(reshape2)
library(dplyr) 
library(magrittr)



ui <- fluidPage(
  
  
  useShinyjs(), 
  #to set modal window size
   tags$head(
     tags$script(HTML("
       function setModalWidth(modalId, width) {
         $('#' + modalId + ' .modal-dialog').css('width', width);
       }
     "))
   ),
  
  # to capture modal close event
  tags$script(HTML("
    $(document).on('hidden.bs.modal', function (e) {
      if ($(e.target).attr('id') == 'shiny-modal') {
        Shiny.setInputValue('modal_closed', true, {priority: 'event'});
      }
    });
    ")), 
  
  tabsetPanel(id = "tabs",
    tabPanel("Preloaded data",  
             br(), hr(),
             
             fluidRow(
               div(style="display: inline-block; margin-left:10px;margin-bottom: -30px",
                   column(12,  align ="left",
                          h5("Expert Answers")),
                   div( style="margin-bottom: 45px"),
                   )
             ),
             
             br(),hr(),
             fluidRow(
               div(style="display: inline-block; margin-left:10px;margin-bottom: -30px",
                   column(12,  align ="left",
                          h5("Hazard Concern Scores")),
                   div( style="margin-bottom: 45px"),
               )
             ),
             
             br(), hr(),
             fluidRow(
               div(style="display: inline-block; margin-left:10px;margin-bottom: -30px",
                   column(12,  align ="left",
                          h5("Gene Expression")),
                   div( style="margin-bottom: 45px"),
               )
             ),
             
             br(),hr(),
             fluidRow(
               div(style="display: inline-block; margin-left:10px;margin-bottom: -30px",
                   column(12,  align ="left",
                          h5("Gene - Material Properties")),
                   div( style="margin-bottom: 45px"),
               ) 
             )  
    ),
    
    
    tabPanel("Tables",  
             fluidRow(br()),
             fluidRow(column(4, align="left",
                             selectizeInput("showtables", "Type of Table", choices = c("Hazard", "Genes", "Properties"), 
                                            selected = NULL, multiple = TRUE,
                                            options = list( 
                                              placeholder = "Select table to show",
                                              maxItems=1)
                             )
             )),
             
             
            conditionalPanel(
               condition = " input.showtables == 'Hazard'",
               hr(),
               fluidRow(column(12, align="left", 
                               DT::dataTableOutput("hazard_main_table")
               ))
            ),
            
            conditionalPanel(
              condition = " input.showtables == 'Genes'",
              hr(),
              fluidRow(column(12, align="left", 
                              DT::dataTableOutput("gene_main_table")
              ))
            ),
            
            conditionalPanel(
              condition = " input.showtables == 'Properties'",
              hr(),
              fluidRow(column(12, align="left", 
                              DT::dataTableOutput("properties_main_table")
              ))
            )
    ),
    
    
    tabPanel("Expert", 
             fluidPage(fluidRow(br()),
             fluidRow(column(12, h4("Check the answers of the experts"))),
             br(),
             fluidRow(column(12, align="left", 
                             selectizeInput("expertID", "Select expert", choices = NULL, 
                                            selected = NULL, multiple = FALSE,
                                            options = list( 
                                              placeholder = "Select the expert ID",
                                              maxItems=1)
                             ))),
             fluidRow(column(12, align="left",  
                    div(style = "vertical-align:top; ",
                             plotOutput("heatmap_expert_questions",height="2920px", width="850px") 
                             )
                      )        
             ))
    ),
    
    
    
    tabPanel("Hazard Visualization", 
             fluidRow(column(3, align= "left",
                             br(),br(),br()
                             
             )),
             fluidRow(column(3, align= "left",
                             shinyBS::bsButton("visualize", "Visualize", style="info" , disabled=FALSE)
             )
             ,column(3,  
                     shinyBS::bsButton("mat_categ", "Material Categories" , disabled=FALSE)
             )
             ,column(5, 
                     div(style = "margin-top: -23px"),
                     textInput("show_mat_class", label="Material Class", value="No material category selected") %>% disabled()
             )),
             
             fluidRow(hr(), div(style = "margin-top: -23px"),
                      column(12, align="left",
                             column(7,checkboxInput("hazard_sparse_select", 
                                                    label = "Sparse selection", value = FALSE) 
                                    ))),
            conditionalPanel(condition='input.hazard_sparse_select == true',
              fluidRow(column(4,align ="left", 
                                shinyWidgets::pickerInput(
                                 inputId = "pick_spare_hazard_rows",
                                 label = "Rows:",
                                 choices = NULL,
                                 multiple=TRUE,
                                 options = list("actions-box" = TRUE)
                             )),
                             column(4, align ="left",
                               shinyWidgets::pickerInput(
                                 inputId = "pick_spare_hazard_cols",
                                 label = "Columns:",
                                 choices = NULL,
                                 multiple=TRUE,
                                 options = list("actions-box" = TRUE)
                             )),
                     column(2, align="left", 
                            div(style = "margin-top: 25px"),
                            shinyBS::bsButton("visualize_sparsemodal", "Visualize Sparsed Selection", style="info" , disabled=FALSE)
                          )
                     
             )),
             fluidRow(column(12, 
                             InteractiveComplexHeatmapOutput(heatmap_id="ht_ui",
                                                             height1 = 2450, 
                                                             width1=850, 
                                                             height2 = 500,
                                                             width2=850,
                                                             width3=1000,
                                                             layout= "1|2|3", 
                                                             compact=FALSE,
                                                             action="click", 
                                                            output_ui = htmlOutput("info"), 
                                                             output_ui_float=FALSE)
                             
             )
             )
    ),
    
    tabPanel("Gene Visualization", 
             fluidRow(column(3, align= "left",
                             br(),
                             shinyBS::bsButton("visualize_model", "Visualize", style="info" , disabled=TRUE)
             )),
           
             fluidRow(hr(), div(style = "margin-top: -23px"),
                      column(12, align="left",
                             column(7,checkboxInput("genes_sparse_select", 
                                                    label = "Sparse selection", value = FALSE) 
                             ))),
             conditionalPanel(condition='input.genes_sparse_select == true', 
                              fluidRow(column(4,align ="left", 
                                              shinyWidgets::pickerInput(
                                                inputId = "pick_spare_genes_rows",
                                                label = "Rows:",
                                                choices = NULL,
                                                multiple=TRUE,
                                                options = list("actions-box" = TRUE)
                                              )),
                                       column(4, align ="left",
                                              shinyWidgets::pickerInput(
                                                inputId = "pick_spare_genes_cols",
                                                label = "Columns:",
                                                choices = NULL,
                                                multiple=TRUE,
                                                options = list("actions-box" = TRUE)
                                              )),
                                       column(2, align="left", 
                                              div(style = "margin-top: 25px"),
                                              shinyBS::bsButton("visualize_gene_sparsemodal", "Visualize Sparsed Selection", style="info" , disabled=FALSE)
                                       ) 
                              )),
      
             fluidRow(column(12, align="left",
                             InteractiveComplexHeatmapOutput(heatmap_id="gene_mod_ui",
                                                             height1 = 2725,
                                                             width1=800,
                                                             height2 = 500,
                                                             width2=800,
                                                             width3=1000,
                                                             layout= "1|2|3", compact=FALSE,
                                                             action="click",
                                                             output_ui = htmlOutput("info_genes"), 
                                                             output_ui_float=FALSE)
             )
             ) 
    ),
    
    tabPanel("Properties Visualization",  
             fluidRow(column(3, align= "left",
                             br(),
                             shinyBS::bsButton("visualize_properties", "Visualize" , disabled=TRUE)
             )
             ), 
             
             fluidRow(hr(), div(style = "margin-top: -23px"),
                      column(12, align="left",
                             column(7,checkboxInput("properties_sparse_select", 
                                                    label = "Sparse selection", value = FALSE) 
                             ))),
             conditionalPanel(condition='input.properties_sparse_select == true', 
                              fluidRow(column(4,align ="left", 
                                              shinyWidgets::pickerInput(
                                                inputId = "pick_spare_properties_rows",
                                                label = "Rows:",
                                                choices = NULL,
                                                multiple=TRUE,
                                                options = list("actions-box" = TRUE)
                                              )),
                                       column(4, align ="left",
                                              shinyWidgets::pickerInput(
                                                inputId = "pick_spare_properties_cols",
                                                label = "Columns:",
                                                choices = NULL,
                                                multiple=TRUE,
                                                options = list("actions-box" = TRUE)
                                              )),
                                       column(2, align="left",
                                              div(style = "margin-top: 25px"),
                                              shinyBS::bsButton("visualize_properties_sparsemodal", "Visualize Sparsed Selection", style="info" , disabled=FALSE)
                                      )
                              )), 
             
             fluidRow(column(12, align="left",
                             InteractiveComplexHeatmapOutput(heatmap_id="propertiesHTM_ui", 
                                                             height1 = 550, 
                                                             width1=800,
                                                             height2 = 550, 
                                                             width2=800,
                                                             width3=1000,
                                                             layout= "1|2|3", compact=FALSE,
                                                             action="click",
                                                             output_ui = htmlOutput("info_properties"), 
                                                             output_ui_float=FALSE)
             )
             ) 
             
    ) 
    
  ) # main panl 
) #fludipage end ui
