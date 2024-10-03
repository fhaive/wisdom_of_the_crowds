library(dplyr) 
library(magrittr)
server <- function(input, output, session) {
  
  #global variables
  gVars <- shiny::reactiveValues( 
    ht_gtable=NULL,
    expert_list=NULL,                  # stores all expert questionnaires in a list, recoded (0=no, 1=NA, 2=yes)
    expert_uniq_df=NULL,               # stores as unique dataframe the list of experts´ questionnaires 
    expert_ID=NULL,                    # stores the IDs of the experts
    original_scoreTable=NULL,          # stores original hazard data table
    subset_scoreTable=NULL,            # stores the data table cleaned by excluded cols  
    click_pos=NULL,                    # stores click on heatmap coordinates
    table_mat_class_subset=NULL,       # stores the table subset according the clicked selected material class
    clicked_mat_class=NULL,            # after clicking a cell in the main heatmap, if clicked on "material cateogory" here stores the material class
    brush_area=NULL,                   # stores row and column indices highlighted by brush area
    original_geneTable_lst=NULL,       # stores original gene model data table
    gene_model_df=NULL,                # store currently selected gene dataframe
    click_pos_gene=NULL,               # stores click on gene model heatmap coordinates
    biomart_ensembl=NULL,              # stores the human ensembl from biomart library
    annotation=NULL,                   # stores all the info from the gene symbols selected
    non_zero_gene_coord=NULL,          # stores rownames(gene) and colnames(endpoint) of the non zero subset identified by sub-heatmap  
    original_propertiesTable_lst=NULL, # stores properties vs endpoint tables
    properties_df=NULL,                # store currently selected properties dataframe
    click_pos_property=NULL,           # stores click on properties heatmap coordinates
    non_zero_prop_coord=NULL,          # stores rownames(property) and colnames(endpoint) of the non zero subset identified by sub-heatmap  
    ceil_max_val_hazard=NULL,          # stores the up-rounded max value in hazard file
    ceil_max_val_genes=NULL,           # stores the up-rounded max value in gene file
    ceil_max_val_properties=NULL,      # stores the up-rounded max value in properties file
    flagHZ=NULL,                       # stores the number of hazard sparse modal windows
    flagHZM=NULL,                      # stores the number of hazard sparse material class modal windows
    flagGE=NULL,                       # stores the number of gene sparse modal windows
    flagPR=NULL                        # stores the number of properties sparse modal windows 
  )
  
  #sources
  source("./functions/woc_explorer_functions.R")
  
  
  #palette
  hts_palette <- c("white", "grey", "black")
  hts_norm_palette<-c("blue", "white", "red")
  
  #global option alignment dt tables
  options(DT.options = list(columnDefs = list(list(className = "dt-center", targets="_all"))))
  
  #cleaning r memory at start
  invisible(gc())
  
  
  ## loading expert answers ##
  expertFilePath <- "./data/Combined_cleaned_responses_anon.xlsx"
  
  df_list <- read_excel_allsheets(filename = expertFilePath,tibble = FALSE)  
  df_list[] <- lapply(names(df_list), function(expertID) {
    tmp_df <- data.frame(df_list[[eval(expertID)]],check.names=FALSE)
    orig <-tmp_df$ENM
    tmp_df[tmp_df != c("y","Y","n","N")| is.na(tmp_df)] <- 1
    tmp_df[tmp_df == "y"] <- 2
    tmp_df[tmp_df == "n"] <- 0
    tmp_df$expertID <- expertID
    tmp_df$ENM <- orig
    tmp_df
  })
  
  updateSelectizeInput(session,
                       "expertID", choices =  names(df_list),
                       selected = character(0), 
                       options = list( 
                         placeholder = "Select the expert ID",
                         maxItems=1)
  )            
  
  gVars$expert_list <- df_list
  gVars$expertID <- names(df_list)
  gVars$expert_uniq_df <- data.table::rbindlist(df_list) 
  
  
  # plotting expert answers
  observeEvent(input$expertID,{
    req(!is.null(input$expertID))
    req(input$expertID != "")
    if (input$expertID != "Overall"){
      single_questionnaire <- gVars$expert_list[[eval(input$expertID)]]
      rownames(single_questionnaire)=single_questionnaire[,1]
      
      m <- as.matrix(single_questionnaire[,-c(1,20)]) 
      
      heatmap_expert_questions <- Heatmap(m, name="expert_tbl",show_row_names = T, show_column_names = T,   
                                          height = nrow(m)*unit(7, "mm"),   #5qui
                                          row_names_side = "right",
                                          column_names_side = "bottom", 
                                          col = c("0" = "blue", "1" = "white", "2" = "red"),
                                          heatmap_legend_param = list(direction = "horizontal", nrow=1, 
                                                                      at = c(0,1,2), labels = c("no toxic", "no answer", "toxic"), 
                                                                      title="Answer", title_position = "topcenter",
                                                                      border = "lightgrey"),
                                          rect_gp = gpar(col = "lightgrey", lwd = 0.5)
      )
      
      output$heatmap_expert_questions <-  renderPlot({
        draw(heatmap_expert_questions, 
             heatmap_legend_side = "top" )
      })
    } 
  })
  
  
  
  ## loading hazard scores ##
  hazardFilePath <- "./data/concern_scores.csv"
  req(!is.null(hazardFilePath))
  df <- read.csv(file=  hazardFilePath, header = TRUE) %>% data.frame(.)
  rownames(df)=df[,1]
  gVars$original_scoreTable <- df[,-1]
  gVars$flagHZ <- gVars$flagHZM <- 0 
  shinyjs::disable("hazard_sparse_select") 
  
  output$hazard_main_table <- DT::renderDataTable({
    req(input$showtables == "Hazard")
    req(!is.null(gVars$original_scoreTable))
    tab = as.data.frame(lapply(gVars$original_scoreTable, 
                               function(x) if(is.numeric(x)) round(x, 2) else x))
    datatable(tab, options = list(paging = FALSE, 
                                  columnDefs = list(
                                    list(className = 'dt-center', targets = '_all'))
    ), rownames = FALSE
    )
  })
  
  
  observeEvent(input$hazard_sparse_select, {
    req(!is.null(input$hazard_sparse_select))
    value<- input$hazard_sparse_select  
    updateCheckboxInput(session, "hazard_sparse_select", value = value)
    if(value) {
      shinyBS::updateButton(session, "mat_categ", style="default",disabled=TRUE )  
    }  
  })
  
  
  ## functions related to make hazard heatmaps interactive ## 
  #' Handle double-click events on template hazard heatmap
  #' @param gVars A list containing global variables.
  #' @param df A DFrame object containing `row_label` and `column_label` from the double-click event.
  #' @return An HTML output displaying the row index, row label, column index, column label, and the value at the clicked position in the hazard dataframe.
  #' 
  #' @details This function processes double-click events on a hazard heatmap. It identifies labels, coordinates and values of the clicked cell from the hazard dataframe.
  dblclick_template_hazard = function(gVars, df) {
    mdf<-gVars$subset_scoreTable
    
    if(!is.null(df)) {
      gVars$click_pos <- reactiveValues(coord= c(df$row_label,df$column_label))
      row_index <- which(rownames(mdf) == df$row_label, arr.ind = TRUE)
      column_index <- which(colnames(mdf) == df$column_label, arr.ind = TRUE)
      
      gVars$clicked_mat_class <- gVars$original_scoreTable[gVars$original_scoreTable$X0 %in% df$row_label,"X1"]
      out <- HTML(paste0("<pre> <p> Row index: ", row_index, "<br/> Row label: ", df$row_label,
                         "<br/> Column index: ", column_index, "<br/> Column label: ", df$column_label,
                         "<br/> value: ", round(mdf[row_index,column_index],3)))
    } else {out<-""}
    if(!is.null(df$row_index)){
      shinyBS::updateButton(session, "mat_categ", style="info",disabled=FALSE ) }
    
    return(out)
  }
  
  #' Update outputUI consequent to double-click events for the whole hazard heatmap
  #' @param df A DFrame object containing `row_label` and `column_label` from the double-click event.
  #' @param input The Shiny input object.
  #' @param output The Shiny output object.
  #' @param session The Shiny session object.
  #'
  #' @return None. This function updates the Shiny UI output `info` with the result of `dblclick_template_hazard`.
  #' @details This function processes double-click events on the entire hazard heatmap and updates the `info` UI output.
  dblclick_action_hazard_whole = function(df, input,output,session) {
    output[["info"]] = renderUI({
      req(isFALSE(input$hazard_sparse_select))
      output <- dblclick_template_hazard(gVars,df) 
      output 
    })
  }
  
  #' Update outputUI consequent to double-click events for the sparse hazard heatmap.
  #' @param df A DFrame object containing `row_label` and `column_label` from the double-click event.
  #' @param input The Shiny input object.
  #' @param output The Shiny output object.
  #' @param session The Shiny session object.
  #'
  #' @return None. This function updates the Shiny UI output (with an ID `info_hazard_ht_sparse` concatenated with a numeric index) with the result of `dblclick_template_hazard`.
  #' @details This function processes double-click events on the sparse hazard heatmap and updates the corresponding UI output.
  dblclick_action_hazard_sparse = function(df, input,output,session) {
    output[[paste0("info_hazard_ht_sparse",gVars$flagHZ)]] = renderUI({
      req(isTRUE(input$hazard_sparse_select))
      output <- dblclick_template_hazard(gVars,df) 
      output 
    })
  }  
  
  
  #' Handle brush events on template hazard heatmap 
  #' @param gVars A list containing global variables.
  #' @param df A DFrame object containing `row_label` and `column_label` from the brush event.
  #' 
  #' @return A list containing updated global variables and a tag list of UI outputs (two scatterplots and two HTML spacer).
  #' @details This function processes brush events on hazard heatmap. It identifies the brushed area, retrieves the corresponding values from the hazard dataframe, and generates a scatterplot showing the mean values.
  brush_template_hazard = function(gVars, df) {
    row_index_ids = unique(unlist(df$row_label))
    column_index_ids = unique(unlist(df$column_label))
    mdf<-gVars$subset_scoreTable 
    
    if(!is.null(df)) { 
      gVars$brush_area <- reactiveValues(coord= list(row_index_ids, column_index_ids))
      upperlimitHZ <- 0.1+ gVars$ceil_max_val_hazard  
      sdf <-  gVars$subset_scoreTable[gVars$brush_area$coord[[1]],gVars$brush_area$coord[[2]], drop=FALSE]        
      col_names <- colnames(sdf)
      row_names <- rownames(sdf)
      brush_area_classes <- lapply(row_names, function(single_mat) gVars$original_scoreTable[gVars$original_scoreTable$X0 %in% single_mat,"X1"])
      sdf_long_mat <- cbind(sdf, material= row_names) %>%
        reshape2::melt(., id.vars="material") %>% setNames(., c("material", "endpoint", "score"))
      
      sdf_long_class <- cbind(sdf, brush_area_classes= as.character(brush_area_classes)) %>%
        reshape2::melt(., id.vars="brush_area_classes") %>% setNames(., c("material_category", "endpoint", "score"))
      
      mat_scatterplot <- sdf_long_mat %>%
        group_by(endpoint, material) %>%
        summarise(mean_val = mean(score, na.rm = T,.groups = "drop_last")) %>% 
        ggplot( aes(x = mean_val, y = endpoint, color = material)) +
        geom_point(position=position_jitter(h=0.08))+ 
        labs(title= "Specific Material", x= "Mean score", y="Endpoint",color="Materials") + 
        scale_x_continuous(breaks = seq(0,1,0.25),limits = c(-0.1,upperlimitHZ) )+
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
              plot.title = element_text(hjust = 0.5))
      
      class_scatterplot <- sdf_long_class %>%
        group_by(endpoint, material_category) %>%
        summarise(mean_val = mean(score, na.rm = T)) %>%
        ggplot( aes(x = mean_val, y = endpoint, color = material_category)) +
        geom_point(position=position_jitter(h=0.08))+ 
        labs(title= "Material Category",x= "Mean score", y="Endpoint",color="Material Categories") + 
        scale_x_continuous(breaks = seq(0,1,0.25),limits = c(-0.1,upperlimitHZ) )+
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
              plot.title = element_text(hjust = 0.5)   
        )
      
      output = tagList()
      output[[1]] <- plotly::renderPlotly({mat_scatterplot})
      output[[2]] <- HTML("<br>")
      output[[3]] <- plotly::renderPlotly({class_scatterplot})
      output[[4]] <- HTML("<br>")
      output
    } else {output <-NULL}
    
    return(list(list(gVars), list(output)))
  }
  
  #' Update outputUI consequent to brush events for the whole hazard heatmap
  #' @param df A DFrame object containing `row_label` and `column_label` from the brush event.
  #' @param input The Shiny input object.
  #' @param output The Shiny output object.
  #' @param session The Shiny session object.
  #'
  #' @return None. This function updates the Shiny UI output `info` with the result of `brush_template_hazard`.
  #' @details This function processes brush events on the entire hazard heatmap and updates the `info` UI output.
  brush_action_hazard_whole = function(df, input,output,session) {
    updateTextInput(session, "show_mat_class",value="No material category selected")
    shinyBS::updateButton(session, "mat_categ" , style= "default" ,disabled=TRUE)
    
    output[["info"]] = renderUI({
      req(isFALSE(input$hazard_sparse_select))
      out_brush <- brush_template_hazard(gVars,df) 
      gVars <- out_brush[[1]][[1]]
      output <- out_brush[[2]]
      output 
    })
  }
  
  #' Update outputUI consequent to brush events for the sparse hazard heatmap
  #' @param df A DFrame object containing `row_label` and `column_label` from the brush event.
  #' @param input The Shiny input object.
  #' @param output The Shiny output object.
  #' @param session The Shiny session object.
  #'
  #' @return None. This function updates the Shiny UI output (with an ID `info_hazard_ht_sparse` concatenated with a numeric index) with the result of `brush_template_hazard`.
  #' @details This function processes brush events on the sparse hazard heatmap and updates the corresponding UI output.
  brush_action_hazard_sparse = function(df, input,output,session) {
    output[[paste0("info_hazard_ht_sparse",gVars$flagHZ)]] = renderUI({
      req(isTRUE(input$hazard_sparse_select))
      out_brush <- brush_template_hazard(gVars,df) 
      gVars <- out_brush[[1]][[1]]
      output <- out_brush[[2]]
      output 
    })
  }
  
  
  observeEvent(input$visualize,{
    req(!is.null(gVars$original_scoreTable))
    shinyjs::enable("hazard_sparse_select")
    
    dataf <- gVars$original_scoreTable 
    rownames(dataf)=dataf[,19]
    gVars$subset_scoreTable <- dataf[,-c(19,20)] 
    gVars$ceil_max_val_hazard <-ceiling(max(gVars$subset_scoreTable))
    
    updatePickerInput(session,
                      inputId = "pick_spare_hazard_rows",
                      choices = rownames(gVars$subset_scoreTable),   
                      options = list(liveSearch = TRUE,
                                     liveSearchPlaceholder = "Search entry")    
    )
    updatePickerInput(session,
                      inputId = "pick_spare_hazard_cols",
                      choices = colnames(gVars$subset_scoreTable),   
                      options = list(liveSearch = TRUE,
                                     liveSearchPlaceholder = "Search entry")    
    )
    
    m <- as.matrix(dataf[,-c(19:20)]) %>% rowsum_ordering(.)
    heatmap_table <- Heatmap(m, show_row_names = TRUE, show_column_names = TRUE,
                             name = "all_materials_tbl",
                             heatmap_height = unit(85, "cm"), 
                             show_row_dend = FALSE,
                             show_column_dend = FALSE, 
                             cluster_rows = FALSE,  
                             heatmap_legend_param = list(title="Score", title_position = "topcenter" )
    )
    
    heatmap_table = draw(heatmap_table)   
    makeInteractiveComplexHeatmap(input, output, session, heatmap_table,   "ht_ui",
                                  click_action = dblclick_action_hazard_whole, 
                                  brush_action = brush_action_hazard_whole 
    )
  })
  
  
  #modal for selection sparsed rows/cols
  observeEvent(input$visualize_sparsemodal,{
    req(isTRUE(input$hazard_sparse_select))
    
    #sets sizes of modal window elements
    height1 <- ifelse(is.null(input$pick_spare_hazard_rows) | 
                                  length(input$pick_spare_hazard_rows) == nrow(gVars$subset_scoreTable) , 
                                2725, length(input$pick_spare_hazard_rows)*20+200)
    width1<- ifelse(is.null(input$pick_spare_hazard_cols) | 
                                         length(input$pick_spare_hazard_cols) == ncol(gVars$subset_scoreTable) , 
                                       850, length(input$pick_spare_hazard_cols)*25+200 )

    width2 <- width1
    height2<-  ifelse(is.null(input$pick_spare_hazard_rows) | 
                                 length(input$pick_spare_hazard_rows) == nrow(gVars$subset_scoreTable) , 
                               (height1)*0.4, (height1)/2 +200)
    
    dims <- c(height1, width1, height2, width2)
     
    if ( gVars$flagHZ > 0 ) {
      old_modal_cleaner("info_hazard_ht_sparse",gVars$flagHZ)
    }
    
    gVars$flagHZ <-gVars$flagHZ+1
    showModal(open_sparse_set_modal(dims, paste0("hazard_ht_sparse",gVars$flagHZ)))
    runjs("setModalWidth('shiny-modal', '1200');")
    
    if (is.null(input$pick_spare_hazard_cols) ){
      cols <- colnames(gVars$subset_scoreTable)
    } else {cols <- input$pick_spare_hazard_cols}
    if (is.null(input$pick_spare_hazard_rows) ){
      rows <- rownames(gVars$subset_scoreTable)
    } else {rows <-input$pick_spare_hazard_rows}
    input$pick_spare_hazard_cols
    input$pick_spare_hazard_rows 
    
    sparsed_dataf <- gVars$subset_scoreTable[rows,cols] %>% data.frame()
    row.names(sparsed_dataf) <-rows
    colnames(sparsed_dataf)<-cols
    m_sparse <- as.matrix(sparsed_dataf, rownames=TRUE) %>% rowsum_ordering(.)
    colnames(m_sparse)<- cols    # needed to fix cases with one row or one column
    
    heatmap_sparse_subset_table <- Heatmap(m_sparse, show_row_names = TRUE, show_column_names = TRUE,
                                           name = "sparse_selection_tbl", 
                                           show_row_dend = FALSE,
                                           show_column_dend = FALSE,  
                                           cluster_rows = FALSE,  
                                           col = colorRamp2(c(0, 0.5,1), hts_norm_palette),
                                           heatmap_legend_param = list(title="Score", title_position = "topcenter" )
    )
    
    heatmap_sparse_subset_table = draw(heatmap_sparse_subset_table) 
    makeInteractiveComplexHeatmap(input, output, session, heatmap_sparse_subset_table, 
                                  heatmap_id= paste0("hazard_ht_sparse",gVars$flagHZ),
                                  click_action = dblclick_action_hazard_sparse, 
                                  brush_action = brush_action_hazard_sparse  
    )
  })
  
  
  
  ## shows the material class if clicked a cell##
  observe({
    gVars$click_pos
    if (is.null(gVars$click_pos)){    #no click on heatmap yet
      updateTextInput(session, "show_mat_class",value="No material category selected")
    } else { 
      updateTextInput(session, "show_mat_class", value= gVars$clicked_mat_class)
    }
  })
  
  
  # only material class modal window subset map for class of materials
  observeEvent(input$mat_categ,{ 
    req(!is.null(gVars$click_pos$coord))         
    clicked_row_lab  <- gVars$click_pos$coord[1]
    gVars$clicked_mat_class <- gVars$original_scoreTable[gVars$original_scoreTable$X0 %in% clicked_row_lab,"X1"]
    gVars$table_mat_class_subset <- gVars$original_scoreTable[gVars$original_scoreTable$X1 %in% gVars$clicked_mat_class,] 
    rownames(gVars$table_mat_class_subset)=gVars$table_mat_class_subset[,19]
    gVars$table_mat_class_subset <- gVars$table_mat_class_subset[,-c(19,20)]
    updateTextInput(session, "selected_class_material", value= gVars$clicked_mat_class)
    
    height1 <- nrow(gVars$table_mat_class_subset)*20+200
    width1<- ncol(gVars$table_mat_class_subset)*35 +200
    width2<- (width1) 
    height2 <- (height1)*0.3 +200
    dims <-c(height1, width1, height2, width2)
    
    if ( gVars$flagHZM > 0 ) {
      old_modal_cleaner("info_mat_categ", gVars$flagHZM)
    }
    
    gVars$flagHZM <- gVars$flagHZM+1
    showModal(subset_mat_class_modal(dims, gVars$flagHZM))
    runjs("setModalWidth('shiny-modal', '1200');")
    
    m_mat_subset <- as.matrix(gVars$table_mat_class_subset) %>% rowsum_ordering(.)
    heatmap_mat_class_subset_table <- Heatmap(m_mat_subset, show_row_names = TRUE, show_column_names = TRUE,
                                              name="material_class_tbl",
                                              show_row_dend = FALSE,
                                              show_column_dend = FALSE, 
                                              cluster_rows = FALSE,  
                                              heatmap_legend_param = list(title="Score", title_position = "topcenter" )
    )
    heatmap_mat_class_subset_table = draw(heatmap_mat_class_subset_table) 
    makeInteractiveComplexHeatmap(input, output, session, heatmap_mat_class_subset_table,  
                                  #"ht_mat_categ", 
                                  paste0("ht_mat_categ",gVars$flagHZM),
                                  click_action = dblclick_action_mat_categ, 
                                  brush_action = brush_action_mat_categ 
    )
  })
  
  
  # Return a modal dialog window to show the material class 
  subset_mat_class_modal <- function(dims, ht_id) {
    modalDialog( 
      div(style= "overflow-x: scroll; overflow-y: scroll; } ",
      fluidRow(br()),
      fluidRow(column(3, disabled(textInput("selected_class_material", label="Material Class", value="No material category selected") ))),
      fluidRow(column(12, 
                      InteractiveComplexHeatmapOutput(heatmap_id=paste0("ht_mat_categ",ht_id),
                                                      height1 = dims[1], 
                                                      width1=dims[2],
                                                      height2 = dims[3], 
                                                      width2=dims[4],
                                                      width3=900, 
                                                      layout= "1|2|3", compact=FALSE,
                                                      action="dblclick",
                                                      response=c("dblclick", "brush"),
                                                      output_ui = htmlOutput(paste0("info_mat_categ",ht_id)), 
                                                      output_ui_float=FALSE 
                      )
      ))
      ),
      
      size="xl",
      easyClose = F, 
      footer = tagList(
        div(style="width: 100%; display: flex; justify-content: space-between;", 
            modalButton("Dismiss")
        )
      ),
      id="customModalmat"
    )
  }
  
  
  #' Handle brush events on categorical material matrix
  #' @param df A DFrame containing `row_label` and `column_index` from the brush event.
  #' @param input The Shiny input object.
  #' @param output The Shiny output object.
  #' @param session The Shiny session object.
  #'
  #' @return None. This function updates the Shiny UI output (with an ID `info_mat_categ` concatenated with a numeric index) with a scatterplot and a HTML spacer.
  #' @details This function processes brush events on a categorical material matrix heatmap. It identifies the brushed area, retrieves the corresponding values from the subset hazard score table, and generates a scatterplot showing the mean values.
  brush_action_mat_categ = function(df, input,output,session) { 
    row_names_ids <- unlist(df$row_label)
    column_indices = unique(unlist(df$column_index))
    mdf<-gVars$subset_scoreTable 
    
    output[[paste0("info_mat_categ",gVars$flagHZM)]] = renderUI({
      if(!is.null(df)) { 
        
        gVars$brush_area <- reactiveValues(coord= list(row_names_ids,column_indices))
        upperlimitHZ <- 0.1+ gVars$ceil_max_val_hazard  
        cat_subset <- mdf[row_names_ids, column_indices, drop = FALSE]
        row_names <- rownames(cat_subset)
        whole_cat_subset_long <- cbind(cat_subset, material= row_names) %>%
          reshape2::melt(., id.vars="material") %>% setNames(., c("material", "endpoint", "score"))
        title_legend <- paste("Type of", input$selected_class_material)
        
        mat_scatterplot <- whole_cat_subset_long %>%
          group_by(endpoint, material) %>%
          summarise(mean_val = mean(score, na.rm = T)) %>%
          ggplot(aes(x = mean_val , y = endpoint, color = material)) +
          geom_point(position=position_jitter(h=0.08))+ 
          labs(title= "Specific Material", color="Material", x= "Mean score", y="Endpoint") + 
          scale_x_continuous(breaks = seq(0,1,0.25),limits = c(-0.1, upperlimitHZ ) )+
          theme(axis.text.x = element_text(angle = 60, hjust = 1),
                plot.title = element_text(hjust = 0.5))
        
        output = tagList()
        output[[1]] <- plotly::renderPlotly({mat_scatterplot})
        output[[2]] <- HTML("<br>")
        output 
      } 
    })
  }
  
  
  #' Handle double-click events on categorical material matrix
  #' @param df A DFrame containing `row_label` and `column_label` from the double-click event.
  #' @param input The Shiny input object.
  #' @param output The Shiny output object.
  #' @param session The Shiny session object.
  #'
  #' @return None. This function updates the Shiny UI output (with an ID `info_mat_categ` concatenated with a numeric index) with the row index, row label, column index, column label, and the value at the clicked position in the subset hazard score table.
  #' @details This function processes double-click events on a categorical material matrix heatmap. It identifies the row and column indices and labels of the clicked cell and retrieves the corresponding value from the subset hazard score table.
  dblclick_action_mat_categ = function(df, input,output,session) { 
    mdf<-gVars$subset_scoreTable
    output[[paste0("info_mat_categ",gVars$flagHZM)]] = renderUI({
      if(!is.null(df)) {
        gVars$click_pos <- reactiveValues(coord= c(df$row_label,df$column_label))
        row_index <- which(rownames(mdf) == df$row_label, arr.ind = TRUE)
        col_index <- which(colnames(mdf) == df$column_label, arr.ind = TRUE)
        HTML(paste0("<pre> <p> Row index: ", row_index, "<br/> Row label: ", df$row_label,
                    "<br/> Column index: ", col_index, "<br/> Column label: ", df$column_label,
                    "<br/> value: ", round(mdf[row_index,col_index],3)))
      } 
    })
  }
  
  
  
  ## loading current biomart for genes present in following file ##
  biomartFilePath <- "./data/biomart_17_05_2024.xlsx"
  req(!is.null(biomartFilePath)) 
  ens_list <- read_excel_allsheets(filename = biomartFilePath,tibble = FALSE) 
  gVars$biomart_ensembl <- ens_list[[1]] %>% 
    mutate (Ensemble_Url = createLink(.$ensembl_gene_id,.$hgnc_symbol))
  
  
  ## loading weighted important genes data ##
  weightedGenePath <- "./data/important_genes_weighted.xlsx"
  req(!is.null(weightedGenePath)) 
  df_list <- read_excel_allsheets(filename = weightedGenePath,tibble = FALSE) 
  
  gVars$original_geneTable_lst <- df_list
  gVars$flagGE <-0
  shinyBS::updateButton(session, "visualize_model" , style= "info" ,disabled=FALSE)
  shinyjs::disable("genes_sparse_select")
  
  output$gene_main_table <- DT::renderDataTable({
    req(input$showtables == "Genes") 
    df <- gVars$original_geneTable_lst[[1]] 
    df = as.data.frame(lapply(df, 
                              function(x) if(is.numeric(x)) round(x, 2) else x))
    datatable(df, options = list(paging = FALSE,autoWidth = TRUE),
              rownames = FALSE)
  })
  
  
  observeEvent(input$visualize_model, {
    req(!is.null(gVars$original_geneTable_lst) & isFALSE(input$genes_sparse_select)) 
    dataf <- gVars$original_geneTable_lst[[1]] 
    
    rownames(dataf)=dataf[,1]
    dataf <- dataf[,-1]   
    gVars$gene_model_df <- dataf
    gVars$ceil_max_val_genes <- ceiling(max(dataf))
    
    shinyjs::enable("genes_sparse_select")
    updatePickerInput(session,
                      inputId = "pick_spare_genes_rows",
                      choices = rownames(gVars$gene_model_df),   
                      options = list(liveSearch = TRUE,
                                     liveSearchPlaceholder = "Search entry")    
    )
    updatePickerInput(session,
                      inputId = "pick_spare_genes_cols",
                      choices = colnames(gVars$gene_model_df),   
                      options = list(liveSearch = TRUE,
                                     liveSearchPlaceholder = "Search entry")    
    )
    
    m <- as.matrix(dataf)  %>% rowsum_ordering(.) %>% as.matrix(., rownames=TRUE)
    heatmap_table <- Heatmap(m, show_row_names = TRUE, show_column_names = TRUE,  
                             name="genes_tbl",
                             heatmap_height = unit(95, "cm"),
                             show_row_dend = FALSE,
                             show_column_dend = FALSE, 
                             cluster_rows = FALSE, 
                             col = colorRamp2(c(0, (gVars$ceil_max_val_genes)/2, gVars$ceil_max_val_genes), hts_palette),
                             heatmap_legend_param = list(title="Score", title_position = "topcenter" ))
    heatmap_table = draw(heatmap_table)  
    makeInteractiveComplexHeatmap(input, output, session, heatmap_table,   "gene_mod_ui" ,
                                  click_action = dblclick_action_gene_whole, 
                                  brush_action = brush_action_gene_whole 
    )
  })
  
  
  #modal for selection sparsed rows/cols for genes
  observeEvent(input$visualize_gene_sparsemodal,{
    req(isTRUE(input$genes_sparse_select))
    
    height1 <- ifelse(is.null(input$pick_spare_genes_rows) | 
                        length(input$pick_spare_genes_rows) == nrow(gVars$gene_model_df) , 
                      2725, length(input$pick_spare_genes_rows)*20+200)
    width1<- ifelse(is.null(input$pick_spare_genes_cols) | 
                      length(input$pick_spare_genes_cols) == ncol(gVars$gene_model_df) , 
                    600, length(input$pick_spare_genes_cols)*25+200 )
    width2 <- width1
    height2<- ifelse(is.null(input$pick_spare_genes_rows) | 
                       length(input$pick_spare_genes_rows) == nrow(gVars$gene_model_df) , 
                     height1*0.2, (height1)/2 +200)
    
    dims <- c(height1, width1, height2, width2)
    
    if ( gVars$flagGE > 0 ) {
      old_modal_cleaner("info_genes_ht_sparse", gVars$flagGE)
    }
    
    gVars$flagGE <- gVars$flagGE +1
    showModal(open_sparse_set_modal(dims, paste0("genes_ht_sparse",gVars$flagGE))) 
    runjs("setModalWidth('shiny-modal', '1200');")
    
    if (is.null(input$pick_spare_genes_cols) ){
      colsg <- colnames(gVars$gene_model_df)
    } else {colsg <- input$pick_spare_genes_cols}
    if (is.null(input$pick_spare_genes_rows) ){
      rowsg <- rownames(gVars$gene_model_df)
    } else {rowsg <-input$pick_spare_genes_rows}
    input$pick_spare_genes_cols
    input$pick_spare_genes_rows 
    
    sparsed_datag <- gVars$gene_model_df[rowsg,colsg] %>% data.frame()
    row.names(sparsed_datag) <- rowsg
    colnames(sparsed_datag)<- colsg
    m_sparseg <- as.matrix(sparsed_datag, rownames=TRUE)  %>% rowsum_ordering(.)# %>% as.matrix(., rownames=TRUE)
    colnames(m_sparseg)<- colsg    # needed to fix cases with one row or one column
    
    heatmap_sparse_subset_tableG <- Heatmap(m_sparseg, show_row_names = TRUE, show_column_names = TRUE,
                                            name = "sparse_selection_tbl",
                                            show_row_dend = FALSE,
                                            show_column_dend = FALSE,
                                            cluster_rows = FALSE, 
                                            col = colorRamp2(c(0, (gVars$ceil_max_val_genes)/2, gVars$ceil_max_val_genes), hts_palette), 
                                            heatmap_legend_param = list(title="Score", title_position = "topcenter" )
    )
    
    heatmap_sparse_genes_subset_tableG = draw(heatmap_sparse_subset_tableG) 
    makeInteractiveComplexHeatmap(input, output, session, heatmap_sparse_genes_subset_tableG,  
                                   paste0("genes_ht_sparse",gVars$flagGE),
                                  click_action = dblclick_action_gene_sparse, 
                                  brush_action = brush_action_gene_sparse  
    )
  })
  
  
  ## functions related to make gene heatmaps interactive ## 
  #' Handle double-click events on template gene heatmap
  #' @param gVars A list containing global variables.
  #' @param df A DFrame object containing `row_label` and `column_label` from the double-click event.
  #' @return A tag list of HTML and DataTable outputs displaying the row index, row label, column index, column label, the value at the clicked position in the gene model dataframe, and gene annotation details if available.
  #' 
  #' @details This function processes double-click events on a gene heatmap. It identifies labels, coordinates and values of the clicked cell from the gene dataframe. It also provides additional gene annotation details if the value is not-null ("Ensembl url", "Ensembl ID", "Gene description" and HGNC id).
  dblclick_template_genes = function(gVars , df) {
    gene_df<-gVars$gene_model_df
    
    if(!is.null(df)) {
      gVars$click_pos_gene <- reactiveValues(coord= list(df$row_label,df$column_label))
      row_index <- which(rownames(gene_df) == gVars$click_pos_gene$coord[[1]], arr.ind = TRUE)
      column_index <- which(colnames(gene_df) == gVars$click_pos_gene$coord[[2]], arr.ind = TRUE)
      
      info_block <- HTML(paste0("<pre> <p> Row index: ", row_index, "<br/> Row label: ", gVars$click_pos_gene$coord[[1]],
                                "<br/> Column index: ", column_index, "<br/> Column label: ", gVars$click_pos_gene$coord[[2]],
                                "<br/> value: ", round(gene_df[row_index,column_index], 3),"</pre>"))
      
      non_null_gene_flag <- gene_df[row_index,column_index] != 0
      if (non_null_gene_flag){
        gVars$annotation <- gVars$biomart_ensembl %>% .[.$hgnc_symbol %in% gVars$click_pos_gene$coord[[1]], ] 
        tmp<- gVars$annotation %>% dplyr::select(Ensemble_Url, ensembl_gene_id, description)
        names(tmp) <- c("Ensembl url", "Ensembl ID","Gene description") 
        rownames(tmp) <- NULL
      } else {tmp <-NULL  }
      
      output = tagList()
      output[[1]] <- info_block
      output[[2]] <- DT::renderDT( tmp, rownames=FALSE,escape=FALSE, server=FALSE,  selection = "single",
                                   options = list(searching=FALSE, paging = FALSE 
                                   ))
      output[[3]] <-  HTML("<br>")  
      output
    } else {output <-""}
    return (output)
  }
  
  
  #' Update outputUI consequent to double-click events for the whole gene heatmap
  #' @param df A DFrame object containing `row_label` and `column_label` from the double-click event.
  #' @param input The Shiny input object.
  #' @param output The Shiny output object.
  #' @param session The Shiny session object.
  #'
  #' @return None. This function updates the Shiny UI output `info_genes` with the result of `dblclick_template_genes`.
  #' @details This function processes double-click events on the entire gene heatmap and updates the `info_genes` UI output.
  dblclick_action_gene_whole = function(df, input,output,session) {
    output[["info_genes"]] = renderUI({
      req(isFALSE(input$genes_sparse_select))
      output <- dblclick_template_genes(gVars,df) 
      output 
    })
  }
  
  
  #' Update outputUI consequent to double-click events for the sparse gene heatmap.
  #' @param df A DFrame object containing `row_label` and `column_label` from the double-click event.
  #' @param input The Shiny input object.
  #' @param output The Shiny output object.
  #' @param session The Shiny session object.
  #'
  #' @return None. This function updates the Shiny UI output (with an ID `info_genes_ht_sparse` concatenated with a numeric index) with the result of `dblclick_template_genes`.
  #' @details This function processes double-click events on the sparse gene heatmap and updates the corresponding UI output.
  dblclick_action_gene_sparse = function(df, input,output,session) {
    output[[paste0("info_genes_ht_sparse",gVars$flagGE)]] = renderUI({
      req(isTRUE(input$genes_sparse_select))
      output <- dblclick_template_genes(gVars,df) 
      output 
    })
  }
  
  
  #' Handle brush events on template gene heatmap 
  #' @param gVars A list containing global variables.
  #' @param df A DFrame object containing `row_label` and `column_label` from the brush event.
  #' 
  #' @return A list containing updated global variables and a tag list of UI outputs (a tDataTable with not-null gene annotation info and HTML contents).
  #' @details This function processes brush events on a gene heatmap. It identifies the brushed area, retrieves the corresponding values from the gene dataframe, and generates a table showing "Ensembl url", "Ensembl ID","Gene description" and HGNC id of non-null genes. 
  brush_template_genes = function(gVars, df) { 
    row_index_ids = unique(unlist(df$row_label))
    column_index_ids = unique(unlist(df$column_label)) 
    
    mdf<-gVars$gene_model_df
    if(!is.null(df)) { 
      gVars$brush_area_genes <- reactiveValues(coord= list(row_index_ids, column_index_ids))
      sub_dfmodel <-  gVars$gene_model_df[gVars$brush_area_genes$coord[[1]],gVars$brush_area_genes$coord[[2]], drop=FALSE]
      
      non_zero_values_ids <- which(sub_dfmodel != 0, arr.ind = TRUE) %>% data.frame(.)
      non_zero_values_ids$rowNMS <- rownames(non_zero_values_ids)
      tmp_col <- colnames(sub_dfmodel)
      non_zero_values_ids$colNMS <-  tmp_col[non_zero_values_ids$col]
      
      gVars$annotation <- gVars$biomart_ensembl %>% .[.$hgnc_symbol %in% non_zero_values_ids$rowNMS, ] %>%
        mutate (Ensemble_Url = createLink(.$ensembl_gene_id,.$hgnc_symbol)) 
      
      gVars$non_zero_gene_coord <- non_zero_values_ids %>% setNames(., c("row", "col", "Gene", "Endpoint"))
      
      tmp<- gVars$annotation %>% dplyr::select(Ensemble_Url, ensembl_gene_id, description)
      names(tmp) <- c("Ensembl url", "Ensembl ID","Gene description") 
      rownames(tmp) <- NULL
      
      output = tagList()
      output[[1]] <- h4(HTML("<br><b>Non-null selected genes</b> <br>"))
      output[[2]] <- DT::renderDT( tmp, rownames=FALSE,escape=FALSE, server=FALSE,  selection="single",
                                   options = list(searching=FALSE, paging = FALSE 
                                   ) )
      
      output[[3]] <- HTML("<br>")  
      output
    } else {output <- NULL} 
    return(list(list(gVars), list(output)))
  }
  
  
  #' Update outputUI consequent to brush events for the whole gene heatmap
  #' @param df A DFrame object containing `row_label` and `column_label` from the brush event.
  #' @param input The Shiny input object.
  #' @param output The Shiny output object.
  #' @param session The Shiny session object.
  #'
  #' @return None. This function updates the Shiny UI output `info_genes` with the result of `brush_template_genes`.
  #' @details This function processes brush events on the entire gene heatmap and updates the `info_genes` UI output.
  brush_action_gene_whole = function(df, input,output,session) {
    output[["info_genes"]] = renderUI({
      req(isFALSE(input$genes_sparse_select))
      out_brush <- brush_template_genes(gVars,df) 
      gVars <- out_brush[[1]][[1]]
      output <- out_brush[[2]]
      output 
    })
  }
  
  #' Update outputUI consequent to brush events for the sparse gene heatmap
  #' @param df A DFrame object containing `row_label` and `column_label` from the brush event.
  #' @param input The Shiny input object.
  #' @param output The Shiny output object.
  #' @param session The Shiny session object.
  #'
  #' @return None. This function updates the Shiny UI output (with an ID `info_genes_ht_sparse` concatenated with a numeric index) with the result of `brush_template_genes`.
  #' @details This function processes brush events on the sparse gene heatmap and updates the corresponding UI output.
  brush_action_gene_sparse = function(df, input,output,session) {
    output[[paste0("info_genes_ht_sparse",gVars$flagGE)]] = renderUI({
      req(isTRUE(input$genes_sparse_select))
      out_brush <- brush_template_genes(gVars,df) 
      gVars <- out_brush[[1]][[1]]
      output <- out_brush[[2]]
      output 
    })
  }
  
  
  
  ## loading properties file and tables ## 
  weightedPropertiesPath <- "./data/important_physchem_weighted.xlsx"
  req(!is.null(weightedPropertiesPath)) 
  dfp_list <- read_excel_allsheets(filename = weightedPropertiesPath,tibble = FALSE) 
  gVars$flagPR <-0
  shinyjs::disable("properties_sparse_select")
  
  gVars$original_propertiesTable_lst <- dfp_list
  shinyBS::updateButton(session, "visualize_properties" , style= "info" ,disabled=FALSE) 
  
  output$properties_main_table <- DT::renderDataTable({
    req(input$showtables == "Properties")
    df <- gVars$original_propertiesTable_lst[[1]]
    df = as.data.frame(lapply(df, 
                              function(x) if(is.numeric(x)) round(x, 2) else x))
    datatable(df, options = list(paging = FALSE),
              rownames = FALSE)
  })
  
  
  #Proeprteis sparse selection
  observeEvent(input$properties_sparse_select, {
    req(!is.null(input$properties_sparse_select))
    valz<- input$properties_sparse_select  
    updateCheckboxInput(session, "properties_sparse_select", value = valz) 
  })
  
  
  observeEvent(input$visualize_properties, {
    req(!is.null(gVars$original_propertiesTable_lst))
    dataf <- gVars$original_propertiesTable_lst[[1]] 
    
    rownames(dataf)=dataf[,1]
    dataf <- dataf[,-1]   
    gVars$properties_df <- dataf 
    gVars$ceil_max_val_properties <-ceiling(max(gVars$properties_df))
    
    shinyjs::enable("properties_sparse_select")
    updatePickerInput(session,
                      inputId = "pick_spare_properties_rows",
                      choices = rownames(gVars$properties_df),   
                      options = list(liveSearch = TRUE,
                                     liveSearchPlaceholder = "Search entry")    
    )
    updatePickerInput(session,
                      inputId = "pick_spare_properties_cols",
                      choices = colnames(gVars$properties_df),   
                      options = list(liveSearch = TRUE,
                                     liveSearchPlaceholder = "Search entry")    
    )
    
    m <- as.matrix(dataf)  %>% rowsum_ordering(.) %>% as.matrix(., rownames=TRUE)
    heatmap_table_properties <- Heatmap(m, show_row_names = TRUE, show_column_names = TRUE, 
                                        name="properties_tbl",
                                        heatmap_height = unit(18, "cm"),
                                        show_row_dend = FALSE,
                                        show_column_dend = FALSE, 
                                        cluster_rows = FALSE, 
                                        col = colorRamp2(c(0, (gVars$ceil_max_val_properties)/2, gVars$ceil_max_val_properties), hts_palette),
                                        heatmap_legend_param = list(title="Score", title_position = "topcenter" ))
    heatmap_table_properties = draw(heatmap_table_properties)  
    makeInteractiveComplexHeatmap(input, output, session, heatmap_table_properties,   "propertiesHTM_ui"                                  ,
                                  click_action = dblclick_action_properties_whole, 
                                  brush_action = brush_action_properties_whole
    )
  })
  
  
  #modal for selection sparsed rows/cols
  observeEvent(input$visualize_properties_sparsemodal,{
    req(isTRUE(input$properties_sparse_select))
    
    height1 <- ifelse(is.null(input$pick_spare_properties_rows) | 
                                  length(input$pick_spare_properties_rows) == nrow(gVars$properties_df) , 
                                550, length(input$pick_spare_properties_rows)*20+200)
    width1<- ifelse(is.null(input$pick_spare_properties_cols) | 
                                length(input$pick_spare_properties_cols) == ncol(gVars$properties_df) , 
                              800, length(input$pick_spare_properties_cols)*25+200 )
    width2 <- width1
    height2<- ifelse(is.null(input$pick_spare_properties_rows) | 
                                 length(input$pick_spare_properties_rows) == nrow(gVars$properties_df) , 
                               height1*0.7, (height1)/2 +200)
    
    dims <- c(height1, width1, height2, width2)
    
    if ( gVars$flagPR > 0 ) {
      old_modal_cleaner("info_properties_ht_sparse",gVars$flagPR)
    }
    
    gVars$flagPR <-gVars$flagPR+1
    showModal(open_sparse_set_modal(dims, paste0("properties_ht_sparse",gVars$flagPR)))
    runjs("setModalWidth('shiny-modal', '1200');")
    
    if (is.null(input$pick_spare_properties_cols) ){
      cols <- colnames(gVars$properties_df)
    } else {cols <- input$pick_spare_properties_cols}
    if (is.null(input$pick_spare_properties_rows) ){
      rows <- rownames(gVars$properties_df)
    } else {rows <-input$pick_spare_properties_rows}
    input$pick_spare_properties_cols
    input$pick_spare_properties_rows 
    
    sparsed_dataf <- gVars$properties_df[rows,cols] %>% data.frame()
    
    row.names(sparsed_dataf) <- rows
    colnames(sparsed_dataf)<- cols
    m_sparse <- as.matrix(sparsed_dataf, rownames=TRUE) %>% rowsum_ordering(.)# %>% as.matrix(., rownames=TRUE)
    colnames(m_sparse)<- cols    # needed to fix cases with one row or one column
    
    heatmap_sparse_subset_table <- Heatmap(m_sparse, show_row_names = TRUE, show_column_names = TRUE,
                                           name = "sparse_selection_tbl2",
                                           show_row_dend = FALSE,
                                           show_column_dend = FALSE, 
                                           cluster_rows = FALSE, 
                                           col = colorRamp2(c(0, (gVars$ceil_max_val_properties)/2, gVars$ceil_max_val_properties), hts_palette),
                                           heatmap_legend_param = list(title="Score", title_position = "topcenter" )
    )
    
    heatmap_sparse_properties_subset_table = draw(heatmap_sparse_subset_table)  
    makeInteractiveComplexHeatmap(input, output, session, heatmap_sparse_properties_subset_table, 
                                  heatmap_id= paste0("properties_ht_sparse",gVars$flagPR),
                                  click_action = dblclick_action_properties_sparse, 
                                  brush_action = brush_action_properties_sparse  
    )   
  })
  
  
  ## functions related to make properties heatmaps interactive ## 
  #' Handle double-click events on template properties heatmap
  #' @param gVars A list containing global variables.
  #' @param df A DFrame object containing `row_label` and `column_label` from the double-click event.
  #'
  #' @return An HTML output displaying the row index, row label, column index, column label, and the value at the clicked position in the properties dataframe.
  #' @details This function processes double-click events on a properties heatmap. It identifies labels, coordinates and values of the clicked cell from the properties dataframe.
  dblclick_template_properties <-function(gVars, df){
    row_index_ids = unique(unlist(df$row_label))
    column_index_ids = unique(unlist(df$column_label))
    properties_df<-gVars$properties_df
    
    if(!is.null(df)) {
      gVars$click_pos_property <- reactiveValues(coord= list(df$row_label,df$column_label))
      row_index <- which(rownames(properties_df) == gVars$click_pos_property$coord[[1]], arr.ind = TRUE)
      column_index<- which(colnames(properties_df) == gVars$click_pos_property$coord[[2]], arr.ind = TRUE)
      out<- HTML(paste0("<pre> <p> Row index: ", row_index, "<br/> Row label: ", gVars$click_pos_property$coord[[1]],
                        "<br/> Column index: ", column_index, "<br/> Column label: ", gVars$click_pos_property$coord[[2]],
                        "<br/> value: ", round(
                          properties_df[row_index,column_index],3)
      ))
    } else {out<-""}
    return(out)
  }
  
  #' Update outputUI consequent to double-click events for the whole properties heatmap.
  #' @param df A DFrame object containing `row_label` and `column_label` from the double-click event.
  #' @param input The Shiny input object.
  #' @param output The Shiny output object.
  #' @param session The Shiny session object.
  #'
  #' @return None. This function updates the Shiny UI output `info_properties` with the result of `dblclick_template_properties`.
  #' @details This function processes double-click events on the entire properties heatmap and updates the `info_properties` UI output.
  dblclick_action_properties_whole = function(df, input,output,session) {
    output[["info_properties"]] = renderUI({
      req(isFALSE(input$properties_sparse_select))
      output <- dblclick_template_properties(gVars,df) 
      output 
    })
  }
  
  
  #' Update outputUI consequent to double-click events for the sparse properties heatmap.
  #' @param df A DFrame object containing `row_label` and `column_label` from the double-click event.
  #' @param input The Shiny input object.
  #' @param output The Shiny output object.
  #' @param session The Shiny session object.
  #'
  #' @return None. This function updates the Shiny UI output (with an ID `info_properties_ht_sparse` concatenated with a numeric index) with the result of `dblclick_template_properties`.
  #' @details This function processes double-click events on the sparse properties heatmap and updates the corresponding UI output.
  dblclick_action_properties_sparse = function(df, input,output,session) {
    output[[paste0("info_properties_ht_sparse",gVars$flagPR)]] = renderUI({
      req(isTRUE(input$properties_sparse_select))
      output <- dblclick_template_properties(gVars,df) 
      output 
    })
  }
  
  
  #' Handle brush events on template properties heatmap 
  #' @param gVars A list containing global variables.
  #' @param df A DFrame containing `row_label` and `column_label` from the brush event.
  #' 
  #' @return A list containing updated global variables and a tag list of UI outputs (a scatterplot and a HTML spacer).
  #' @details This function processes brush events on a properties heatmap. It identifies the brushed area, retrieves the corresponding values from the properties dataframe, and generates a scatterplot showing the mean values.
  brush_template_properties <-function(gVars, df){
    row_index_ids = unique(unlist(df$row_label))
    column_index_ids = unique(unlist(df$column_label))
    properties_df<-gVars$properties_df
    
    if(!is.null(df)) {  
      gVars$brush_area_properties <- reactiveValues(coord= list(row_index_ids, column_index_ids))
      upperlimitPR <- 0.1+ gVars$ceil_max_val_properties 
      sub_dfprop <-  gVars$properties_df[gVars$brush_area_properties$coord[[1]],gVars$brush_area_properties$coord[[2]], drop=FALSE]
      
      col_names <- colnames(sub_dfprop)
      row_names <- rownames(sub_dfprop)
      prop_long_mat <- cbind(sub_dfprop, property= row_names) %>%
        reshape2::melt(., id.vars="property") %>% setNames(., c("property", "endpoint", "score"))
      
      prop_scatterplot <- prop_long_mat %>%
        group_by(endpoint, property) %>%
        summarise(mean_val = mean(score, na.rm = T)) %>%
        ggplot(aes(x = mean_val, y = endpoint, color = property)) +
        geom_point(position=position_jitter(h=0.05,w=0))+ 
        labs(title= "Properties", x= "Mean score",y="Endpoint", color="Properties") + 
        scale_x_continuous(breaks = seq(0,gVars$ceil_max_val_properties ,0.5),
                           limits = c(-0.1,upperlimitPR) )+
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
              plot.title = element_text(hjust = 0.5))
      
      output = tagList()
      output[[1]] <- plotly::renderPlotly({prop_scatterplot})
      output[[2]] <- HTML("<br>")  
      output
    } else {output <-NULL}
    return(list(list(gVars), list(output)))
  }
  
  
  #' Update outputUI consequent to brush events for the whole properties heatmap
  #' @param df A DFrame object containing `row_label` and `column_label` from the brush event.
  #' @param input The Shiny input object.
  #' @param output The Shiny output object.
  #' @param session The Shiny session object.
  #'
  #' @return None. This function updates the Shiny UI output `info_properties` with the result of `brush_template_properties`.
  #' @details This function processes brush events on the entire properties heatmap and updates the `info_properties` UI output.
  brush_action_properties_whole = function(df, input,output,session) {
    output[["info_properties"]] = renderUI({
      req(isFALSE(input$properties_sparse_select))
      out_brush <- brush_template_properties(gVars,df) 
      gVars <- out_brush[[1]][[1]]
      output <- out_brush[[2]]
      output 
    })
  }
  
  
  #' Update outputUI consequent to brush events for the sparse properties heatmap
  #' @param df A DFrame object containing `row_label` and `column_label` from the brush event.
  #' @param input The Shiny input object.
  #' @param output The Shiny output object.
  #' @param session The Shiny session object.
  #'
  #' @return None. This function updates the Shiny UI output (with an ID `info_properties_ht_sparse` concatenated with a numeric index) with the result of `brush_template_properties`.
  #' @details This function processes brush events on the sparse properties heatmap and updates the corresponding UI output.
  brush_action_properties_sparse = function(df, input,output,session) {
    output[[paste0("info_properties_ht_sparse",gVars$flagPR)]] = renderUI({
      req(isTRUE(input$properties_sparse_select))
      out_brush <- brush_template_properties(gVars,df) 
      gVars <- out_brush[[1]][[1]]
      output <- out_brush[[2]]
      output 
    })
  }
  
  
  
  ## modal dialog window for sparse selection of hazard, gene, properties  
  open_sparse_set_modal <- function(dims, ht_id) { 
    modalDialog(
      div(style= "overflow-x: scroll; overflow-y: scroll;",
          fluidRow(br()),
          fluidRow(column(12,  
                          InteractiveComplexHeatmapOutput(heatmap_id=ht_id,
                                                          height1 =dims[1], 
                                                          width1=dims[2], 
                                                          height2 = dims[3], 
                                                          width2=dims[4],
                                                          width3=1000, 
                                                          layout= "1|2|3", compact=FALSE,
                                                          action="dblclick",
                                                          response=c("dblclick", "brush"),
                                                          output_ui = htmlOutput(paste0("info_",ht_id)), 
                                                          output_ui_float=FALSE)
          ))
      ),  
      size="xl", 
      easyClose = FALSE, 
      footer = tagList(
        div(style="width: 100%; display: flex; justify-content: space-between;", 
            modalButton("Dismiss")
        )
      ),
      id="customModal"
    ) 
  }
  
  
  
  # in ui jscode senses closing of modals  
   observeEvent(input$modal_closed, {
     req(!is.null(input$tabs))
     if ((input$tabs == "Hazard Visualization" &  isTRUE(input$hazard_sparse_select)) |
         (input$tabs == "Gene Visualization" &  isTRUE(input$genes_sparse_select)) |
         (input$tabs == "Properties Visualization" &  isTRUE(input$properties_sparse_select))
         ) {
           showNotification(paste("Message:", 
                                  "Please untick 'Sparse select' to display the output of main heatmap selection"),
                            duration = 5, type="message")
       }
   })
  
  
  ## other functions related to modal window 
  #' Clean old modal UI elements 
  #' This function removes UI elements and clears corresponding outputs from a modal window with a specified `id_prefix` and `flag`.
  #'
  #' @param id_prefix A character string. The prefix for the ID of the modal window to clean.
  #' @param flag A character string. A flag to append to the prefix to form the complete ID of the modal window to clean.
  #'
  #' @return This function does not return a value but it modifies the modal window UI.
  #' @details The function constructs the full ID of the modal window. Through 'removeUI', it removes the modal window and sets corresponding output to NULL.
  old_modal_cleaner = function(id_prefix, flag ){
    tmpID<-paste0(id_prefix, flag)
    removeUI(selector=paste0("#",tmpID))
    output[[tmpID]] <- NULL
  } 
  
} #srv
