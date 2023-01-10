# server.R
time <- Sys.time()
cat(file = stderr(), 'Server Go!\n')
#options(shiny.trace=TRUE)
options(shiny.sanitize.errors = FALSE)

# load stuff for server
library(plotly)
library(readr)
library(shiny)
library(shinyjs)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(visNetwork)
library(tibble)
library(ggiraph)
library(RSQLite)
library(DBI)
library(pool)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(shadowtext)
library(htmltools)
library(pool)
library(RSQLite)
library(ggtext)
library(stringr)

# pools for sqlite DBs ------------
gene_pool_2022 <- dbPool(drv = SQLite(), dbname = "./www/2022/eyeIntegration_2022_human.sqlite", idleTimeout = 3600000)
gene_pool_2017 <- dbPool(drv = SQLite(), dbname = "./www/2017/eyeIntegration_human_2017_01.sqlite", idleTimeout = 3600000)
gene_pool_2019 <- dbPool(drv = SQLite(), dbname = "./www/2019/EiaD_human_expression_2019_04.sqlite", idleTimeout = 3600000)
DNTx_pool_2019 <- dbPool(drv = SQLite(), dbname = "./www/2019/EiaD_human_expression_2020_02.DNTx01.sqlite", idleTimeout = 3600000)
SC_pool <- dbPool(drv = SQLite(), dbname = "./www/single_cell_retina_info_04.sqlite", idleTimeout = 3600000)
scEiaD_pool <- dbPool(drv = SQLite(), dbname = "./www/2022/scEiaD.sqlite", idleTimeout = 3600000)

#source('./www/cowplot::theme_cowplot.R')
gene_names_2022 <- gene_pool_2022 %>% tbl('gene_IDs') %>% pull(ID) %>% unique()
gene_names_2017 <- gene_pool_2017 %>% tbl('gene_IDs') %>% pull(ID)
gene_names_2019 <- gene_pool_2019 %>% tbl('gene_IDs') %>% pull(ID)
#geneTX_names_2022 <- gene_pool_2022 %>% tbl('tx_IDs') %>% pull(ID) %>% unique()
geneTX_names_2017 <- gene_pool_2017 %>% tbl('tx_IDs') %>% pull(ID)
geneTX_names_2019 <- gene_pool_2019 %>% tbl('tx_IDs') %>% pull(ID)
geneTX_names_2019_DNTx <- DNTx_pool_2019 %>% tbl('tx_IDs') %>% pull(ID)
core_tight_2022 <- gene_pool_2022 %>% tbl('metadata') %>% as_tibble() %>% select(sample_accession:sample_attribute, region:Comment, Sample_comment, Perturbation)
core_tight_2017 <- gene_pool_2017 %>% tbl('metadata') %>% as_tibble()
core_tight_2019 <- gene_pool_2019 %>% tbl('metadata') %>% as_tibble()

# Data for PCA Visualization - created by the EiaD_build/scripts/pca_workup.Rmd script
load('./www/2022/eye_pca_data.Rdata')
load('./www/2022/eye_percentVar_data.Rdata')

load('./www/2017/retina_module_network_lists.Rdata') # NOTE THESE ARE PRECOMPUTED htmlwidgets 
load('./www/2017/rpe_module_network_lists.Rdata') # NOTE THESE ARE PRECOMPUTED htmlwidgets 
#load('./www/go_heatmap.Rdata')
load('./www/basic_stats.Rdata')
cat(file=stderr(), 'Data loaded in ')
cat(file=stderr(), Sys.time() - time)
cat(file=stderr(), ' seconds.\n')

onStop(function() {
  poolClose(gene_pool_2022)
})
onStop(function() {
  poolClose(gene_pool_2017)
})
onStop(function() {
  poolClose(gene_pool_2019)
})
onStop(function() {
  poolClose(SC_pool)
})
onStop(function() {
  poolClose(scEiaD_pool)
})

core_tight_2017$sample_accession<-gsub('E-MTAB-','E.MTAB.',core_tight_2017$sample_accession)
core_tight_2017$Sub_Tissue <- gsub('_',' - ',core_tight_2017$Sub_Tissue)
core_tight_2019$sample_accession<-gsub('E-MTAB-','E.MTAB.',core_tight_2019$sample_accession)
core_tight_2019$Sub_Tissue <- gsub('_',' - ',core_tight_2019$Sub_Tissue)
core_tight_2022$sample_accession<-gsub('E-MTAB-','E.MTAB.',core_tight_2022$sample_accession)
core_tight_2022$Sub_Tissue <- gsub('_',' - ',core_tight_2022$Sub_Tissue)

metasc <- scEiaD_pool %>% tbl("scEiaD_CT_table") %>% select(CellType_predict) %>% as_tibble() %>% unique()
CellType_predict_val <- setNames(c(pals::glasbey(n = 32), pals::kelly(n = scEiaD_pool %>% tbl('cell_types') %>% pull(1) %>% length() - 32)) %>% colorspace::lighten(0.3), metasc %>% pull(CellType_predict) %>% sort())
CellType_predict_col <- scale_colour_manual(values = CellType_predict_val)
CellType_predict_fill <- scale_fill_manual(values = CellType_predict_val)

# make global tissue vals
tissues <- c(core_tight_2017$Tissue, core_tight_2019$Tissue, core_tight_2022$Tissue)%>% unique() %>% sort() 
tissue_val <- setNames(c( pals::kelly(n = tissues %>% length() - 32), pals::glasbey(n = 32)) %>% colorspace::darken(0.2), tissues %>% sort())



# site begins! ---------
shinyServer(function(input, output, session) {
  # Observe: update fields for pan tissue plots --------
  # also handles custom url
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query[['Dataset']])) {
      updateTextInput(session, "Database", value = gsub('_', ' ', as.character(query['Dataset'])))
    }
    db = input$Database # c("Gene 2017", "Gene 2019", "Transcript 2017", "Transcript 2019", "Gene 2022")
    if (db == 'Gene 2017'){ID_names = gene_names_2017 %>% sort()
    } else if (db == 'Gene 2019'){ID_names = gene_names_2019 %>% sort()
    } else if (db == 'Transcript 2017'){ID_names = geneTX_names_2017 %>% sort()
    } else if (db == 'Gene 2022'){ID_names = gene_names_2022 %>% sort()
    } else if (db == 'DNTx v01'){ID_names = geneTX_names_2019_DNTx %>% sort()
    } else if (db == 'scEiaD_pool'){ID_names = Gene %>% sort()
    } else {ID_names = geneTX_names_2019 %>% sort()}
    # gene / tx lists
    if (is.null(query[['scmaturity']])){
      updateSelectizeInput(session, 'scmaturity',
                           choices = c("Early Development", "Maturing", "Mature", "Late Development"),
                           options = list(placeholder = 'Type to search'),
                           server = TRUE)
    }
    if (is.null(query[['scGene']])){
      updateSelectizeInput(session, 'scGene',
                           choices = scEiaD_pool %>% tbl("scEiaD_CT_table") %>% pull(Gene) %>% unique(),
                           options = list(placeholder = 'Type to search'),
                           server = TRUE)
    }
    if (is.null(query[['scplot_tissue_gene']])){
      updateSelectizeInput(session, 'scplot_tissue_gene',
                           choices = scEiaD_pool %>% tbl("scEiaD_CT_table") %>% pull(CellType_predict) %>% unique(),
                           options = list(placeholder = 'Type to search'),
                           server = TRUE)
    }
    if (is.null(query[['ID']])){
      updateSelectizeInput(session, 'ID',
                           choices = ID_names,
                           options = list(placeholder = 'Type to search'),
                           server = TRUE)
    }
    if (!is.null(query[['ID']])) {
      select_gene <- strsplit(toupper(as.character(query['ID'])),split = ',')[[1]] %>% 
        gsub(pattern = '_', ' ', .)
      #select_gene <- c(select_gene, input$ID)
      updateSelectizeInput(session, 'ID',
                           choices = ID_names,
                           selected = select_gene,
                           server = TRUE)
    }
    if (is.null(query[['Tissue']])){
      # tissue choices
      if (grepl('2017', db)){tissues <- unique(sort(core_tight_2017$Sub_Tissue))}
      else if (grepl('2022', db)){tissues <- unique(sort(core_tight_2022 %>% filter(Tissue != 'EyeLid') %>% pull(Tissue)))
      } else {tissues <- unique(sort(core_tight_2019 %>% filter(Sub_Tissue != 'Choroid Plexus - Adult Tissue') %>% pull(Sub_Tissue)))}
      updateSelectizeInput(session, 'plot_tissue_gene',
                           choices= tissues,
                           options = list(placeholder = 'Type to search'),
                           server = TRUE)
    }
    if (!is.null(query[['Tissue']])){
      select_tissue <- gsub(pattern = '_',replacement = ' ', x = as.character(query['Tissue']))
      select_tissue <- strsplit(select_tissue, split = ',')[[1]]
      #select_tissue <- c(select_tissue, input$plot_tissue_gene)
      if (grepl('2017', db)){tissues <- unique(sort(core_tight_2017$Sub_Tissue))}
      else if (grepl('2022', db)){tissues <- unique(sort(core_tight_2022 %>% filter(Sub_Tissue != 'Choroid Plexus' & Sub_Tissue != 'WIBR3 hESC Choroid plexus Organoids') %>% pull(Sub_Tissue)))
      } else {tissues <- unique(sort(core_tight_2019 %>% filter(Sub_Tissue != 'Choroid Plexus - Adult Tissue') %>% pull(Sub_Tissue)))}
      updateSelectizeInput(session, 'plot_tissue_gene',
                           choices= tissues,
                           selected = select_tissue,
                           server = TRUE)
    }
  })
  # Observe: differential expression -----
  # update database (2017 and 2019) based on user input
  observe({
    db = input$diff_database # c("Gene 2017", "Gene 2019", "Transcript 2019")
    if (grepl('2017', db)){
      pool <- 'gene_pool_2017'
      ids <- 'gene_IDs'
    } else if (db == 'Gene 2019') {
      pool <- 'gene_pool_2019'
      ids <- 'gene_IDs'
    } else if (db == 'Transcript 2019'){
      pool <- 'gene_pool_2019'
      ids <- 'tx_ids'
    } else if (db == 'DNTx v01'){
      pool <- 'DNTx_pool_2019'
      ids <- 'tx_ids'
    } else if (db == 'scEiaD_pool'){
      pool <- 'scEiaD_pool'
      ids <- 'gene_IDs'
    }
    
    de_comparison_contrast_names <- get(pool) %>% tbl('limma_DE_tests') %>% pull(Comparison)
    names(de_comparison_contrast_names) <- get(pool) %>% tbl('limma_DE_tests') %>% pull(Name)
    de_comparison_contrast_names <- sort(de_comparison_contrast_names)
    gene_tx_types = get(pool) %>% tbl(ids) %>% pull(2) %>% unique() %>% sort()
    
    updateSelectizeInput(session, 'gene_tx_type',
                         choices = gene_tx_types,
                         selected = 'protein_coding',
                         server = TRUE)
    
    # DE tests
    updateSelectizeInput(session, 'de_comparison',
                         choices = de_comparison_contrast_names,
                         options = list(placeholder = 'Type to search'),
                         server = TRUE)
    
  })
  
  # Observe: Gene Expression -> Mouse SC Expression observe -----
  # change mouse gene name list depending on what dataset is selected (Clark or Macosko)
  # also change the cell types
  observe({
    #req(input$SC_dataset)
    SC_dataset <- (input$SC_dataset %>% strsplit(' '))[[1]][1] %>% tolower()
    table_name <- paste(SC_dataset, 'gene_names', sep='__')
    metadata_name <- paste(SC_dataset, 'SC_metadata_long', sep = '__')
    mouse_gene_names <- dbGetQuery(SC_pool, paste('SELECT * FROM ', table_name)) %>% pull(Gene)
    SC_cell_types <- dbGetQuery(SC_pool, paste('SELECT * from ', metadata_name)) %>% pull(`Cell Type`) %>% unique() %>% sort()
    SC_cell_types <- SC_cell_types[!grepl('Red|Doub', SC_cell_types)]
    
    updateSelectizeInput(session, 'mGene',
                         choices = sort(mouse_gene_names),
                         selected = c(mouse_gene_names[grepl('nrl', mouse_gene_names, ignore.case = T)],
                                      mouse_gene_names[grepl('crx$', mouse_gene_names, ignore.case = T)]),
                         server = TRUE)
    updateSelectizeInput(session, 'SC_cell_types',
                         choices = SC_cell_types,
                         selected = c(SC_cell_types[1], SC_cell_types[3]),
                         server = TRUE)
  })
  
  # Observe: 2D clustering -> Mouse SC observe -----
  # like just above, change mouse gene name list depending on what dataset is selected (Clark or Macosko)
  # also change the ages displayed (if multiple available)
  observe({
    SC_dataset <- (input$SC_dataset_tsne %>% strsplit(' '))[[1]][1] %>% tolower()
    table_name <- paste(SC_dataset, 'gene_names', sep='__')
    metadata_name <- paste(SC_dataset, 'SC_metadata_long', sep = '__')
    mouse_gene_names <- dbGetQuery(SC_pool, paste('SELECT * FROM ', table_name)) %>% pull(Gene)
    if (SC_dataset == 'macosko'){
      tsne_age <- c('P14')
    } else {
      tsne_age <- dbGetQuery(SC_pool, paste('SELECT * from ', metadata_name)) %>% pull(Age) %>% unique() %>% sort()
    }
    
    updateSelectizeInput(session, 'mGene_tsne',
                         choices = sort(mouse_gene_names),
                         selected = c(mouse_gene_names[grepl('abca4$', mouse_gene_names, ignore.case = T)],
                                      mouse_gene_names[grepl('crx$', mouse_gene_names, ignore.case = T)]),
                         server = TRUE)
    updateSelectizeInput(session, 'age_tsne',
                         choices = tsne_age,
                         selected = c(tsne_age[3]),
                         server = TRUE)
  })
  
  # Observe: update SC datatable  -------
  # sort of like above, change tissue list depending on what dataset is selected (Clark or Macosko)
  observe({
    SC_dataset <- (input$SC_datatable_dataset %>% strsplit(' '))[[1]][1] %>% tolower()
    table_name <- paste(SC_dataset, 'SC_metadata_long', sep='__')
    metadata_name <- paste(SC_dataset, 'SC_metadata_long', sep = '__')
    SC_cell_types <- dbGetQuery(SC_pool, paste('SELECT * from ', metadata_name)) %>% pull(`Cell Type`) %>% unique() %>% sort()
    SC_cell_types <- SC_cell_types[!grepl('Red|Doub', SC_cell_types)]
    
    updateSelectizeInput(session, 'sc_datatable_tissue',
                         choices = SC_cell_types,
                         selected = c(SC_cell_types[grepl('Rod',SC_cell_types, ignore.case = T)],
                                      SC_cell_types[grepl('Cone$', SC_cell_types, ignore.case = T)]),
                         server = TRUE)
  })
  
  # Observe: update data table inputs ----------
  observe({
    table_db <- input$table_db
    if (grepl('2017', table_db)){
      updateSelectizeInput(session, 'table_tissue',
                           choices= unique(sort(core_tight_2017$Tissue)),
                           selected= 'Retina',
                           server = TRUE)
      updateSelectizeInput(session, 'table_columns',
                           choices = sort(colnames(core_tight_2017)),
                           selected = colnames(core_tight_2017) %>%
                             grep("study_abstract|sample_attribute|Kept", ., value = T, invert = T),
                           server = TRUE)
      if (grepl('Gene', table_db)){
        updateSelectizeInput(session, 'table_gene',
                             choices = gene_names_2017,
                             selected= 'TYRP1',
                             server = TRUE)
      } else {
        updateSelectizeInput(session, 'table_gene',
                             choices = geneTX_names_2017,
                             selected= 'TYRP1 (ENST00000381136.2)',
                             server = TRUE)
      }
    }
    
    else if (grepl('2022', table_db)){
      updateSelectizeInput(session, 'table_tissue',
                           choices= unique(sort(core_tight_2022$Tissue)),
                           selected= 'Retina',
                           server = TRUE)
      updateSelectizeInput(session, 'table_columns',
                           choices = sort(colnames(core_tight_2022)),
                           selected = colnames(core_tight_2022) %>%
                             grep("study_abstract|sample_attribute|Kept", ., value = T, invert = T),
                           server = TRUE)
      if (grepl('Gene', table_db)){
        updateSelectizeInput(session, 'table_gene',
                             choices = gene_names_2022,
                             selected= 'TYRP1',
                             server = TRUE)
      } else {
        updateSelectizeInput(session, 'table_gene',
                             choices = geneTX_names_2022,
                             selected= 'TYRP1 (ENST00000381136.2)',
                             server = TRUE)
      }
    }
    
    else {
      updateSelectizeInput(session, 'table_columns',
                           choices = sort(colnames(core_tight_2019)),
                           selected = colnames(core_tight_2019) %>%
                             grep("study_abstract|sample_attribute|Kept", ., value = T, invert = T),
                           server = TRUE)
      updateSelectizeInput(session, 'table_tissue',
                           choices= unique(sort(core_tight_2019$Tissue)),
                           selected= 'Retina',
                           server = TRUE)
      if (grepl('Gene', table_db)){
        updateSelectizeInput(session, 'table_gene',
                             choices = gene_names_2019,
                             selected= 'TYRP1',
                             server = TRUE)
      } else {
        updateSelectizeInput(session, 'table_gene',
                             choices = geneTX_names_2019,
                             selected= 'TYRP1 (ENST00000381136.2)',
                             server = TRUE)
      }
    }
  })
  
  observe({
    table_db <- input$temporal_retina_heatmap_table
    if (grepl('Transcript', table_db)){
      updateSelectizeInput(session, 'temporal_retina_heatmap_ID',
                           choices = geneTX_names_2019,
                           selected= c('OTX2 (ENST00000339475.9)','OTX2 (ENST00000408990.7)',
                                       'CRX (ENST00000539067.5)','CRX (ENST00000602001.1)',
                                       'CRX (ENST00000613299.1)', 'CRX (ENST00000221996.11)'),
                           server = TRUE) }
    else {
      updateSelectizeInput(session, 'temporal_retina_heatmap_ID',
                           choices = gene_names_2019,
                           selected= c('OTX2','NRL'),
                           server = TRUE)
    }
  })
  
  updateSelectizeInput(session, 'FaF_ID',
                       choices = gene_names_2019,
                       selected= 'OTX2',
                       server = TRUE)
  
  updateSelectizeInput(session, 'retina_gene',
                       choices = gene_pool_2017 %>% tbl('retina_gene_name_colors') %>%  pull(id),
                       selected='ABCA4',
                       server = TRUE)
  updateSelectizeInput(session, 'retina_gene_edges',
                       choices= gene_pool_2017 %>% tbl('retina_gene_name_colors') %>%  pull(id),
                       selected='ABCA4',
                       server = TRUE)
  updateSelectizeInput(session, 'rpe_gene',
                       choices = gene_pool_2017 %>% tbl('rpe_gene_name_colors') %>% pull(id),
                       selected='RPL24',
                       server = TRUE)
  updateSelectizeInput(session, 'rpe_gene_edges',
                       choices= gene_pool_2017 %>% tbl('rpe_gene_name_colors') %>% pull(id),
                       selected='RPL24',
                       server = TRUE)
  
  # Observe: GENE update tissues for fold change -----
  observe({
    selected_tissue <- input$plot_tissue_gene
    updateSelectizeInput(session, 'Bench_gene',
                         choices= selected_tissue,
                         selected = selected_tissue,
                         server = TRUE)
    
  })
  
  # Pan - Eye PCA -------
  
  visualize_pca_function <- eventReactive(input$pca_button, {
    pcFirst <- input$pca_component_one
    pcSecond <- input$pca_component_two
    
    validate(
      need(input$pca_component_one != input$pca_component_two, 
           "Please select two distinct PCA components and click the (RE)Draw PCA Plot! button. It may take a few seconds for the plot to appear.")
    )
    
      p <- eye_pca_data %>%
        #filter(Source != 'scRNA') %>%
        as_tibble() %>%
        ggplot(., aes(.data[[pcFirst]], .data[[pcSecond]])) +
        geom_point(size=3, aes(color=Tissue, shape = Source,
                               text = paste("Study: ", study_accession, "\n", "Sample: ", sample_accession, "\n",
                                            "Tissue: ", Tissue, "\n", "Sub-Tissue: ", Sub_Tissue, "\n", "Source: ", 
                                            Source, "\n", "Age: ", Age, "\n", "Count: ", 
                                            Count, "\n", "Nearest Neighboring Tissues: ", Tissue2))) +
        xlab(paste0(pcFirst, ": ",eye_percentVar_data[str_extract(pcFirst, '\\d+') %>% as.integer()],"% variance")) +
        ylab(paste0(pcSecond, ": ",eye_percentVar_data[str_extract(pcSecond, '\\d+') %>% as.integer()],"% variance")) +
        cowplot::theme_cowplot() +
        ggtitle(label = "Ocular Sample PCA Visualization for the top 1000 protein coding genes in eyeIntegration") +
        scale_color_manual(values = c(pals::glasbey(), pals::alphabet2(), pals::alphabet2()) %>% unname()) +
        scale_fill_manual(values = c(pals::glasbey(), pals::alphabet2(), pals::alphabet2()) %>% unname()) +
        scale_shape_manual(values = 0:10)
      
      ggplotly(p, tooltip = 'text')
    })
    
    output$eye_pca_plot <- renderPlotly({
      visualize_pca_function()
    })
  
  # Pan - Tissue Boxplot -------
  boxPlot_gene_func <- eventReactive(input$pan_button_gene, {
    cat(file=stderr(), 'boxPlot Gene call\n')
    db <- input$Database
    gene <- input$ID
    tissue <- input$plot_tissue_gene
    col_num <- input$num_gene
    if (length(db) < 1 || length(gene) < 1 || length(tissue) < 1){
      showModal(modalDialog(title = "Box Plot Warning",
                            "Have you specified at least one gene or tissue?",
                            easyClose = T,
                            footer = NULL))
    }
    if (db == 'Gene 2017') {
      query = paste0('select * from lsTPM_gene where ID in ("',paste(gene, collapse='","'),'")')
      p <- dbGetQuery(gene_pool_2017, query) %>%
        left_join(.,core_tight_2017) %>%
        left_join(., gene_pool_2017 %>% tbl('gene_IDs') %>% as_tibble()) %>%
        as_tibble()
    } else if (db == 'Transcript 2017'){
      query = paste0('select * from lsTPM_TX where ID in ("',paste(gene, collapse='","'),'")')
      p <- dbGetQuery(gene_pool_2017, query) %>% left_join(.,core_tight_2017) %>%
        left_join(., gene_pool_2017 %>% tbl('tx_IDs') %>% as_tibble()) %>%
        as_tibble()
    } else if (db == 'Gene 2019'){
      query = paste0('select * from lsTPM_gene where ID in ("',paste(gene, collapse='","'),'")')
      p <- dbGetQuery(gene_pool_2019, query) %>% left_join(.,core_tight_2019) %>%
        left_join(., gene_pool_2019 %>% tbl('gene_IDs') %>% as_tibble()) %>%
        as_tibble()
    } else if (db == 'Transcript 2019'){
      query = paste0('select * from lsTPM_tx where ID in ("',paste(gene, collapse='","'),'")')
      p <- dbGetQuery(gene_pool_2019, query) %>% left_join(.,core_tight_2019) %>%
        left_join(., gene_pool_2019 %>% tbl('tx_IDs') %>% as_tibble()) %>%
        as_tibble()
    } else if (db == 'DNTx v01'){
      query = paste0('select * from lsTPM_tx where ID in ("',paste(gene, collapse='","'),'")')
      p <- dbGetQuery(DNTx_pool_2019, query) %>% left_join(.,core_tight_2019) %>%
        left_join(., DNTx_pool_2019 %>% tbl('tx_IDs') %>% as_tibble()) %>%
        as_tibble()
    } else if (db == 'Gene 2022'){
      query = paste0('select * from lsTPM_gene where ID in ("',paste(gene, collapse='","'),'")')
      p <- dbGetQuery(gene_pool_2022, query) %>% left_join(.,core_tight_2022) %>%
        left_join(., gene_pool_2022 %>% tbl('gene_IDs') %>% as_tibble()) %>%
        as_tibble()
    }
    p$Type <- p %>% select(contains('type')) %>% pull(1)
    
    if (grepl('2022', db)){
      plot_data <- p %>%
        filter(Tissue %in% tissue) 
    } else {
      plot_data <- p %>%
        filter(Sub_Tissue %in% tissue) 
    }
    
    # fix tissue <-> color
    tissue_val <- tissue_val[plot_data$Tissue %>% unique()]
    tissue_col <- scale_colour_manual(values = tissue_val)
    tissue_fill <- scale_fill_manual(values = tissue_val)
    
    if (!grepl('2022', db)){
      p <- plot_data %>%
        mutate(Info = paste('SRA: ',
                            sample_accession,
                            '\nStudy: ',
                            study_title, '\n',
                            gsub('\\|\\|', '\n',
                                 sample_attribute),
                            sep =''),
               ID = gsub(' \\(', '\n(', ID)) %>%
        ggplot(data=.,aes(x=Sub_Tissue,y=log2(value+1), color = Tissue, fill = Tissue)) +
        #geom_violin(alpha=0.5, scale = 'width') +
        geom_boxplot(alpha=0.7, outlier.shape = NA, width = 0.6, fill = 'black', color = 'white') +
        geom_point_interactive(size=1, position = 'jitter', alpha=0.25, stroke = 3, aes(tooltip=htmlEscape(Info, TRUE), shape = Type)) +
        xlab('') +
        facet_wrap(~ID, ncol=col_num, scales = 'free_x') +
        cowplot::theme_cowplot(font_size = 15) + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.2)) +
        ggtitle('Box Plot of Pan-Human Gene Expression') +
        ylab("log2(TPM + 1)") +
        scale_shape_manual(values=c(0:2,5,6,15:50)) +
        theme(plot.margin=grid::unit(c(0,0,0,0.1), "cm"),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.key.size= unit(0.2, "cm"),
              legend.spacing = unit(0.2, "cm")) +
        tissue_col + tissue_fill
    } else {
      p <- plot_data %>% 
        #mutate(Perturbation = case_when(grepl('MGS', Source_details) ~ Source_details)) %>% 
        mutate(Sub_Tissue = case_when(is.na(Sub_Tissue) ~ '', TRUE ~ Sub_Tissue), 
               Source = case_when(is.na(Source) ~ '', TRUE ~ Source), 
               Age = case_when(is.na(Age) ~ '', TRUE ~ Age),
               Perturbation = case_when(is.na(Perturbation) ~ '', TRUE ~ Perturbation)) %>% 
        mutate(Sub_Tissue = glue::glue("<span style='color:#E41A1C'>{Sub_Tissue}</span>"),
               Source = glue::glue("<span style='color:#377EB8'>{Source}</span>"),
               Age = glue::glue("<span style='color:#4DAF4A'>{Age}</span>"),
               Perturbation = glue::glue("<span style='color:#984EA3'>{Perturbation}</span>")
        ) %>% 
        mutate(Info = paste('SRA: ',
                            sample_accession,
                            '\nStudy: ',
                            study_title, '\n',
                            gsub('\\|\\|', '\n',
                                 sample_attribute),
                            sep =''),
               ID = gsub(' \\(', '\n(', ID)) %>%
        ggplot(data=.,aes(x=interaction(Source, Sub_Tissue, Age, Perturbation, sep = ' | '),y=log2(value+1), 
                          color = Tissue, 
                          fill = Tissue)) +
        #geom_violin(alpha=0.5, scale = 'width') +
        geom_boxplot(alpha=0.7, outlier.shape = NA, width = 0.6, fill = 'black', color = 'white') +
        cowplot::theme_cowplot(font_size = 15) + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.2)) +
        ggtitle('Box Plot of Pan-Human Gene Expression') +
        ylab("log2(TPM + 1)") +
        scale_shape_manual(values=c(0:2,5,6,15:50)) +
        theme(strip.background = element_rect(fill = 'black'),
              strip.text = element_text(color = 'white'),
              panel.background = element_rect(fill = 'gray80'),
              plot.margin=grid::unit(c(0,0,0,0.1), "cm"),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.key.size= unit(0.2, "cm"),
              legend.spacing = unit(0.2, "cm"))  +
        tissue_col + 
        tissue_fill 
      
      if (input$rotation == 1){
        p <- p + 
          coord_flip() + 
          facet_grid(rows = vars(Tissue), cols = vars(ID), scales = 'free_y', space = 'free') +
          theme(strip.text.y.right = element_text(angle = 0)) +
          theme(
            axis.text.y = element_markdown(),
            axis.title.y = element_markdown()) +
          labs(x = "<span style='color:#377EB8'>Source</span> | 
       <span style='color:#E41A1C'>Sub Tissue</span> |
       <span style='color:#4DAF4A'>Age</span> |
       <span style='color:#984EA3'>Perturbation</span>")
      } else {
        p <- p + 
          facet_grid(cols = vars(Tissue), rows = vars(ID), 
                     scales = 'free_x', space = 
                       'free', labeller = labeller(Tissue = label_wrap_gen(7))) +
          theme(
            axis.text.x = element_markdown(),
            axis.title.x = element_markdown()) +
          labs(x = "<span style='color:#377EB8'>Source</span> | 
       <span style='color:#E41A1C'>Sub Tissue</span> |
       <span style='color:#4DAF4A'>Age</span> |
       <span style='color:#984EA3'>Perturbation</span>")
      }
      
      if (input$points){
        p <- p + geom_point_interactive(size=1, position = 'jitter', alpha=0.25, stroke = 3, aes(tooltip=htmlEscape(Info, TRUE), shape = Type)) 
      }
      
    }
    output <- list()
    if (!grepl('2022',db)){
      output$plot <- girafe(ggobj = p,
                            width_svg = 14,
                            height_svg= max(12, (6 * (length(gene)/min(col_num,length(gene)))))) %>%
        girafe_options(., opts_toolbar(position = "top") )
    } else {
      output$plot <- girafe(ggobj = p,
                            width_svg = 16,
                            height_svg= max(12, (2 * (length(tissue)/min(col_num,length(tissue)))))) %>%
        girafe_options(., opts_toolbar(position = "top") )
    }
    output$data <- plot_data
    output
    
  })
  output$boxPlot_gene <- renderGirafe({
    boxPlot_gene_func()$plot
  })
  
  # output plot data -------
  output$plot_table <- DT::renderDataTable({
    boxPlot_gene_func()$data %>%  DT::datatable(extensions = 'Buttons',
                                                filter = list(position = 'top', clear = TRUE, plain = TRUE),
                                                options = list(pageLength = 10, scrollX = T, searchHighlight = TRUE, dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
  }, server=FALSE)
  
  
  # Gene heatmap -------
  bulk_tissue_heatmap_func <- eventReactive(input$pan_button_gene, {
    #cat(file=stderr(), 'Gene heatmap call\n')
    db = input$Database
    gene <- input$ID
    tissue <- input$plot_tissue_gene
    if (length(db) < 1 || length(gene) < 2 || length(tissue) < 2){
      showModal(modalDialog(title = "Heatmap Warning",
                            "Have you specified at least two genes and two tissues?",
                            easyClose = T,
                            footer = NULL))
    }
    label_size = 0.2
    if (db == 'Gene 2017') {
      pool <- 'gene_pool_2017'
      table <- 'mean_rank_decile_gene'
      meta <- 'core_tight_2017'
    } else if (db == 'Transcript 2017'){
      pool <- 'gene_pool_2017'
      table <- 'mean_rank_decile_tx'
      meta <- 'core_tight_2017'
      label_size = 0.8
    } else if (db == 'Gene 2022'){
      pool <- 'gene_pool_2022'
      table <- 'mean_rank_decile_gene'
      meta <- 'core_tight_2022'
      label_size = 0.8
    } else if (db == 'Gene 2019') {
      tissue <- trimws(tissue)
      pool <- 'gene_pool_2019'
      meta <- 'core_tight_2019'
      table <- 'mean_rank_decile_gene'
    } else if (db == 'Transcript 2019'){
      tissue <- trimws(tissue)
      pool <- 'gene_pool_2019'
      meta <- 'core_tight_2019'
      table <- 'mean_rank_decile_tx'
      label_size = 0.8
    } else if (db == 'DNTx v01'){
      tissue <- trimws(tissue)
      pool <- 'DNTx_pool_2019'
      meta <- 'core_tight_2019'
      table <- 'mean_rank_decile_tx'
      label_size = 0.8
    }
    
    id_matrix <- get(pool) %>%
      tbl(table) %>%
      filter(ID %in% gene, Sub_Tissue %in% tissue) %>%
      data.frame() %>%
      mutate(lsTPM = log2(meanlsTPM+1)) %>%
      select(-meanlsTPM, -Rank, -Decile) %>%
      spread(Sub_Tissue, lsTPM)
    gene_IDs <- id_matrix %>% pull(ID)
    id_matrix <- id_matrix %>% select(-ID) %>% as.matrix()
    
    row.names(id_matrix) <- gene_IDs
    text_col <- NA
    
    ha = HeatmapAnnotation(Tissue = get(meta) %>%
                             select(Tissue, Sub_Tissue) %>%
                             unique() %>%
                             mutate(Sub_Tissue = trimws(Sub_Tissue)) %>%
                             right_join(colnames(id_matrix) %>%
                                          trimws() %>%
                                          enframe(),
                                        by = c('Sub_Tissue' = 'value')) %>%
                             pull(Tissue),
                           col = list(Tissue = tissue_val),
                           show_annotation_name = TRUE,
                           which = 'column')
    
    breaks = c(0,5,10)
    show_row_names = TRUE
    if (1 %in% input$heatmap_clustering_checkbox){
      cluster_rows = TRUE
    } else {cluster_rows = FALSE}
    if (2 %in% input$heatmap_clustering_checkbox){
      cluster_cols = TRUE
    } else {cluster_cols = FALSE}
    
    Heatmap(id_matrix,
            top_annotation = ha,
            cluster_columns = cluster_cols,
            #column_title = title,
            cluster_rows = cluster_rows,
            col = colorRamp2(breaks = breaks, colors = viridis(length(breaks))),
            rect_gp = gpar(col= "white"),
            show_row_names = show_row_names,
            name = 'log2(TPM+1)',
            #show_heatmap_legend = show_heatmap_legend,
            column_names_max_height = unit(8, "cm"),
            row_names_max_width = unit(8, "cm"),
            clustering_distance_rows = "pearson",
            clustering_distance_columns = "euclidean")
    
    
    
  })
  output$bulk_tissue_heatmap <- renderPlot({
    bulk_tissue_heatmap_func()
  },  height=eventReactive(input$pan_button_gene, {max(500, (40*length(input$ID))/min(input$num_gene,length(input$ID)))}))
  
  
  
  # Gene info table ---------
  gene_info_maker <- eventReactive(input$pan_button_gene, {
    gene <- gsub(' (.*)','', input$ID) %>% unique()
    ensembl_base <- "https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g="
    genecards_base <- "https://www.genecards.org/cgi-bin/carddisp.pl?gene="
    omim_base <- "https://www.omim.org/search/?index=entry&sort=score+desc%2C+prefix_sort+desc&start=1&limit=10&search="
    gene_info <- gene %>% as_tibble() %>% rename(ID = value) %>%
      mutate(Ensembl = paste0('<a href="', ensembl_base, ID, '", target="_blank">Ensembl</a>'),
             GeneCards = paste0('<a href="', genecards_base, ID, '", target="_blank">GeneCards</a>'),
             OMIM = paste0('<a href="', omim_base, ID, '", target="_blank">OMIM</a>'))
    
    gene_info %>% DT::datatable(rownames = F, escape = F,
                                options = list(pageLength = 1000, dom = ''))
    
  })
  output$gene_info <- DT::renderDataTable({
    gene_info_maker()
  })
  
  # Gene boxplot stats ----
  rankStats_gene_func <- eventReactive(input$pan_button_gene, {
    cat(file=stderr(), 'Gene rankStats call\n')
    db <- input$Database
    gene <- input$ID
    tissue <- input$plot_tissue_gene
    if (db == 'Gene 2017'){
      pool <- 'gene_pool_2017'
      table <- 'mean_rank_decile_gene'
    } else if (db == 'Transcript 2017'){
      pool <- 'gene_pool_2017'
      table <- 'mean_rank_decile_tx'
    } else if (db == 'Gene 2022'){
      pool <- 'gene_pool_2022'
      table <- 'mean_rank_decile_gene'
    } else if (db == 'Gene 2019'){
      tissue <- trimws(tissue)
      pool <- 'gene_pool_2019'
      table <- 'mean_rank_decile_gene'
    } else if (db == 'Transcript 2019'){
      tissue <- trimws(tissue)
      pool <- 'gene_pool_2019'
      table <- 'mean_rank_decile_tx'
    } else if (db == 'DNTx v01'){
      tissue <- trimws(tissue)
      pool <- 'DNTx_pool_2019'
      table <- 'mean_rank_decile_tx'
    }
    get(pool) %>% tbl(table) %>%
      filter(ID %in% gene, Sub_Tissue %in% tissue) %>%
      mutate(ID, Tissue = `Sub_Tissue`) %>%
      ungroup() %>%
      dplyr::select(ID, Tissue, meanlsTPM, Rank, Decile) %>%
      arrange(ID, Tissue) %>%
      as_tibble() %>%
      mutate(`log2(TPM + 1)` = log2(meanlsTPM + 1)) %>%
      dplyr::select(-meanlsTPM) %>%
      DT::datatable(extensions = 'Buttons', rownames = F, options = list(
        pageLength = 20, dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatRound(c('log2(TPM + 1)'), digits=2)
  })
  
  output$rankStats_gene <- DT::renderDataTable(server = TRUE, {
    rankStats_gene_func()
  })
  
  # FC table stats ----------
  output$basicStats_gene <- DT::renderDataTable(server = TRUE, {
    input$pan_button_gene
    isolate({
      db <- input$Database
      gene <- input$ID
      tissue <- input$plot_tissue_gene
      bench <- input$Bench_gene
    })
    if (db == 'Gene 2017'){
      query = paste0('select * from lsTPM_gene where ID in ("',paste(gene, collapse='","'),'")')
      plot_data <- dbGetQuery(gene_pool_2017,  query) %>%
        left_join(.,core_tight_2017)
    } else if (db == 'Transcript 2017'){
      query = paste0('select * from lsTPM_tx where ID in ("',paste(gene, collapse='","'),'")')
      plot_data <- dbGetQuery(gene_pool_2017,  query) %>%
        left_join(.,core_tight_2017)
    } else if (db == 'Gene 2022'){
      query = paste0('select * from lsTPM_gene where ID in ("',paste(gene, collapse='","'),'")')
      plot_data <- dbGetQuery(gene_pool_2022,  query) %>%
        left_join(.,core_tight_2022)
    } else if (db == 'Gene 2019'){
      query = paste0('select * from lsTPM_gene where ID in ("',paste(gene, collapse='","'),'")')
      plot_data <- dbGetQuery(gene_pool_2019,  query) %>%
        left_join(.,core_tight_2019)
    } else if (db == 'Transcript 2019'){
      query = paste0('select * from lsTPM_tx where ID in ("',paste(gene, collapse='","'),'")')
      plot_data <- dbGetQuery(gene_pool_2019,  query) %>%
        left_join(.,core_tight_2019)
    } else if (db == 'DNTx v01'){
      query = paste0('select * from lsTPM_tx where ID in ("',paste(gene, collapse='","'),'")')
      plot_data <- dbGetQuery(DNTx_pool_2019,  query) %>%
        left_join(.,core_tight_2019)
    }
    base_stats <- plot_data %>%
      filter(Sub_Tissue %in% tissue) %>%
      group_by(ID) %>%
      mutate(Bench=ifelse(Sub_Tissue %in% bench, 1, 0), BenchValue=mean(log2(value[Bench==1]+1))) %>%
      group_by(ID, Sub_Tissue) %>%
      summarise(log2DeltaFC=mean(log2(value+1)) - mean(BenchValue), mean=mean(log2(value+1)))
    
    # does t.test against a user-defined reference
    # corrects for number of tests
    tissue_subset <- plot_data %>% filter(Sub_Tissue %in% bench)
    pvals <- plot_data %>%
      group_by(Sub_Tissue, ID) %>%
      do(broom::tidy(t.test(.$value, tissue_subset$value))) %>%
      # multiple test correction
      mutate(`t test p` = signif(min(1,p.value * length(unique(plot_data$Sub_Tissue)))),3) %>%
      dplyr::select(ID, Sub_Tissue, `t test p`)
    stat_join <- left_join(base_stats, pvals) %>%
      mutate(`Gene Name` = ID, Tissue = `Sub_Tissue`, `log2 Fold Change` = log2DeltaFC, `Fold Change` = 2^log2DeltaFC, `Mean Expression` = mean) %>%
      ungroup() %>%
      dplyr::select(`Gene Name`, Tissue, `log2 Fold Change`, `Fold Change`, `Mean Expression`, `t test p`) %>%
      arrange(`Gene Name`, Tissue) %>%
      DT::datatable(rownames = F, options = list(pageLength = 20))  %>%
      DT::formatRound(c('log2 Fold Change','Mean Expression'), digits=2) %>%
      DT::formatRound('Fold Change', digits=6)
    stat_join
  })
  
  # temporal retina heatmap -------
  temporal_retina_heatmap_func <- eventReactive(input$pan_button_temporal_heatmap, {
    cat(file=stderr(), 'Temporal heatmap call\n')
    table <- input$temporal_retina_heatmap_table
    
    if (table == 'Gene 2019') {
      pool <- 'gene_pool_2019'
    } else {
      pool <- 'gene_pool_2022'
    }
    
    if (table == 'Transcript 2019') {
      table <- 'lsTPM_tx'
    } else {
      table <- 'lsTPM_gene'
    }
    gene <- input$temporal_retina_heatmap_ID
    if (1 %in% input$temporal_retina_heatmap_clustering){
      cluster_row <- TRUE
    } else {cluster_row <- FALSE}
    make_heatmap <- function(title,
                             matrix,
                             breaks = c(0,5,10,15),
                             cluster_row,
                             show_row_names = FALSE,
                             show_heatmap_legend = FALSE){
      Heatmap(log2(matrix+1),
              cluster_columns = F,
              column_title = title,
              cluster_rows = cluster_row,
              col = colorRamp2(breaks = breaks, colors = viridis(length(breaks))),
              rect_gp = gpar(col= "white"),
              show_row_names = show_row_names,
              name = 'log2(TPM+1)',
              show_heatmap_legend = show_heatmap_legend,
              clustering_distance_rows = "pearson",
              clustering_distance_columns = "euclidean")
    }
    
    query = paste0('select * from ', table, ' where ID in ("',paste(gene, collapse='","'),'")')
    p <- dbGetQuery(gene_pool_2019, query) %>% left_join(.,core_tight_2019) %>%
      left_join(., gene_pool_2019 %>% tbl('gene_IDs') %>% as_tibble()) %>%
      as_tibble() %>% filter(Tissue %in% c('ESC','Retina'))
    
    ESC <- p %>%
      filter(Tissue == 'ESC') %>%
      mutate(Days = 0, Type = 'ESC') %>%
      group_by(ID, Days) %>%
      summarise(value = mean(value)) %>%
      mutate(Days = as.integer(Days))
    organoid_swaroop_GFP <- p %>%
      filter(Sub_Tissue == 'Retina - 3D Organoid Stem Cell', !grepl('GFP negative', sample_attribute), study_accession != 'SRP159246') %>%
      group_by(ID, Age_Days) %>%
      summarise(value = mean(value)) %>%
      mutate(Days = as.integer(Age_Days), Type = 'GFP+ 3D Organoid') %>%
      select(-Age_Days)
    organoid_swaroop_GFPneg <-  p %>%
      filter(Sub_Tissue == 'Retina - 3D Organoid Stem Cell', grepl('GFP negative', sample_attribute), study_accession != 'SRP159246') %>%
      group_by(ID, Age_Days) %>%
      summarise(value = mean(value)) %>%
      mutate(Days = as.integer(Age_Days), Type = 'Kaewkhaw GFP- 3D Retina')%>%
      select(-Age_Days)
    organoid_johnston <-  p %>%
      filter(study_accession == 'SRP159246') %>%
      group_by(ID, Age_Days) %>%
      summarise(value = mean(value)) %>%
      mutate(Days = as.integer(Age_Days), Type = 'Kaewkhaw GFP+ 3D Retina') %>%
      select(-Age_Days)
    fetal_tissue <- p %>%
      filter(Sub_Tissue == 'Retina - Fetal Tissue') %>%
      group_by(ID, Age_Days) %>%
      summarise(value = mean(value)) %>%
      mutate(Days = as.integer(Age_Days), Type = 'Fetal Tissue') %>%
      select(-Age_Days)
    adult_tissue <- p %>%
      filter(Sub_Tissue == 'Retina - Adult Tissue') %>%
      group_by(ID) %>%
      summarise(value = mean(value), Type = 'Adult Tissue') %>%
      mutate(Days = 1000)
    if (input$temporal_retina_heatmap_viz == 'Split by type'){
      
      # tissue
      tissue <- bind_rows(fetal_tissue, adult_tissue)
      matrix <- tissue %>% select(-Type) %>% spread(ID, value) %>% t()
      colnames(matrix) <- matrix['Days',]
      colnames(matrix)[ncol(matrix)] <- 'Adult'
      matrix <- matrix[-1,]
      one <- make_heatmap('Retina Tissue', matrix, show_heatmap_legend = T, cluster_row = cluster_row)
      
      # swaroop GFP+
      x <- rbind(organoid_swaroop_GFP, ESC)
      y <- x %>% select(-Type) %>% spread(ID, value) %>% t()
      colnames(y) <- y['Days',]
      colnames(y)[1] <- 'ESC'
      y <- y[-1,]
      two <- make_heatmap(title = 'GFP+ 3D\nRetina\n(Kaewkhaw)', y, cluster_row = cluster_row)
      
      # swaroop GFP-
      x <- rbind(organoid_swaroop_GFPneg, ESC)
      y <- x %>% select(-Type) %>% spread(ID, value) %>% t()
      colnames(y) <- y['Days',]
      colnames(y)[1] <- 'ESC'
      y <- y[-1,]
      three <- make_heatmap('GFP- 3D\nRetina\n(Kaewkhaw)', y, show_row_names = T, cluster_row = cluster_row)
      
      # johnston
      x <- rbind(organoid_johnston, ESC)
      y <- x %>% select(-Type) %>% spread(ID, value) %>% t()
      colnames(y) <- y['Days',]
      colnames(y)[1] <- 'ESC'
      y <- y[-1,]
      four <- make_heatmap('3D Retina (Eldred)', y, cluster_row = cluster_row)
      
      one + four + two + three
    } else {
      breaks = c(0,5,10,15)
      tissue <- bind_rows(fetal_tissue %>% mutate(Type = 'Tissue'),
                          adult_tissue %>% mutate(Type = 'Tissue'),
                          ESC %>% mutate(Type = 'ESC'),
                          organoid_johnston %>% mutate(Type = 'Eldred'),
                          organoid_swaroop_GFP %>% mutate(Type = 'Kaewkhaw GFP+'),
                          organoid_swaroop_GFPneg %>% mutate(Type = 'Kaewkhaw GFP-')) %>%
        mutate(Range = case_when(Days == 0 ~ 0,
                                 Days <= 10 ~ 10,
                                 Days <= 20 ~ 20,
                                 Days < 40 ~ 35,
                                 Days < 65 ~ 55,
                                 Days < 85 ~ 75,
                                 Days < 105 ~ 100,
                                 Days < 125 ~ 120,
                                 Days < 145 ~ 140,
                                 Days < 185 ~ 180,
                                 Days < 255 ~ 250,
                                 TRUE ~ 1000)) %>%
        group_by(ID, Range) %>% summarise(value = mean(value))
      matrix <- tissue %>% spread(ID, value) %>% select(-Range) %>% t()
      days <- (tissue %>% spread(ID, value) %>% t())[1,]
      days[1] <- 'ESC'
      days[length(days)] <- 'Adult'
      colnames(matrix) <- days
      Heatmap(log2(matrix+1), cluster_columns = F,
              col = colorRamp2(breaks = breaks, colors = viridis(length(breaks))),
              clustering_distance_rows = "euclidean",
              rect_gp = gpar(col= "white"),
              name = 'log2(TPM + 1)',
              show_heatmap_legend = T)
    }
  })
  temporal_retina_heatmap_height <- eventReactive(input$pan_button_temporal_heatmap, {
    if (input$temporal_retina_heatmap_viz == 'Split by type'){
      max(28*length(input$temporal_retina_heatmap_ID), 150)
    } else {
      max(25*length(input$temporal_retina_heatmap_ID), 75)
    }
  })
  output$temporal_retina_heatmap <- renderPlot({
    temporal_retina_heatmap_func()
  },  height=temporal_retina_heatmap_height)
  
  # modal for custom url ---------
  observeEvent(input$build_pan_url, {
    url = paste0("http://eyeIntegration.nei.nih.gov/?Dataset=",
                 gsub(' ' , '_', paste(input$Database, collapse = ',')),
                 "&ID=", gsub(' ', '_', paste(input$ID, collapse = ',')),
                 "&Tissue=", gsub(' ', '_', paste(input$plot_tissue_gene, collapse = ',')),
                 "&num=", as.character(input$num_gene))
    showModal(modalDialog(
      title = 'Paste this into your browser to jump to these parameters in the future:',
      url,
      easyClose = T,
      footer = "Warning! With the custom URL you will be unable to change the Dataset"
    ))
  })
  
  # single cell info page -----
  # first set (see ui.R)
  # table
  
  scboxPlot_gene_func <- eventReactive(input$scpan_button_gene, {
    cat(file=stderr(), 'scboxPlot Gene call\n')
    db <- scEiaD_pool
    pool <- scEiaD_pool
    maturity <- input$scmaturity
    gene <- input$scGene
    tissue <- input$scplot_tissue_gene
    col_num <- input$scnum_gene
    if (length(db) < 1 || length(gene) < 1 || length(tissue) < 1){
      showModal(modalDialog(title = "Box Plot Warning",
                            "Have you specified at least one gene or tissue?", 
                            easyClose = T,
                            footer = NULL))
    }
    
    query = paste0('select * from scEiaD_CT_table where Gene in ("',paste(gene, collapse='","'),'")')
    p <- dbGetQuery(scEiaD_pool, query) %>% 
      as_tibble()
    
    p$Type <- p %>% select(contains('type')) %>% pull(1)
    
    if (!is.null(maturity)){
      p <- p %>% filter(Stage %in% maturity)
    }
    
    p <- p %>% 
      filter(CellType_predict %in% tissue) %>% 
      mutate(Info = paste('Study: ', 
                          study_accession,
                          sep ='')) %>% 
      ggplot(data=.,aes(x=CellType_predict,y=`Meanlog2(Counts+1)`,colour=CellType_predict)) +
      geom_boxplot(alpha=0.5, outlier.shape = NA) + 
      geom_point_interactive(size=2, position = 'jitter', alpha=0.25, stroke = 3, aes(tooltip=htmlEscape(Info, TRUE), fill = CellType_predict)) + 
      xlab('') + 
      facet_wrap(~Gene + Stage, ncol=col_num) +
      cowplot::theme_cowplot(font_size = 15) + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.2)) +
      ggtitle('Box Plot of Single Cell Gene Expression') +
      ylab("Mean log2(Counts+1)") +
      scale_shape_manual(values=c(0:2,5,6,15:50)) +
      theme(plot.margin=grid::unit(c(0,0,0,0.1), "cm")) +
      CellType_predict_col
    girafe(ggobj = p,width_svg = 12, 
           height_svg= max(10, (6 * (length(gene)/min(col_num,length(gene)))))) %>% 
      girafe_options(., opts_toolbar(position = "top") )
    
  })
  output$scboxPlot_gene <- renderGirafe({
    scboxPlot_gene_func()
  })
  
  # table
  output$SC_gene_means_by_type_table <- DT::renderDataTable(server = TRUE, {
    input$SC_density_pan_button
    isolate({
      req(input$mGene)
      mouse_genes <- c(input$mGene)
      SC_dataset <- (input$SC_dataset %>% strsplit(' '))[[1]][1] %>% tolower()
    })
    table_name <- paste(SC_dataset, 'gene_cell_type_stats', sep='__')
    if (input$single_cell_stat_type == 'Mean Gene Expression (across all cells, split by cell type)') {
      dec_column <- 'Decile_mean'
      dec_name <- 'Decile'
      rank_column <- 'Rank_mean'
    } else {
      dec_column <- 'Percentage Cells'
      dec_name <- 'Percentage'
      rank_column <- 'Rank_cells'
    }
    gene_type_counts <- dbGetQuery(SC_pool,paste0('SELECT * FROM ', table_name, ' WHERE Gene IN ("', 
                                                  paste(mouse_genes, collapse='","'), '") '))
    DT::datatable(gene_type_counts %>% 
                    filter(`Total Count` > 10) %>% 
                    dplyr::select(Gene:`Total Count`, 'Mean Expression Count of Gene in Tissue' = Mean, 'Rank' = !!(rank_column), !!dec_name := !!(dec_column)),
                  extensions = 'Buttons', rownames = FALSE, options = list(
                    dom = 'frtBip', searchable = TRUE, buttons = c('pageLength', 'copy', 'csv'), pageLength=12)) %>% 
      DT::formatRound(c('Mean Expression Count of Gene in Tissue'), digits=2)
    
  })
  # heatmap for single cell ------------
  output$SC_gene_means_by_type_heatmap <- renderPlot({
    input$SC_density_pan_button
    isolate({
      req(input$mGene)
      gene = c(input$mGene)[1]
      SC_dataset <- (input$SC_dataset %>% strsplit(' '))[[1]][1] %>% tolower()
    })
    
    metadata_name <- paste(SC_dataset, 'SC_metadata_long', sep = '__')
    table_name <- paste(SC_dataset, 'gene_cell_type_stats', sep='__')
    if (input$single_cell_stat_type == 'Decile of Mean Gene Expression (across all cells, split by cell type)') {
      dec_column <- 'Decile_mean'
      rank_column <- 'Rank_mean'
      legend_name <- 'Decile Mean\nGene Expression'
    } else if (input$single_cell_stat_type == 'Percentage Cells Expressing Gene (across all cells, split by cell type and age (if available))') {
      dec_column <- 'Percentage Cells'
      rank_column <- 'Rank_cells'
      legend_name <- '% cells\nExpressing Gene'
    }
    # replace missing time points and cell types
    add_missing <- function(mat,
                            replace = 0,
                            expected_cols = c("E11", "E12", "E14", "E16", "E18", "P0", "P2", "P5", "P8", "P14")){
      expected_rows = dbGetQuery(SC_pool, paste('SELECT * from ', metadata_name)) %>%
        filter(`Cell Type` != 'Doublets', `Cell Type` != 'Red Blood Cells') %>%
        pull(`Cell Type`) %>%
        unique() %>%
        sort()
      out <- mat %>% data.frame()
      row.names(out) <- mat %>% pull(1)
      out <- out %>% select(2:ncol(out))  %>% replace(is.na(.), replace)
      missing_col <- setdiff(expected_cols, names(out))
      missing_row <- setdiff(expected_rows,  row.names(out))
      out[,missing_col] <- replace
      out[missing_row,] <- replace
      out <- out[expected_rows, expected_cols, drop=FALSE]
      out
    }
    if (SC_dataset == 'clark') {
      mat <- SC_pool %>%
        tbl(table_name) %>%
        filter(Gene == gene, `Total Count` > 10, `Cell Type` != 'Doublets', `Cell Type` != 'Red Blood Cells') %>%
        select(Age, `Cell Type`, !!(dec_column)) %>%
        as_tibble() %>%
        spread(Age, !!(dec_column))
      vals <- SC_pool %>%
        tbl(table_name) %>%
        filter(Gene == gene, `Total Count` > 10, `Cell Type` != 'Doublets', `Cell Type` != 'Red Blood Cells') %>%
        select(Age, `Cell Type`, !!(rank_column)) %>%
        as_tibble() %>%
        spread(Age, !!(rank_column))
      
      mat_CH <- add_missing(mat) %>% as.matrix()
      vals_CH <- add_missing(vals, 'ND') %>% as.matrix()
    } else {
      # macosko
      mat_CH <- add_missing(SC_pool %>% tbl(table_name) %>% filter(Gene == gene, `Total Count` > 10) %>% select(`Cell Type`, !!(dec_column)) %>% data.frame(), expected_cols = 'Decile_mean') %>% as.matrix()
      colnames(mat_CH)[1] <- 'P14'
      mat_CH <- cbind(mat_CH, NA)
      
      vals_CH <- add_missing(SC_pool %>% tbl(table_name) %>% filter(Gene == gene, `Total Count` > 10) %>% select(`Cell Type`, !!(rank_column)) %>% data.frame(), expected_cols = 'Rank_mean') %>% as.matrix()
      colnames(vals_CH)[1] <- 'P14'
      vals_CH <- cbind(vals_CH, 'ND')
    }
    row.names(mat_CH) <- gsub('r M', 'r\nM', row.names(mat_CH))
    row.names(vals_CH) <- gsub('r M', 'r\nM', row.names(vals_CH))
    # colors_CH <- mat_CH
    # colors_CH <- ifelse(colors_CH < 3, 'gray', 'black')
    if (input$heatmap_overlay == 'None') {
      layer_fun <- NA
    } else{
      layer_fun <- function(j, i, x, y, width, height, fill) {
        # since grid.text can also be vectorized
        v = pindex(vals_CH, i, j)
        grid.shadowtext(sprintf("%s",v), x, y, gp =
                          gpar(fontsize = 12,
                               col = 'white'),
                        bg.r = 0.07)
      }
    }
    
    breaks = seq(0,max(mat_CH), by = (max(mat_CH)/5) )
    
    Heatmap(mat_CH,
            cluster_columns = FALSE,
            #column_title = title,
            cluster_rows = FALSE,
            col = colorRamp2(breaks = breaks, colors = viridis(length(breaks))),
            rect_gp = gpar(col= "white"),
            show_row_names = TRUE,
            name = legend_name,
            show_heatmap_legend = TRUE,
            clustering_distance_rows = "pearson",
            clustering_distance_columns = "euclidean",
            layer_fun = layer_fun)
    
  })
  
  # SC tSNE t-SNE (mouse retina for now) --------
  output$single_cell_tsne_plot <- renderGirafe({
    input$SC_tsne_pan_button
    isolate({
      req(input$mGene_tsne)
      req(input$age_tsne)
      req(input$SC_dataset_tsne)
      req(input$min_single_cell_gene_count)
      mouse_gene <- input$mGene_tsne
      input_age <- c(input$age_tsne)
      min_single_cell_gene_count <- input$min_single_cell_gene_count
      SC_dataset <- (input$SC_dataset_tsne %>% strsplit(' '))[[1]][1] %>% tolower()
      table_name_tsne <- paste(SC_dataset, 'tsne_coords', sep='__')
      table_name_gc <- paste(SC_dataset, 'SC_gene_counts', sep='__')
    })
    tsne_coords <- dbGetQuery(SC_pool, paste('SELECT * FROM ', table_name_tsne)) %>%
      filter(!grepl('Red|Doub', `Cell Type`)) # remove Red Blood Cells and Doublets
    # identify cells with gene expression count above min_single_cell_gene_count
    samples_gene_up <- dbGetQuery(SC_pool, paste0('SELECT `Cell ID` FROM ', table_name_gc, ' WHERE Gene is "', mouse_gene, '" AND `Gene Count` > ',
                                                  min_single_cell_gene_count))
    
    if (nrow(samples_gene_up) == 0) {
      interactive <- NULL
    } else {
      if (SC_dataset == 'macosko') {
        interactive <- geom_point_interactive(data = tsne_coords %>%
                                                filter(`Cell ID` %in% samples_gene_up$`Cell ID`),
                                              aes(tooltip=htmlEscape(paste0('Cell Type: ', `Cell Type`, '\nCell ID: ', `Cell ID`), TRUE),
                                                  fill = `Cell Type`),
                                              alpha = 0.6,
                                              size = 1.4,
                                              shape = 21,
                                              color = 'white')
      } else{
        interactive <- geom_point_interactive(data = tsne_coords %>%
                                                filter(age %in% input_age) %>%
                                                filter(`Cell ID` %in% samples_gene_up$`Cell ID`),
                                              aes(tooltip=htmlEscape(paste0('Cell Type: ', `Cell Type`, '\nCell ID: ', `Cell ID`), TRUE),
                                                  fill = `Cell Type`),
                                              alpha = 0.6,
                                              size = 1.4,
                                              shape = 21,
                                              color = 'white')
      }
    }
    if (SC_dataset == 'macosko'){
      p <- ggplot(tsne_coords, aes(x=tSNE_1,y=tSNE_2, fill = `Cell Type`))  +
        geom_point(alpha=0.1, size=1, shape = 21, aes(fill = `Cell Type`, colour = `Cell Type`), stroke = 0) +
        xlab('tSNE 1') +
        ylab('tSNE 2') +
        guides(colour = FALSE) +
        guides(fill = guide_legend(override.aes = list(alpha = 1))) +
        cowplot::theme_cowplot(font_size=8) +
        guides(fill=guide_legend(nrow = 4,byrow=TRUE)) + interactive
      girafe(code = print(p)) %>% girafe_options(., opts_toolbar(position = "top") )
      
    } else {
      p <- ggplot(tsne_coords %>%
                    filter(age %in% input_age),
                  aes(x=umap_coord1,y=umap_coord2))  +
        geom_point(alpha=0.1, size=1, shape = 21, aes(fill = `Cell Type`, colour = `Cell Type`), stroke = 0) +
        facet_wrap(~age) +
        xlab('UMAP 1') +
        ylab('UMAP 2') +
        guides(colour = FALSE) +
        guides(fill = guide_legend(override.aes = list(alpha = 1))) +
        cowplot::theme_cowplot(font_size=8) +
        guides(fill=guide_legend(nrow = 4,byrow=TRUE)) + interactive
      girafe(code = print(p)) %>% girafe_options(., opts_toolbar(position = "top") )
    }
  })
  
  # Bulk t-SNE -------------
  bulk_tsne <- eventReactive(input$bulk_tsne_button, {
    isolate({
      req(input$bulk_tsne_dataset)
      req(input$perplexity)
      dataSET = input$bulk_tsne_dataset
      perp = input$perplexity
    })
    if (dataSET == '2017'){
      load('./www/all_tsne_plot_prepped__2017_02.Rdata') # tsne 5->50 perplexity for bulk RNA. script to create called 'dbscan_interactive_page_calculate.R'
      perplexity_level <- perp
      tsne_plot<- all_tsne_plot_prepped %>% filter(perplexity==perplexity_level)
      
      p <- ggplot(tsne_plot) +
        ggtitle('Pan tissue t-SNE') +
        geom_point_interactive(size=20, alpha=0.2, aes(x=X1,y=X2,colour=Cluster, tooltip=htmlEscape(Label, TRUE))) +
        geom_point(data=tsne_plot %>% dplyr::select(X1,X2), size=5, alpha=0.2, colour='black', aes(x=X1,y=X2)) +
        xlab('t-SNE 1') + ylab('t-SNE 2') +
        theme_minimal() +
        guides(colour = guide_legend(ncol = 3,
                                     override.aes = list(size=10, alpha = 1)))
    } else if (dataSET == '2019'){
      all_tsne_plot_prepped <- gene_pool_2019 %>%
        tbl('tSNE_bulk_RNA') %>%
        filter(perplexity == perp)
      p <- ggplot(all_tsne_plot_prepped %>% as_tibble() %>%
                    mutate(Label = case_when(grepl('^0', Label) ~ '', TRUE ~ Label),
                           Info = paste('SRA: ',
                                        sample_accession,
                                        '\nStudy: ',
                                        study_title, '\n',
                                        gsub('\\|\\|', '\n',
                                             sample_attribute),
                                        sep ='')),
                  aes(x=X1,y=X2)) +
        scale_shape_manual(values=c(0:2,5,6,15:20)) +
        ggtitle('Pan tissue t-SNE') +
        geom_point(size=20, alpha=0.1, aes(colour=Tissue)) +
        geom_point_interactive(size=5, alpha=0.6, aes(shape=Origin, tooltip = htmlEscape(Info, TRUE))) +
        geom_label_repel(aes(label=Label), alpha=0.8, size=4, box.padding = unit(0.3, "lines")) +
        theme_minimal() +
        xlab('t-SNE 1') +
        ylab('t-SNE 2') +
        tissue_col +
        guides(colour = guide_legend(ncol = 3,
                                     override.aes = list(size=10, alpha = 1)))
    }
    ggiraph(code = print(p),selection_type = 'none', zoom_max=3, width_svg = 14, height_svg = 10)
    
  })
  output$tsne <- renderGirafe({
    bulk_tsne()
  })
  
  
  # diff expression ------------------------------------
  output$diff.exp <- DT::renderDataTable(server = TRUE, {
    req(input$de_comparison)
    if (input$diff_database == 'Gene 2017') {
      de_data <- gene_pool_2017 %>%
        tbl('limma_DE_gene') %>%
        filter(Comparison == local(input$de_comparison)) %>%
        as_tibble()
    } else if (input$diff_database == 'Transcript 2019') {
      de_data <- gene_pool_2019 %>%
        tbl('limma_DE_tx') %>%
        filter(Comparison == local(input$de_comparison)) %>%
        as_tibble() %>%
        # left join to filter on gene or tx type (e.g. protein coding / miRNA / etc)
        left_join(., gene_pool_2019 %>% tbl('tx_IDs') %>% as_tibble())  %>%
        filter(transcript_type %in% input$gene_tx_type) %>%
        rename(Class = transcript_type)
    } else if (input$diff_database == 'Gene 2019'){
      de_data <- gene_pool_2019 %>%
        tbl('limma_DE_gene') %>%
        filter(Comparison == local(input$de_comparison)) %>%
        as_tibble() %>%
        # left join to filter on gene or tx type (e.g. protein coding / miRNA / etc)
        left_join(., gene_pool_2019 %>% tbl('gene_IDs') %>% as_tibble())  %>%
        filter(gene_type %in% input$gene_tx_type) %>%
        rename(Class = gene_type)
    } else if (input$diff_database == 'DNTx v01'){
      de_data <- DNTx_pool_2019 %>%
        tbl('limma_DE_tx') %>%
        filter(Comparison == local(input$de_comparison)) %>%
        as_tibble() %>%
        # left join to filter on gene or tx type (e.g. protein coding / miRNA / etc)
        left_join(., DNTx_pool_2019 %>% tbl('tx_IDs') %>% as_tibble())  %>%
        filter(transcript_type %in% input$gene_tx_type) %>%
        rename(Class = transcript_type)
    }
    de_data %>%
      dplyr::select(-Comparison) %>%
      mutate(`P.Value` = format(`P.Value`, digits = 3),
             `adj.P.Val` = format(`adj.P.Val`, digits = 3)) %>%
      DT::datatable(extensions = 'Buttons', options = list(
        dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatRound(c('logFC','AveExpr','t','B'), digits=2)
  })
  
  # word clouds ---------
  output$word_cloud_image_up <- renderUI({
    req(input$de_comparison)
    
    if (grepl('2017', input$diff_database)){
      cloud_path <- './2017/word_cloud_png/'
    } else {
      cloud_path <- './2019/word_cloud_png/'
    }
    src = paste0(cloud_path, input$de_comparison, '__Up.png')
    print(src)
    tags$img(src = src)})
  
  output$word_cloud_image_down <- renderUI({
    req(input$de_comparison)
    if (grepl('2017', input$diff_database)){
      cloud_path <- './2017/word_cloud_png/'
    } else {
      cloud_path <- './2019/word_cloud_png/'
    }
    src = paste0(cloud_path, input$de_comparison, '__Down.png')
    tags$img(src = src)})
  
  output$comparison_up1 <- renderText({
    req(input$de_comparison)
    if (grepl('2017', input$diff_database)){
      load('./www/2017/de_comparison_name_list.Rdata')
    } else {
      de_comparison_contrast_names <- gene_pool_2019 %>% tbl('limma_DE_tests') %>% pull(Comparison)
      names(de_comparison_contrast_names) <- gene_pool_2019 %>% tbl('limma_DE_tests') %>% pull(Name)
    }
    gsub('vs', '>', names(de_comparison_contrast_names[de_comparison_contrast_names==input$de_comparison]))
    print(gsub('vs', '>', names(de_comparison_contrast_names[de_comparison_contrast_names==input$de_comparison])))
  })
  output$comparison_down1 <- renderText({
    req(input$de_comparison)
    if (grepl('2017', input$diff_database)){
      load('./www/2017/de_comparison_name_list.Rdata')
    } else {
      de_comparison_contrast_names <- gene_pool_2019 %>% tbl('limma_DE_tests') %>% pull(Comparison)
      names(de_comparison_contrast_names) <- gene_pool_2019 %>% tbl('limma_DE_tests') %>% pull(Name)
    }
    gsub('vs', '<', names(de_comparison_contrast_names[de_comparison_contrast_names==input$de_comparison]))
  })
  output$comparison_up2 <- renderText({
    req(input$de_comparison)
    if (grepl('2017', input$diff_database)){
      load('./www/2017/de_comparison_name_list.Rdata')
    } else {
      de_comparison_contrast_names <- gene_pool_2019 %>% tbl('limma_DE_tests') %>% pull(Comparison)
      names(de_comparison_contrast_names) <- gene_pool_2019 %>% tbl('limma_DE_tests') %>% pull(Name)
    }
    gsub('vs', '>', names(de_comparison_contrast_names[de_comparison_contrast_names==input$de_comparison]))
  })
  output$comparison_down2 <- renderText({
    req(input$de_comparison)
    if (grepl('2017', input$diff_database)){
      load('./www/2017/de_comparison_name_list.Rdata')
    } else {
      de_comparison_contrast_names <- gene_pool_2019 %>% tbl('limma_DE_tests') %>% pull(Comparison)
      names(de_comparison_contrast_names) <- gene_pool_2019 %>% tbl('limma_DE_tests') %>% pull(Name)
    }
    gsub('vs', '<', names(de_comparison_contrast_names[de_comparison_contrast_names==input$de_comparison]))
  })
  # diff expression GO tables ------
  output$go.table.up <- DT::renderDataTable(server = TRUE, {
    req(input$de_comparison)
    if (input$diff_database == 'Gene 2017') {
      go <- gene_pool_2017 %>%
        tbl('all_vs_all_go') %>%
        filter(Set == local(input$de_comparison), Test=='Up') %>%
        dplyr::select(`GO ID`:Term,Ontology)
    } else {
      go <- gene_pool_2019 %>%
        tbl('all_vs_all_go') %>%
        filter(Set == local(input$de_comparison), Test=='Up') %>%
        dplyr::select(ONTOLOGY:`p.adjust`)
    }
    go %>%
      as_tibble()} %>%
      DT::datatable(extensions = 'Buttons', options = list(
        dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))))
  output$go.table.down <- DT::renderDataTable(server = TRUE, {
    req(input$de_comparison)
    if (input$diff_database == 'Gene 2017') {
      go <- gene_pool_2017 %>%
        tbl('all_vs_all_go') %>%
        filter(Set == local(input$de_comparison), Test=='Down') %>%
        dplyr::select(`GO ID`:Term,Ontology)
    } else {
      go <- gene_pool_2019 %>%
        tbl('all_vs_all_go') %>%
        filter(Set == local(input$de_comparison), Test=='Down') %>%
        dplyr::select(ONTOLOGY:`p.adjust`)
    }
    go %>%
      as_tibble()} %>%
      DT::datatable(extensions = 'Buttons', options = list(
        dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))))
  
  # Find a Friend ------------------
  output$FaF_euc_dist <- renderDataTable(server = TRUE, {
    DT::datatable(gene_pool_2019 %>% tbl('Euc_Dist_Top_100') %>% filter(ID1 == input$FaF_ID) %>% as_tibble() %>%
                    mutate(Distance = as.integer(Distance)),
                  rownames = F, extensions = 'Buttons', options = list(
                    dom = 'rtBip', buttons = c('pageLength','copy', 'csv')))
  })
  
  # retina network -----------------------
  output$retina_network_mod <- renderVisNetwork({
    networks.retina.list[[paste0(input$retina_mod_color_vis_mod, "_k", input$retina_kNN)]]
  })
  output$retina_network_gene <- renderVisNetwork({
    req(input$retina_gene)
    new_mod_color = gene_pool_2017 %>% tbl('retina_gene_name_colors') %>% filter(id == local(input$retina_gene)) %>% pull(Module_Color)
    networks.retina.list[[paste0(new_mod_color, "_k", input$retina_kNN)]] %>%
      visOptions(nodesIdSelection = list(enabled = T, selected = input$retina_gene),
                 highlightNearest = list(enabled = T, degree = 0, hover = F, algorithm = "all"))
  })
  output$retina_table_mod <- renderDataTable(server = TRUE, {
    DT::datatable(gene_pool_2017 %>% tbl('edges_retina') %>% filter(Color == local(input$retina_mod_color_vis_mod)) %>% as_tibble(), rownames = F, extensions = 'Buttons', options = list(
      dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('length'), digits = 5)
  })
  output$retina_table_gene <- renderDataTable(server = TRUE, {
    req(input$retina_gene)
    new_mod_color = gene_pool_2017 %>% tbl('retina_gene_name_colors') %>% filter(id == local(input$retina_gene)) %>% pull(Module_Color)
    #edges.retina.mod = edges.retina.list[[new_mod_color]]
    DT::datatable(gene_pool_2017 %>%
                    tbl('edges_retina') %>%
                    filter(Color == new_mod_color) %>%
                    filter(Gene1 == local(input$retina_gene) | Gene2 == local(input$retina_gene)) %>%
                    as_tibble(),
                  rownames = F, extensions = 'Buttons', options = list(
                    dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('length'), digits = 5)
  })
  output$retina_connect_table_mod <- renderDataTable(server = TRUE, {
    DT::datatable(gene_pool_2017 %>% tbl('retina_mod_connect') %>%
                    filter(`Module.Color` == local(input$retina_mod_color_vis_mod)) %>%
                    as_tibble(), extensions = 'Buttons', options = list(
                      dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('kWithin'), digits = 5)
  })
  output$retina_connect_table_gene <- renderDataTable(server = TRUE, {
    req(input$retina_gene)
    new_mod_color = gene_pool_2017 %>% tbl('retina_gene_name_colors') %>% filter(id == local(input$retina_gene)) %>% pull(Module_Color)
    DT::datatable(gene_pool_2017 %>% tbl('retina_mod_connect') %>%
                    filter(`Module.Color` == new_mod_color) %>%
                    as_tibble(), extensions = 'Buttons', options = list(
                      dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('kWithin'), digits = 5)
  })
  output$retina_GO_table_mod <- renderDataTable(server = TRUE, {
    DT::datatable(gene_pool_2017 %>% tbl('retina_network_GO') %>% filter(Color == local(input$retina_mod_color_vis_mod)) %>% as_tibble(),
                  extensions = 'Buttons', options = list(
                    dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
  })
  output$retina_GO_table_gene <- renderDataTable(server = TRUE, {
    req(input$retina_gene)
    new_mod_color = gene_pool_2017 %>% tbl('retina_gene_name_colors') %>% filter(id == local(input$retina_gene)) %>% pull(Module_Color)
    DT::datatable(gene_pool_2017 %>% tbl('retina_network_GO') %>% filter(Color == new_mod_color) %>% as_tibble(),
                  extensions = 'Buttons', options = list(
                    dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
  })
  output$retina_full_edge_table = DT::renderDataTable(server = TRUE, {
    
    req(input$retina_gene_edges)
    
    if(input$retina_edge_show == 'Only Outside Module'){
      table_display <- dbGetQuery(gene_pool_2017, paste0('select * from retina_network_edges WHERE Gene1 == "', input$retina_gene_edges,'" OR Gene2 == "', input$retina_gene_edges, '"')) %>%
        filter(same.module == 0) %>% mutate(Decile = ntile(length, 10)) %>% dplyr::select(-same.module) %>% arrange(desc(length))
    }else if(input$retina_edge_show == 'Only Within Module'){
      table_display <- dbGetQuery(gene_pool_2017, paste0('select * from retina_network_edges WHERE Gene1 == "', input$retina_gene_edges,'" OR Gene2 == "', input$retina_gene_edges, '"')) %>%
        filter(same.module == 1) %>% mutate(Decile = ntile(length, 10)) %>% dplyr::select(-same.module) %>% arrange(desc(length))
    }else{
      table_display <- dbGetQuery(gene_pool_2017, paste0('select * from retina_network_edges WHERE Gene1 == "', input$retina_gene_edges,'" OR Gene2 == "', input$retina_gene_edges, '"')) %>%
        mutate(Decile = ntile(length, 10)) %>% dplyr::select(-same.module) %>% arrange(desc(length))
    }
    table_display %>%   DT::datatable(extensions = 'Buttons', options = list(
      pageLength = 20,   lengthMenu = c(20, 100, 1000, 10000), dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('length'), digits = 5)
  })
  
  # rpe network --------------
  output$rpe_network_mod <- renderVisNetwork({
    networks.rpe.list[[paste0(input$rpe_mod_color_vis_mod, "_k", input$rpe_kNN)]]
  })
  output$rpe_network_gene <- renderVisNetwork({
    req(input$rpe_gene)
    new_mod_color = gene_pool_2017 %>% tbl('rpe_gene_name_colors') %>% filter(id == local(input$rpe_gene)) %>% pull(Module_Color)
    networks.rpe.list[[paste0(new_mod_color, "_k", input$rpe_kNN)]] %>%
      visOptions(nodesIdSelection = list(enabled = T, selected = input$rpe_gene),
                 highlightNearest = list(enabled = T, degree = 0, hover = F, algorithm = "all"))
  })
  output$rpe_table_mod <- renderDataTable(server = TRUE, {
    DT::datatable(gene_pool_2017 %>% tbl('edges_rpe') %>% filter(Color == local(input$rpe_mod_color_vis_mod)) %>% as_tibble(), rownames = F, extensions = 'Buttons', options = list(
      dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('length'), digits = 5)
  })
  output$rpe_table_gene <- renderDataTable(server = TRUE, {
    req(input$rpe_gene)
    new_mod_color = gene_pool_2017 %>% tbl('rpe_gene_name_colors') %>% filter(id == local(input$rpe_gene)) %>% pull(Module_Color)
    #edges.rpe.mod = edges.rpe.list[[new_mod_color]]
    DT::datatable(gene_pool_2017 %>%
                    tbl('edges_rpe') %>%
                    filter(Color == new_mod_color) %>%
                    filter(Gene1 == local(input$rpe_gene) | Gene2 == local(input$rpe_gene)) %>%
                    as_tibble(),
                  rownames = F, extensions = 'Buttons', options = list(
                    dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('length'), digits = 5)
  })
  output$rpe_connect_table_mod <- renderDataTable(server = TRUE, {
    DT::datatable(gene_pool_2017 %>% tbl('rpe_mod_connect') %>%
                    filter(`Module.Color` == local(input$rpe_mod_color_vis_mod)) %>%
                    as_tibble(), extensions = 'Buttons', options = list(
                      dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('kWithin'), digits = 5)
  })
  output$rpe_connect_table_gene <- renderDataTable(server = TRUE, {
    req(input$rpe_gene)
    new_mod_color = gene_pool_2017 %>% tbl('rpe_gene_name_colors') %>% filter(id == local(input$rpe_gene)) %>% pull(Module_Color)
    DT::datatable(gene_pool_2017 %>% tbl('rpe_mod_connect') %>%
                    filter(`Module.Color` == new_mod_color) %>%
                    as_tibble(), extensions = 'Buttons', options = list(
                      dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('kWithin'), digits = 5)
  })
  output$rpe_GO_table_mod <- renderDataTable(server = TRUE, {
    DT::datatable(gene_pool_2017 %>% tbl('rpe_network_GO') %>% filter(Color == local(input$rpe_mod_color_vis_mod)) %>% as_tibble(),
                  extensions = 'Buttons', options = list(
                    dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
  })
  output$rpe_GO_table_gene <- renderDataTable(server = TRUE, {
    req(input$rpe_gene)
    new_mod_color = gene_pool_2017 %>% tbl('rpe_gene_name_colors') %>% filter(id == local(input$rpe_gene)) %>% pull(Module_Color)
    DT::datatable(gene_pool_2017 %>% tbl('rpe_network_GO') %>% filter(Color == new_mod_color) %>% as_tibble(),
                  extensions = 'Buttons', options = list(
                    dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
  })
  output$rpe_full_edge_table = DT::renderDataTable(server = TRUE, {
    req(input$rpe_gene_edges)
    if(input$rpe_edge_show == 'Only Outside Module'){
      table_display <- dbGetQuery(gene_pool_2017, paste0('select * from rpe_network_edges WHERE Gene1 == "', input$rpe_gene_edges,'" OR Gene2 == "', input$rpe_gene_edges, '"')) %>%
        filter(same.module == 0) %>% mutate(Decile = ntile(length, 10)) %>% dplyr::select(-same.module) %>% arrange(desc(length))
    }else if(input$rpe_edge_show == 'Only Within Module'){
      table_display <- dbGetQuery(gene_pool_2017, paste0('select * from rpe_network_edges WHERE Gene1 == "', input$rpe_gene_edges,'" OR Gene2 == "', input$rpe_gene_edges, '"')) %>%
        filter(same.module == 1) %>% mutate(Decile = ntile(length, 10)) %>% dplyr::select(-same.module) %>% arrange(desc(length))
    }else{
      table_display <- dbGetQuery(gene_pool_2017, paste0('select * from rpe_network_edges WHERE Gene1 == "', input$rpe_gene_edges,'" OR Gene2 == "', input$rpe_gene_edges, '"')) %>%
        mutate(Decile = ntile(length, 10)) %>% dplyr::select(-same.module) %>% arrange(desc(length))
    }
    table_display %>%   DT::datatable(extensions = 'Buttons', options = list(
      pageLength = 20,   lengthMenu = c(20, 100, 1000, 10000), dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('length'), digits = 5)
    
  })
  
  # bulk RNA data table --------
  output$table =
    DT::renderDataTable(server = TRUE, {
      db <- input$table_db
      core_cols <- c('ID','value')
      if (db == 'Gene 2017'){
        table <- dbGetQuery(gene_pool_2017,  paste0('select * from lsTPM_gene where ID in ("',paste(input$table_gene, collapse='","'),'")') ) %>%
          left_join(.,core_tight_2017)
      } else if (db == 'Gene 2019') {
        table <- dbGetQuery(gene_pool_2019,  paste0('select * from lsTPM_gene where ID in ("',paste(input$table_gene, collapse='","'),'")') ) %>%
          left_join(.,core_tight_2019)
      } else if (db == 'Transcript 2017'){
        table <- dbGetQuery(gene_pool_2017,  paste0('select * from lsTPM_tx where ID in ("',paste(input$table_gene, collapse='","'),'")') ) %>%
          left_join(.,core_tight_2017)
      } else if (db == 'Transcript 2019'){
        table <- dbGetQuery(gene_pool_2019,  paste0('select * from lsTPM_tx where ID in ("',paste(input$table_gene, collapse='","'),'")') ) %>%
          left_join(.,core_tight_2019)
      } else if (db == 'DNTx v01'){
        table <- dbGetQuery(DNTx_pool_2019,  paste0('select * from lsTPM_tx where ID in ("',paste(input$table_gene, collapse='","'),'")') ) %>%
          left_join(.,core_tight_2019)
      }
      table %>% filter(Tissue == local(input$table_tissue)) %>%
        dplyr::select(one_of(c(core_cols, input$table_columns))) %>%
        DT::datatable(extensions = 'Buttons', rownames = F,
                      options = list(
                        pageLength = 20,  lengthMenu = c(5, 10, 20, 100, 1000, 5000),
                        dom = 'frtBip',
                        buttons = c('pageLength','copy', 'csv'))) %>%
        DT::formatRound(c('value'), digits=2)
    })
  
  # SC RNA data table ------------
  output$SCtable_1 =
    DT::renderDataTable(server = TRUE, {
      SC_dataset <- (input$SC_datatable_dataset %>% strsplit(' '))[[1]][1] %>% tolower()
      table_name <- paste(SC_dataset, 'gene_cell_type_stats', sep='__')
      if (SC_dataset == 'Macosko'){
        table <- dbGetQuery(SC_pool, paste0('select Gene, "Cell Count",  "Total Count", "Mean", "Rank_mean", "Decile_mean" from ', table_name, ' WHERE "Cell Type"="', input$sc_datatable_tissue,'"')) %>%
          dplyr::select(Gene:`Total Count`, 'Mean Expression Count of Gene in Tissue' = Mean, 'Rank' = Rank_mean, 'Decile' = Decile_mean)
      } else {
        table <- dbGetQuery(SC_pool, paste0('select Gene, Age, "Cell Count",  "Total Count", "Mean", "Rank_mean", "Decile_mean" from ', table_name, ' WHERE "Cell Type"="', input$sc_datatable_tissue,'"')) %>%
          dplyr::select(Gene:`Total Count`, 'Mean Expression Count of Gene in Tissue' = Mean, 'Rank' = Rank_mean, 'Decile' = Decile_mean)
      }
      
      table %>% DT::datatable(extensions = 'Buttons',
                              rownames = F,
                              options = list(
                                pageLength = 20, dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
        DT::formatRound(c('Mean Expression Count of Gene in Tissue'), digits=2)
    })
  output$SCtable_2 =
    DT::renderDataTable(server = TRUE, {
      SC_dataset <- (input$SC_datatable_dataset %>% strsplit(' '))[[1]][1] %>% tolower()
      table_name <- paste(SC_dataset, 'gene_cell_type_stats', sep='__')
      if (SC_dataset == 'Macosko'){
        table <- dbGetQuery(SC_pool, paste0('select Gene, "Cell Count",  "Total Count", "Percentage Cells", "Rank_cells", "Decile_cells" from ', table_name, ' WHERE "Cell Type"="', input$sc_datatable_tissue,'"')) %>%
          dplyr::select(Gene:`Total Count`, 'Percentage of Cells Expressing Gene' = `Percentage Cells`, 'Rank' = Rank_cells, 'Decile' = Decile_cells )
      } else {
        table <- dbGetQuery(SC_pool, paste0('select Gene, Age, "Cell Count",  "Total Count", "Percentage Cells", "Rank_cells", "Decile_cells" from ', table_name, ' WHERE "Cell Type"="', input$sc_datatable_tissue,'"')) %>%
          dplyr::select(Gene:`Total Count`, 'Percentage of Cells Expressing Gene' = `Percentage Cells`, 'Rank' = Rank_cells, 'Decile' = Decile_cells )
      }
      table %>% DT::datatable(extensions = 'Buttons',
                              rownames = F,
                              options = list(
                                pageLength = 20, dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
        DT::formatRound(c('Percentage of Cells Expressing Gene'), digits=2)
    })
  output$SCtable_3 =
    DT::renderDataTable(server = TRUE, {
      SC_dataset <- (input$SC_datatable_dataset %>% strsplit(' '))[[1]][1] %>% tolower()
      table_name <- paste(SC_dataset, 'gene_cell_type_stats', sep='__')
      if (SC_dataset == 'Macosko'){
        table <- dbGetQuery(SC_pool, paste0('select Gene, "Cell Count",  "Total Count", "Percentage Cell Types", "Rank_cell_types", "Decile_cell_types" from ', table_name, ' WHERE "Cell Type"="', input$sc_datatable_tissue,'"')) %>%
          dplyr::select(Gene:`Total Count`, `Percentage Cell Types`, 'Rank' = Rank_cell_types, 'Decile' = Decile_cell_types )
      } else {
        table <- dbGetQuery(SC_pool, paste0('select Gene, Age, "Cell Count",  "Total Count", "Percentage Cell Types", "Rank_cell_types", "Decile_cell_types" from ', table_name, ' WHERE "Cell Type"="', input$sc_datatable_tissue,'"')) %>%
          dplyr::select(Gene:`Total Count`, `Percentage Cell Types`, 'Rank' = Rank_cell_types, 'Decile' = Decile_cell_types )
      }
      
      table %>% DT::datatable(extensions = 'Buttons',
                              rownames = F,
                              options = list(
                                pageLength = 20, dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
        DT::formatRound(c('Percentage Cell Types'), digits=2)
    })
  output$basic_stats <- renderTable({basic_stats}, striped = FALSE, rownames = F, align = 'rl')
})