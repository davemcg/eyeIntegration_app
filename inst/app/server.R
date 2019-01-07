# server.R
time <- Sys.time()
cat(file = stderr(), 'Server Go!\n')
#options(shiny.trace=TRUE)
options(shiny.sanitize.errors = FALSE)

# load stuff for server
#library(plotly)
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

# pools for sqlite DBs ------------
gene_pool_2017 <- dbPool(drv = SQLite(), dbname = "./www/2017/eyeIntegration_human_2017_01.sqlite", idleTimeout = 3600000)
gene_pool_2019 <- dbPool(drv = SQLite(), dbname = "./www/2019/eyeIntegration_human_expression_2019_01.sqlite", idleTimeout = 3600000)
SC_pool <- dbPool(drv = SQLite(), dbname = "./www/single_cell_retina_info_04.sqlite", idleTimeout = 3600000)

source('./www/theme_Publication.R')
gene_names_2017 <- gene_pool_2017 %>% tbl('gene_IDs') %>% pull(ID)
gene_names_2019 <- gene_pool_2019 %>% tbl('gene_IDs') %>% pull(ID)
geneTX_names_2017 <- gene_pool_2017 %>% tbl('tx_IDs') %>% pull(ID)
geneTX_names_2019 <- gene_pool_2019 %>% tbl('tx_IDs') %>% pull(ID)
core_tight_2017 <- gene_pool_2017 %>% tbl('metadata') %>% as.tibble()
core_tight_2019 <- gene_pool_2019 %>% tbl('metadata') %>% as.tibble()

load('./www/2017/retina_module_network_lists.Rdata') # NOPE THESE ARE PRECOMPUTED htmlwidgets 
load('./www/2017/rpe_module_network_lists.Rdata') # NOPE THESE ARE PRECOMPUTED htmlwidgets 
#load('./www/go_heatmap.Rdata')
load('./www/basic_stats.Rdata')
cat(file=stderr(), 'Data loaded in ')
cat(file=stderr(), Sys.time() - time)
cat(file=stderr(), ' seconds.\n')

onStop(function() {
  poolClose(gene_pool_2017)
})
onStop(function() {
  poolClose(gene_pool_2019)
})
onStop(function() {
  poolClose(SC_pool)
})

core_tight_2017$sample_accession<-gsub('E-MTAB-','E.MTAB.',core_tight_2017$sample_accession)
core_tight_2017$Sub_Tissue <- gsub('_',' - ',core_tight_2017$Sub_Tissue)
core_tight_2019$sample_accession<-gsub('E-MTAB-','E.MTAB.',core_tight_2019$sample_accession)
core_tight_2019$Sub_Tissue <- gsub('_',' - ',core_tight_2019$Sub_Tissue)

# site begins! ---------
shinyServer(function(input, output, session) {
  # Observe: update fields for pan tissue plots --------
  # also handles custom url
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query[['Dataset']])) {
      updateTextInput(session, "Database", value = gsub('_', ' ', as.character(query['Dataset'])))
    }
    db = input$Database # c("Gene 2017", "Gene 2019", "Transcript 2017", "Transcript 2019")
    if (db == 'Gene 2017'){ID_names = gene_names_2017 %>% sort()
    } else if (db == 'Gene 2019'){ID_names = gene_names_2019 %>% sort()
    } else if (db == 'Transcript 2017'){ID_names = geneTX_names_2017 %>% sort()
    } else {ID_names = geneTX_names_2019 %>% sort()}
    # gene / tx lists
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
      if (grepl('2017', db)){tissues <- unique(sort(core_tight_2017$Sub_Tissue))
      } else {tissues <- unique(sort(core_tight_2019$Sub_Tissue))}
      updateSelectizeInput(session, 'plot_tissue_gene',
                           choices= tissues,
                           options = list(placeholder = 'Type to search'),
                           server = TRUE)
    }
    if (!is.null(query[['Tissue']])){
      select_tissue <- gsub(pattern = '_',replacement = ' ', x = as.character(query['Tissue']))
      select_tissue <- strsplit(select_tissue, split = ',')[[1]]
      #select_tissue <- c(select_tissue, input$plot_tissue_gene)
      if (grepl('2017', db)){tissues <- unique(sort(core_tight_2017$Sub_Tissue))
      } else {tissues <- unique(sort(core_tight_2019$Sub_Tissue))}
      updateSelectizeInput(session, 'plot_tissue_gene',
                           choices= tissues,
                           selected = select_tissue,
                           server = TRUE)
    }
  })
  # differential expression
  # update database (2017 and 2019) based on user input
  observe({
    db = input$diff_database # c("Gene 2017", "Gene 2019", "Transcript 2019")
    if (grepl('2017', db)){
      pool <- 'gene_pool_2017'
      ids <- 'gene_IDs'
    } else if (db == 'Gene 2019') {
      pool <- 'gene_pool_2019'
      ids <- 'gene_IDs'
    } else {
      pool <- 'gene_pool_2019'
      ids <- 'tx_ids'
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
  
  # pan - tissue boxplot -------
  boxPlot_gene_func <- eventReactive(input$pan_button_gene, {
    cat(file=stderr(), 'boxPlot Gene call\n')
    db <- input$Database
    gene <- input$ID
    tissue <- input$plot_tissue_gene
    col_num <- input$num_gene
    if (length(db) < 1 | length(gene) < 1 | length(tissue) < 1){
      showModal(modalDialog(title = "Box Plot Warning",
                            "Have you specified at least one gene or tissue?", 
                            easyClose = T,
                            footer = NULL))
    }
    if (db == 'Gene 2017') {
      query = paste0('select * from lsTPM_gene where ID in ("',paste(gene, collapse='","'),'")')
      p <- dbGetQuery(gene_pool_2017, query) %>% 
        left_join(.,core_tight_2017) %>% 
        left_join(., gene_pool_2017 %>% tbl('gene_IDs') %>% as.tibble()) %>% 
        as.tibble()
    } else if (db == 'Transcript 2017'){
      query = paste0('select * from lsTPM_TX where ID in ("',paste(gene, collapse='","'),'")')
      p <- dbGetQuery(gene_pool_2017, query) %>% left_join(.,core_tight_2017) %>% 
        left_join(., gene_pool_2017 %>% tbl('tx_IDs') %>% as.tibble()) %>% 
        as.tibble()
    } else if (db == 'Gene 2019'){
      query = paste0('select * from lsTPM_gene where ID in ("',paste(gene, collapse='","'),'")')
      p <- dbGetQuery(gene_pool_2019, query) %>% left_join(.,core_tight_2019) %>% 
        left_join(., gene_pool_2019 %>% tbl('gene_IDs') %>% as.tibble()) %>% 
        as.tibble()
    } else {
      query = paste0('select * from lsTPM_tx where ID in ("',paste(gene, collapse='","'),'")')
      p <- dbGetQuery(gene_pool_2019, query) %>% left_join(.,core_tight_2019) %>% 
        left_join(., gene_pool_2019 %>% tbl('tx_IDs') %>% as.tibble()) %>% 
        as.tibble()
    }
    p$Type <- p %>% select(contains('type')) %>% pull(1)
    p <- p %>% 
      filter(Sub_Tissue %in% tissue) %>% 
      mutate(Info = paste('SRA: ', 
                          sample_accession,
                          '\nStudy: ', 
                          study_title, '\n', 
                          gsub('\\|\\|', '\n', 
                               sample_attribute), 
                          sep =''),
             ID = gsub(' \\(', '\n(', ID)) %>% 
      ggplot(data=.,aes(x=Sub_Tissue,y=log2(value+1),colour=Tissue)) +
      geom_boxplot(alpha=0.5, outlier.shape = NA) + 
      geom_point_interactive(size=3, position = 'jitter', alpha=0.6, aes(tooltip=Info, shape = Type)) + 
      xlab('') + 
      facet_wrap(~ID, ncol=col_num) +
      theme_Publication(base_size = 15) + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.2)) +
      ggtitle('Box Plot of Pan-Human Gene Expression') +
      ylab("log2(TPM +1)") +
      scale_shape_manual(values=c(0:2,5,6,15:50))
    girafe(ggobj = p,width_svg = 12, 
           height_svg= max(10, (6 * (length(gene)/min(col_num,length(gene)))))) %>% 
      girafe_options(., opts_toolbar(position = NULL) )
    
  })
  output$boxPlot_gene <- renderGirafe({
    boxPlot_gene_func()
  })  
  
  # Gene fold change -----
  FC_gene_func <- eventReactive(input$pan_button_gene, {
    cat(file=stderr(), 'Gene FC call\n')
    db = input$Database
    gene <- input$ID
    tissue <- input$plot_tissue_gene
    col_num <- input$num_gene
    bench <- input$Bench_gene
    # })
    if (length(db) < 1 | length(gene) < 1 | length(tissue) < 1){
      showModal(modalDialog(title = "Fold Change Box Plot Warning",
                            "Have you specified at least one gene or tissue?", 
                            easyClose = T,
                            footer = NULL))
    }
    if (db == 'Gene 2017') {
      query = paste0('select * from lsTPM_gene where ID in ("',paste(gene, collapse='","'),'")')
      p <- dbGetQuery(gene_pool_2017, query) %>% left_join(.,core_tight_2017)
    } else if (db == 'Transcript 2017'){
      query = paste0('select * from lsTPM_TX where ID in ("',paste(gene, collapse='","'),'")')
      p <- dbGetQuery(gene_pool_2017, query) %>% left_join(.,core_tight_2017)
    } else if (db == 'Gene 2019'){
      query = paste0('select * from lsTPM_gene where ID in ("',paste(gene, collapse='","'),'")')
      p <- dbGetQuery(gene_pool_2019, query) %>% left_join(.,core_tight_2019)
    } else {
      query = paste0('select * from lsTPM_tx where ID in ("',paste(gene, collapse='","'),'")')
      p <- dbGetQuery(gene_pool_2019, query) %>% left_join(.,core_tight_2019)
    }
    plot_data <- p %>% 
      filter(Sub_Tissue %in% tissue) %>%
      group_by(ID) %>%
      mutate(Bench=ifelse(Sub_Tissue %in% bench, 1, 0), BenchValue=mean(log2(value[Bench==1]+1))) %>%
      group_by(ID, Sub_Tissue) %>%
      summarise(log2FC=mean(log2(value+1)) - mean(BenchValue))
    p <- ggplot(data=(plot_data),aes(x=Sub_Tissue,y=log2FC,fill=Sub_Tissue)) +
      geom_bar(stat = 'identity') + xlab('') + facet_wrap(~ID, ncol=col_num) +
      theme_Publication() + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.2)) +
      geom_hline(aes(yintercept=0,colour='Red')) +
      ggtitle('Fold Change (log2) of pan-human gene expression') +
      ylab("log2 Fold Change of Gene Expression")
    p
  })
  output$FC_gene <- renderPlot({
    FC_gene_func()
  }, height=function(){(400*length(input$ID))/min(input$num_gene,length(input$ID))})
  
  
  # Gene heatmap -------
  bulk_tissue_heatmap_func <- eventReactive(input$pan_button_gene, {
    cat(file=stderr(), 'Gene heatmap call\n')
    db = input$Database
    gene <- input$ID
    tissue <- input$plot_tissue_gene
    if (length(db) < 1 | length(gene) < 1 | length(tissue) < 2){
      showModal(modalDialog(title = "Heatmap Warning",
                            "Have you specified at least one gene and two tissues?", 
                            easyClose = T,
                            footer = NULL))
    }
    if (db == 'Gene 2017') {
      pool <- 'gene_pool_2017'
      table <- 'mean_rank_decile_gene'
    } else if (db == 'Transcript 2017'){
      pool <- 'gene_pool_2017'
      table <- 'mean_rank_decile_tx' 
    }
    else if (db == 'Gene 2019') {
      tissue <- trimws(tissue)
      pool <- 'gene_pool_2019'
      table <- 'mean_rank_decile_gene' 
    } else {
      tissue <- trimws(tissue)
      pool <- 'gene_pool_2019'
      table <- 'mean_rank_decile_tx' 
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
    superheat::superheat(id_matrix, heat.na.col = 'white',
                         scale = F, smooth.heat = F, #X.text = vals_id_matrix,
                         heat.pal = viridisLite::viridis(10),
                         #heat.pal.values = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                         X.text.col = text_col, #NA to leave blank
                         X.text.size = 4,
                         grid.hline.col = "white",
                         grid.vline.col = "white",
                         left.label.size = 0.8,
                         left.label.col = 'white',
                         left.label.text.alignment = 'right',
                         bottom.label.col = 'white',
                         bottom.label.text.alignment = 'right',
                         bottom.label.size = 2.4,
                         bottom.label.text.angle = 90,
                         legend = F)
  })
  output$bulk_tissue_heatmap <- renderPlot({
    bulk_tissue_heatmap_func()
  },  height=function(){(300*length(input$ID))/min(input$num_gene,length(input$ID))})
  
  # Gene info table ---------
  gene_info_maker <- eventReactive(input$pan_button_gene, {
    gene <- gsub(' (.*)','', input$ID) %>% unique()
    ensembl_base <- "https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g="
    genecards_base <- "https://www.genecards.org/cgi-bin/carddisp.pl?gene="
    omim_base <- "https://www.omim.org/search/?index=entry&sort=score+desc%2C+prefix_sort+desc&start=1&limit=10&search="
    gene_info <- gene %>% as.tibble() %>% rename(ID = value) %>% 
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
    } else if (db == 'Gene 2019'){
      tissue <- trimws(tissue)
      pool <- 'gene_pool_2019'
      table <- 'mean_rank_decile_gene'
    } else {
      tissue <- trimws(tissue)
      pool <- 'gene_pool_2019'
      table <- 'mean_rank_decile_tx'
    }
    get(pool) %>% tbl(table) %>% 
      filter(ID %in% gene, Sub_Tissue %in% tissue) %>%
      mutate(ID, Tissue = `Sub_Tissue`) %>%
      ungroup() %>%
      dplyr::select(ID, Tissue, Rank, Decile) %>%
      arrange(ID, Tissue) %>%
      as.tibble() %>% 
      DT::datatable(extensions = 'Buttons', rownames = F, options = list(
        pageLength = 20, dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
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
    } else if (db == 'Gene 2019'){
      query = paste0('select * from lsTPM_gene where ID in ("',paste(gene, collapse='","'),'")')
      plot_data <- dbGetQuery(gene_pool_2019,  query) %>%
        left_join(.,core_tight_2019) 
    } else {
      query = paste0('select * from lsTPM_tx where ID in ("',paste(gene, collapse='","'),'")')
      plot_data <- dbGetQuery(gene_pool_2019,  query) %>%
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
      do(tidy(t.test(.$value, tissue_subset$value))) %>%
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
  
  # modal for custom url ---------
  observeEvent(input$build_pan_url, {
    url = paste0("http://cyclops.nei.nih.gov/shiny/eyeIntegration/?Dataset=",
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
    if (input$single_cell_stat_type == 'Mean Gene Expression (across all cells, split by cell type)') {
      dec_column <- 'Decile_mean'
      rank_column <- 'Rank_mean'
    } else if (input$single_cell_stat_type == 'Percentage Cells Expressing Gene (across all cells, split by cell type and age (if available))') {
      dec_column <- 'Percentage Cells'
      rank_column <- 'Rank_cells'
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
        as.tibble() %>% 
        spread(Age, !!(dec_column)) 
      vals <- SC_pool %>% 
        tbl(table_name) %>% 
        filter(Gene == gene, `Total Count` > 10, `Cell Type` != 'Doublets', `Cell Type` != 'Red Blood Cells') %>%  
        select(Age, `Cell Type`, !!(rank_column)) %>% 
        as.tibble() %>% 
        spread(Age, !!(rank_column))
      
      mat_CH <- add_missing(mat) %>% as.matrix()
      vals_CH <- add_missing(vals, 'ND') %>% as.matrix()
    } else {
      mat_CH <- add_missing(SC_pool %>% tbl(table_name) %>% filter(Gene == gene, `Total Count` > 10) %>% select(`Cell Type`, !!(dec_column)) %>% data.frame(), expected_cols = 'Decile_mean') %>% as.matrix()
      colnames(mat_CH)[1] <- 'P14'
      mat_CH <- cbind(mat_CH, NA)
      
      vals_CH <- add_missing(SC_pool %>% tbl(table_name) %>% filter(Gene == gene, `Total Count` > 10) %>% select(`Cell Type`, !!(rank_column)) %>% data.frame(), expected_cols = 'Rank_mean') %>% as.matrix()
      colnames(vals_CH)[1] <- 'P14'
      vals_CH <- cbind(vals_CH, 'ND')
    }
    row.names(mat_CH) <- gsub('r M', 'r\nM', row.names(mat_CH)) 
    row.names(vals_CH) <- gsub('r M', 'r\nM', row.names(vals_CH)) 
    colors_CH <- mat_CH
    colors_CH <- ifelse(colors_CH < 3, 'gray', 'black')
    if (input$heatmap_overlay == 'None') {
      text_col <- NA
    } else{
      text_col <- colors_CH
    }
    superheat::superheat(mat_CH, heat.na.col = 'white',
                         scale = F, smooth.heat = F, X.text = vals_CH,
                         heat.pal = viridisLite::viridis(10),
                         heat.pal.values = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                         X.text.col = text_col, #NA to leave blank
                         X.text.size = 4,
                         grid.hline.col = "white",
                         grid.vline.col = "white",
                         left.label.size = 0.55,
                         left.label.col = 'white',
                         left.label.text.alignment = 'right',
                         bottom.label.col = 'white',
                         bottom.label.size = 0.05
    )
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
                                              aes(tooltip=paste0('Cell Type: ', `Cell Type`, '\nCell ID: ', `Cell ID`),
                                                  fill = `Cell Type`), 
                                              alpha = 0.6,
                                              size = 1.4,
                                              shape = 21,
                                              color = 'white') 
      } else{
        interactive <- geom_point_interactive(data = tsne_coords %>% 
                                                filter(age %in% input_age) %>%
                                                filter(`Cell ID` %in% samples_gene_up$`Cell ID`),
                                              aes(tooltip=paste0('Cell Type: ', `Cell Type`, '\nCell ID: ', `Cell ID`),
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
        theme_Publication(base_size=8) + 
        guides(fill=guide_legend(nrow = 4,byrow=TRUE)) + interactive
      girafe(code = print(p)) %>% girafe_options(., opts_toolbar(position = NULL) )
      
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
        theme_Publication(base_size=8) +
        guides(fill=guide_legend(nrow = 4,byrow=TRUE)) + interactive
      girafe(code = print(p)) %>% girafe_options(., opts_toolbar(position = NULL) )
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
        geom_point_interactive(size=20, alpha=0.2, aes(x=X1,y=X2,colour=Cluster, tooltip=Label)) +
        geom_point(data=tsne_plot %>% dplyr::select(X1,X2), size=5, alpha=0.2, colour='black', aes(x=X1,y=X2)) + 
        xlab('t-SNE 1') + ylab('t-SNE 2') +
        theme_minimal() +
        guides(colour = guide_legend(ncol = 3, 
                                     override.aes = list(size=10, alpha = 1)))
    } else if (dataSET == '2019'){
      all_tsne_plot_prepped <- gene_pool_2019 %>% 
        tbl('tSNE_bulk_RNA') %>% 
        filter(perplexity == perp)
      p <- ggplot(all_tsne_plot_prepped %>% as.tibble() %>% 
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
        geom_point_interactive(size=5, alpha=0.6, aes(shape=Origin, tooltip = Info)) +
        geom_label_repel(aes(label=Label), alpha=0.8, size=4, box.padding = unit(0.3, "lines")) + 
        theme_minimal() + 
        xlab('t-SNE 1') + 
        ylab('t-SNE 2') +
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
        filter(Comparison == input$de_comparison) %>% 
        as.tibble()
    } else if (input$diff_database == 'Transcript 2019') {
      de_data <- gene_pool_2019 %>%
        tbl('limma_DE_tx') %>%
        filter(Comparison == input$de_comparison) %>%
        as.tibble() %>%
        # left join to filter on gene or tx type (e.g. protein coding / miRNA / etc)
        left_join(., gene_pool_2019 %>% tbl('tx_IDs') %>% as.tibble())  %>%
        filter(transcript_type %in% input$gene_tx_type) %>% 
        rename(Class = transcript_type)
    } else {
      de_data <- gene_pool_2019 %>%
        tbl('limma_DE_gene') %>%
        filter(Comparison == input$de_comparison) %>%
        as.tibble() %>%
        # left join to filter on gene or tx type (e.g. protein coding / miRNA / etc)
        left_join(., gene_pool_2019 %>% tbl('gene_IDs') %>% as.tibble())  %>%
        filter(gene_type %in% input$gene_tx_type) %>% 
        rename(Class = gene_type)
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
        filter(Set == input$de_comparison, Test=='Up') %>% 
        dplyr::select(`GO ID`:Term,Ontology)
    } else {
      go <- gene_pool_2019 %>% 
        tbl('all_vs_all_go') %>% 
        filter(Set == input$de_comparison, Test=='Up') %>% 
        dplyr::select(ONTOLOGY:`p.adjust`)
    }
    go %>%  
      as.tibble()} %>%  
      DT::datatable(extensions = 'Buttons', options = list(
        dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))))
  output$go.table.down <- DT::renderDataTable(server = TRUE, {
    req(input$de_comparison)
    if (input$diff_database == 'Gene 2017') {
      go <- gene_pool_2017 %>% 
        tbl('all_vs_all_go') %>% 
        filter(Set == input$de_comparison, Test=='Down') %>% 
        dplyr::select(`GO ID`:Term,Ontology)
    } else {
      go <- gene_pool_2019 %>% 
        tbl('all_vs_all_go') %>% 
        filter(Set == input$de_comparison, Test=='Down') %>% 
        dplyr::select(ONTOLOGY:`p.adjust`)
    }
    go %>%  
      as.tibble()} %>%  
      DT::datatable(extensions = 'Buttons', options = list(
        dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))))
  
  # retina network -----------------------
  
  output$retina_network_mod <- renderVisNetwork({
    networks.retina.list[[paste0(input$retina_mod_color_vis_mod, "_k", input$retina_kNN)]]
  })
  output$retina_network_gene <- renderVisNetwork({
    req(input$retina_gene)
    new_mod_color = gene_pool_2017 %>% tbl('retina_gene_name_colors') %>% filter(id == input$retina_gene) %>% pull(Module_Color) 
    networks.retina.list[[paste0(new_mod_color, "_k", input$retina_kNN)]] %>%
      visOptions(nodesIdSelection = list(enabled = T, selected = input$retina_gene),
                 highlightNearest = list(enabled = T, degree = 0, hover = F, algorithm = "all"))
  })
  output$retina_table_mod <- renderDataTable(server = TRUE, {
    DT::datatable(gene_pool_2017 %>% tbl('edges_retina') %>% filter(Color == input$retina_mod_color_vis_mod) %>% as.tibble(), rownames = F, extensions = 'Buttons', options = list(
      dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('length'), digits = 5)
  })
  output$retina_table_gene <- renderDataTable(server = TRUE, {
    req(input$retina_gene)
    new_mod_color = gene_pool_2017 %>% tbl('retina_gene_name_colors') %>% filter(id == input$retina_gene) %>% pull(Module_Color) 
    #edges.retina.mod = edges.retina.list[[new_mod_color]]
    DT::datatable(gene_pool_2017 %>% 
                    tbl('edges_retina') %>% 
                    filter(Color == new_mod_color) %>% 
                    filter(Gene1 == input$retina_gene | Gene2 == input$retina_gene) %>% 
                    as.tibble(), 
                  rownames = F, extensions = 'Buttons', options = list(
                    dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('length'), digits = 5)
  })
  output$retina_connect_table_mod <- renderDataTable(server = TRUE, {
    DT::datatable(gene_pool_2017 %>% tbl('retina_mod_connect') %>% 
                    filter(`Module.Color` == input$retina_mod_color_vis_mod) %>% 
                    as.tibble(), extensions = 'Buttons', options = list(
                      dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('kWithin'), digits = 5)
  })
  output$retina_connect_table_gene <- renderDataTable(server = TRUE, {
    req(input$retina_gene)
    new_mod_color = gene_pool_2017 %>% tbl('retina_gene_name_colors') %>% filter(id == input$retina_gene) %>% pull(Module_Color) 
    DT::datatable(gene_pool_2017 %>% tbl('retina_mod_connect') %>% 
                    filter(`Module.Color` == new_mod_color) %>% 
                    as.tibble(), extensions = 'Buttons', options = list(
                      dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('kWithin'), digits = 5)
  })
  output$retina_GO_table_mod <- renderDataTable(server = TRUE, {
    DT::datatable(gene_pool_2017 %>% tbl('retina_network_GO') %>% filter(Color == input$retina_mod_color_vis_mod) %>% as.tibble(), 
                  extensions = 'Buttons', options = list(
                    dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
  })
  output$retina_GO_table_gene <- renderDataTable(server = TRUE, {
    req(input$retina_gene)
    new_mod_color = gene_pool_2017 %>% tbl('retina_gene_name_colors') %>% filter(id == input$retina_gene) %>% pull(Module_Color) 
    DT::datatable(gene_pool_2017 %>% tbl('retina_network_GO') %>% filter(Color == new_mod_color) %>% as.tibble(), 
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
    new_mod_color = gene_pool_2017 %>% tbl('rpe_gene_name_colors') %>% filter(id == input$rpe_gene) %>% pull(Module_Color) 
    networks.rpe.list[[paste0(new_mod_color, "_k", input$rpe_kNN)]] %>%
      visOptions(nodesIdSelection = list(enabled = T, selected = input$rpe_gene),
                 highlightNearest = list(enabled = T, degree = 0, hover = F, algorithm = "all"))
  })
  output$rpe_table_mod <- renderDataTable(server = TRUE, {
    DT::datatable(gene_pool_2017 %>% tbl('edges_rpe') %>% filter(Color == input$rpe_mod_color_vis_mod) %>% as.tibble(), rownames = F, extensions = 'Buttons', options = list(
      dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('length'), digits = 5)
  })
  output$rpe_table_gene <- renderDataTable(server = TRUE, {
    req(input$rpe_gene)
    new_mod_color = gene_pool_2017 %>% tbl('rpe_gene_name_colors') %>% filter(id == input$rpe_gene) %>% pull(Module_Color) 
    #edges.rpe.mod = edges.rpe.list[[new_mod_color]]
    DT::datatable(gene_pool_2017 %>% 
                    tbl('edges_rpe') %>% 
                    filter(Color == new_mod_color) %>% 
                    filter(Gene1 == input$rpe_gene | Gene2 == input$rpe_gene) %>% 
                    as.tibble(), 
                  rownames = F, extensions = 'Buttons', options = list(
                    dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('length'), digits = 5)
  })
  output$rpe_connect_table_mod <- renderDataTable(server = TRUE, {
    DT::datatable(gene_pool_2017 %>% tbl('rpe_mod_connect') %>% 
                    filter(`Module.Color` == input$rpe_mod_color_vis_mod) %>% 
                    as.tibble(), extensions = 'Buttons', options = list(
                      dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('kWithin'), digits = 5)
  })
  output$rpe_connect_table_gene <- renderDataTable(server = TRUE, {
    req(input$rpe_gene)
    new_mod_color = gene_pool_2017 %>% tbl('rpe_gene_name_colors') %>% filter(id == input$rpe_gene) %>% pull(Module_Color) 
    DT::datatable(gene_pool_2017 %>% tbl('rpe_mod_connect') %>% 
                    filter(`Module.Color` == new_mod_color) %>% 
                    as.tibble(), extensions = 'Buttons', options = list(
                      dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatSignif(c('kWithin'), digits = 5)
  })
  output$rpe_GO_table_mod <- renderDataTable(server = TRUE, {
    DT::datatable(gene_pool_2017 %>% tbl('rpe_network_GO') %>% filter(Color == input$rpe_mod_color_vis_mod) %>% as.tibble(), 
                  extensions = 'Buttons', options = list(
                    dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
  })
  output$rpe_GO_table_gene <- renderDataTable(server = TRUE, {
    req(input$rpe_gene)
    new_mod_color = gene_pool_2017 %>% tbl('rpe_gene_name_colors') %>% filter(id == input$rpe_gene) %>% pull(Module_Color) %>% as.tibble()
    DT::datatable(gene_pool_2017 %>% tbl('rpe_network_GO') %>% filter(Color == new_mod_color), 
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
      } else {
        table <- dbGetQuery(gene_pool_2019,  paste0('select * from lsTPM_tx where ID in ("',paste(input$table_gene, collapse='","'),'")') ) %>%
          left_join(.,core_tight_2019) 
      }
      table %>% filter(Tissue == input$table_tissue) %>%
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
