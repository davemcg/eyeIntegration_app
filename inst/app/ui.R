print('UI Start')
print(Sys.time())

library(shiny)
library(shinyjs)
library(ggplot2)
library(visNetwork)
library(colourpicker)
library(ggiraph)
library(shinythemes)

load('./www/2017/retina_colors.Rdata')
load('./www/2017/rpe_colors.Rdata')

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  useShinyjs(),
  div(
    id = "app-content",
    navbarPage(title = div(img(src = "eye.png", height = "43px", 
                               style = "position: relative; top: -10px; left: 5px")),
               windowTitle = 'eyeIntegration',
               footer=list(img(src='NEI_logo_mine.svg', height='75px')), 
               theme = shinytheme('cosmo'),
               selected = 'Overview',
               navbarMenu('Expression',
                          tabPanel('Pan-Tissue Plots',
                                   fluidPage(
                                     fluidRow(
                                       column(2,
                                              actionButton('pan_button_gene','(Re)Draw Plot!', 
                                                           style='background-color: #3399ff; color: #ffffff'), br(), br(),
                                              radioButtons('plot_type_gene',strong('Visualization:'),
                                                           choices = c('Box Plot','Fold Change', 'Heatmap'),
                                                           selected = 'Box Plot'),
                                              selectizeInput('Database', strong('Dataset:'),
                                                             choices = c("Gene 2017", "Transcript 2017", "Gene 2019", "Transcript 2019"), 
                                                             multiple = FALSE, selected = "Gene 2019"),
                                              selectizeInput('ID', strong('ID:'),
                                                             choices=NULL, multiple=TRUE),
                                              selectizeInput('plot_tissue_gene',strong('Tissues:'),
                                                             choices=NULL,multiple=TRUE),
                                              conditionalPanel(condition = "input.plot_type_gene != 'Heatmap'",
                                                               numericInput('num_gene', strong('Number of columns:'),
                                                                            value = 2, min = 1, max = 8)), br(), br(),
                                              actionButton('build_pan_url','Build URL Shortcut', 
                                                           style='background-color: #808080; color: #ffffff'), br(), br(), br(), br(), br()
                                       ),
                                       column(7,
                                              conditionalPanel(condition = "input.pan_button_gene == 0", 
                                                               "Fill in the Gene and Tissue values and click the (RE)Draw Plot! button. 
                                                               It may take a few seconds for the plot and table to appear."), 
                                              conditionalPanel(condition = "input.plot_type_gene == 'Box Plot'",
                                                               girafeOutput('boxPlot_gene', height = '1000px', width='100%')
                                              ),
                                              conditionalPanel(condition = "input.plot_type_gene == 'Fold Change'",
                                                               selectizeInput('Bench_gene','Select Reference Tissue(s):',
                                                                              choices = NULL ,multiple = TRUE),
                                                               plotOutput('FC_gene')
                                              ),
                                              conditionalPanel(condition = "input.plot_type_gene == 'Heatmap'",
                                                               plotOutput('bulk_tissue_heatmap')
                                              )
                                       ),
                                       column(3,
                                              div(DT::dataTableOutput('gene_info'),style='font-size:75%'), br(), br(),
                                              conditionalPanel(condition = "input.plot_type_gene == 'Fold Change'",
                                                               div(DT::dataTableOutput('basicStats_gene'),style='font-size:75%')),
                                              conditionalPanel(condition = "input.plot_type_gene != 'Fold Change'",
                                                               div(DT::dataTableOutput('rankStats_gene'),style='font-size:75%'))
                                       )
                                     )
                                   )
                          ),
                          tabPanel('Differential',
                                   fluidPage(
                                     fluidRow(selectizeInput('diff_database', strong('Dataset:'),
                                                             choices = c('Gene 2017', 'Gene 2019', 'Transcript 2019'),
                                                             selected = 'Gene 2019', 
                                                             multiple = FALSE)),
                                     conditionalPanel(condition = "input.diff_database != 'Gene 2017'",
                                                      fluidRow(selectizeInput('gene_tx_type', 
                                                                              strong('Class:'),
                                                                              choices = NULL,
                                                                              multiple = T
                                                      ))),
                                     fluidRow(selectizeInput('de_comparison',strong('Search for Differential Expression Comparison: '),
                                                             choices = NULL, #selected = NULL,
                                                             multiple = F,
                                                             options = list(placeholder = 'Type to search'))),
                                     fluidRow(DT::dataTableOutput('diff.exp')), br(),
                                     conditionalPanel(condition = "input.diff_database != 'Transcript 2019'",
                                                      fluidRow(strong('Word Clouds of most common enriched GO terms')), br(),
                                                      fluidRow(
                                                        column(uiOutput('comparison_up1'), uiOutput('word_cloud_image_up'), width=6, align='center'),
                                                        column(uiOutput('comparison_down1'), uiOutput('word_cloud_image_down'), width=6, align='center')
                                                      ),
                                                      fluidRow(strong('GO Term Enrichment')), br(),
                                                      fluidRow(uiOutput('comparison_up2')), br(),
                                                      fluidRow(DT::dataTableOutput('go.table.up')), br(),
                                                      fluidRow(uiOutput('comparison_down2')), br(),
                                                      fluidRow(DT::dataTableOutput('go.table.down'))
                                     ))),
                          tabPanel('Mouse Single Cell Retina Expression',
                                   fluidPage(
                                     fluidRow(br()),
                                     fluidRow(strong('Mouse retina Single Cell Gene Expression Statistics by User Selected Gene(s)')),
                                     fluidRow('Please be patient with this section as calculating the density plots for all genes in across all cell types can take several seconds.'), 
                                     fluidRow(br()),
                                     fluidRow(actionButton('SC_density_pan_button','(Re)Draw Plot and Table!', style='background-color: #3399ff; color: #ffffff'), br(), br()),
                                     fluidRow(selectizeInput('SC_dataset', strong('Dataset:'),
                                                             choices = c('Clark et al.', 'Macosko et al.'))),
                                     fluidRow(selectizeInput('mGene', strong('Mouse Genes:'),
                                                             choices=NULL, multiple=FALSE)),
                                     fluidRow(radioButtons('single_cell_stat_type',strong('Summary Statistic Representation:'),
                                                           choices = c('Mean Gene Expression (across all cells, split by cell type)',
                                                                       'Percentage Cells Expressing Gene (across all cells, split by cell type and age (if available))'), 
                                                           selected = 'Mean Gene Expression (across all cells, split by cell type)')),
                                     fluidRow(radioButtons('heatmap_overlay', strong('Overlay: '),
                                                           choices = c('None', 'Rank'),
                                                           selected = 'None')),
                                     fluidRow(br()),
                                     fluidRow(
                                       column(7,strong('Heatmap'),
                                              plotOutput('SC_gene_means_by_type_heatmap', height = '700px')),
                                       column(4,strong('Table Statistics'),
                                              div(DT::dataTableOutput('SC_gene_means_by_type_table'),style='font-size:75%')),
                                       fluidRow(br()),
                                       fluidRow(column(7, strong('ND = Not Detected'), offset = 1)),
                                       fluidRow(br())
                                     )
                                   )
                          )
               ),
               navbarMenu('2D Clustering',
                          tabPanel('Bulk RNA-Seq by Tissue',
                                   fluidPage(
                                     fluidRow(strong('tSNE Clustering of bulk human RNA-seq samples')),
                                     fluidRow('Each point is an RNA-seq(uenced) tissue. Similar tissues cluster together.', br(), br()),
                                     fluidRow('The perplexity may be viewed as a knob that sets the number of effective nearest neighbors. It is comparable with the number of nearest neighbors k that is employed in many manifold learners.', br(), br()),
                                     fluidRow('Please be patient, as over 1200 interactive data points are being loaded.', br(), br()),
                                     fluidRow(actionButton('bulk_tsne_button','(Re)Draw!', style='background-color: #3399ff; color: #ffffff'), br(), br()),
                                     fluidRow(
                                       column(2, selectizeInput('bulk_tsne_dataset', strong('Dataset:'), choices = c('2017', '2019'), selected = '2019')),
                                       column(4, numericInput('perplexity','Perplexity (5 - 50):', value=50, min=5, max=50, step = 5)
                                       ),
                                       column(10,
                                              ggiraphOutput('tsne', height = '1200px', width = '100%')))
                                   )
                          ),
                          tabPanel('Mouse Retina Single Cell RNA-Seq',
                                   fluidPage(
                                     fluidRow(strong('Visualization of single cell mouse retina sample expression patterning (tSNE or UMAP based)')),
                                     fluidRow(column(10, 'Each point is a individual cell with dimensionality reduction and cell labelling done by the publishing scientists. Similar cell types cluster together. Darker opacity are cells which express detectable levels of the user selected gene. The user can also select the minimum level of expression of the gene. Most gene\'s expression level ranges from 0 to 5 per cell. ', br(),br())),
                                     fluidRow(
                                       column(2,
                                              actionButton('SC_tsne_pan_button','(Re)Draw!', style='background-color: #3399ff; color: #ffffff'), br(), br(),
                                              selectizeInput('SC_dataset_tsne', strong('Dataset:'),
                                                             choices = c('Clark et al.', 'Macosko et al.')),
                                              selectizeInput('mGene_tsne', strong('Mouse Gene:'),
                                                             choices=NULL, multiple=FALSE),
                                              selectizeInput('age_tsne', strong('Up to 4 Time Points:'),
                                                             choices=NULL, multiple=TRUE, options = list(maxItems = 4)),
                                              numericInput('min_single_cell_gene_count', strong('Gene Count Greater Than:'), value=0, min=0, max=15)),
                                       column(9,
                                              girafeOutput('single_cell_tsne_plot', height = '800px'))
                                     )))),
               
               navbarMenu("Eye Networks",
                          tabPanel('Retina Network',
                                   fluidPage(
                                     fluidRow(
                                       column(1,
                                              radioButtons('retina_search_type',strong('Search type:'),
                                                           choices = c('By Gene','By Module'),
                                                           selected = 'By Gene')
                                       ),
                                       column(3,
                                              conditionalPanel(condition = "input.retina_search_type == 'By Module'",
                                                               colourInput(inputId = "retina_mod_color_vis_mod",
                                                                           label = strong("Module to visualize:"),
                                                                           NULL,
                                                                           returnName = TRUE,
                                                                           showColour = 'background',
                                                                           palette = "limited",
                                                                           allowedCols = retina_colors
                                                               )),
                                              conditionalPanel(condition = "input.retina_search_type == 'By Gene'",
                                                               selectizeInput('retina_gene',div(strong('Search for a specific gene to load a network:')), #style="color:red"),
                                                                              choices=NULL)
                                              )),
                                       column(3,
                                              selectInput(inputId = 'retina_kNN', 
                                                          label = strong('K-nearest genes to display (1, 3, or 5):'), 
                                                          choices = c(1, 3, 5), 
                                                          selected=3, 
                                                          multiple = FALSE)
                                       )),
                                     fluidRow(column(6, strong('Retina Network: '))),br(),
                                     fluidRow(
                                       column(1),
                                       column(6,
                                              conditionalPanel(condition = "input.retina_search_type == 'By Module'",
                                                               visNetwork::visNetworkOutput('retina_network_mod', width = "auto", height = "600px")
                                              ),
                                              conditionalPanel(condition = "input.retina_search_type == 'By Gene'",
                                                               visNetwork::visNetworkOutput('retina_network_gene', width = "auto", height = "600px")
                                              ))),
                                     fluidRow(column(8,img(src='2017/retina_module_num.svg', align='middle'))
                                     ),
                                     fluidRow(),
                                     fluidRow(
                                       column(4,
                                              conditionalPanel(strong('Table of pair-wise gene connection strength'),
                                                               condition = "input.retina_search_type == 'By Module'",
                                                               div(DT::dataTableOutput('retina_table_mod'),style='font-size:100%')
                                              ),
                                              conditionalPanel(strong('Table of pair-wise gene connection strength'),
                                                               condition = "input.retina_search_type == 'By Gene'",
                                                               div(DT::dataTableOutput('retina_table_gene'),style='font-size:100%')
                                              )),
                                       column(4,
                                              conditionalPanel(strong('Table of gene connectivities within module'),
                                                               condition = "input.retina_search_type == 'By Module'",
                                                               div(DT::dataTableOutput('retina_connect_table_mod'),style='font-size:100%')
                                              ),
                                              conditionalPanel(strong('Table of gene connectivities within module'),
                                                               condition = "input.retina_search_type == 'By Gene'",
                                                               div(DT::dataTableOutput('retina_connect_table_gene'),style='font-size:100%')
                                              ))),
                                     fluidRow(),
                                     fluidRow(
                                       column(8,
                                              conditionalPanel(strong('Table of enriched GO terms'),
                                                               condition = "input.retina_search_type == 'By Module'",
                                                               div(DT::dataTableOutput('retina_GO_table_mod'),style='font-size:100%')
                                              ),
                                              conditionalPanel(strong('Table of enriched GO terms'),
                                                               condition = "input.retina_search_type == 'By Gene'",
                                                               div(DT::dataTableOutput('retina_GO_table_gene'),style='font-size:100%')
                                              ))))),
                          tabPanel('Retina Network Edge Table',
                                   fluidPage(
                                     
                                     fluidRow(column(3,
                                                     selectizeInput('retina_gene_edges',div(strong('Search for a specific gene to load its edges:')), #style="color:red"),
                                                                    choices=NULL)),
                                              column(3,
                                                     radioButtons('retina_edge_show',strong('Connections to show:'),
                                                                  choices = c('Only Within Module',
                                                                              'Only Outside Module',
                                                                              'Both'), selected = 'Both'))
                                     ),
                                     fluidRow(DT::dataTableOutput('retina_full_edge_table')
                                     ))),
                          tabPanel('RPE Network',
                                   fluidPage(
                                     
                                     fluidRow(
                                       column(1,
                                              radioButtons('rpe_search_type',strong('Search type:'),
                                                           choices = c('By Gene','By Module'),
                                                           selected = 'By Gene')
                                       ),
                                       column(3,
                                              conditionalPanel(condition = "input.rpe_search_type == 'By Module'",
                                                               colourInput(inputId = "rpe_mod_color_vis_mod",
                                                                           label = strong("Module to visualize:"),
                                                                           NULL,
                                                                           returnName = TRUE,
                                                                           showColour = 'background',
                                                                           palette = "limited",
                                                                           allowedCols = rpe_colors
                                                               )),
                                              conditionalPanel(condition = "input.rpe_search_type == 'By Gene'",
                                                               selectizeInput('rpe_gene',div(strong('Search for a specific gene to load a network:')),# style="color:red"),
                                                                              choices=NULL)
                                              )),
                                       column(3,
                                              selectInput(inputId = 'rpe_kNN', label = strong('K-nearest genes to display (1, 3, or 5):'), choices = c(1, 3, 5), selected=3, multiple = FALSE)
                                       )),
                                     fluidRow(column(6, strong('RPE Network'))),br(),
                                     fluidRow(
                                       column(1),
                                       column(6,
                                              conditionalPanel(condition = "input.rpe_search_type == 'By Module'",
                                                               visNetwork::visNetworkOutput('rpe_network_mod', width = "auto", height = "600px")
                                              ),
                                              conditionalPanel(condition = "input.rpe_search_type == 'By Gene'",
                                                               visNetwork::visNetworkOutput('rpe_network_gene', width = "auto", height = "600px")
                                              ))),
                                     fluidRow(column(8,img(src='2017/rpe_module_num.svg', align='middle'))),
                                     fluidRow(),
                                     fluidRow(
                                       column(4,
                                              conditionalPanel(strong('Table of pair-wise gene connection strength'),
                                                               condition = "input.rpe_search_type == 'By Module'",
                                                               div(DT::dataTableOutput('rpe_table_mod'),style='font-size:100%')
                                              ),
                                              conditionalPanel(strong('Table of pair-wise gene connection strength'),
                                                               condition = "input.rpe_search_type == 'By Gene'",
                                                               div(DT::dataTableOutput('rpe_table_gene'),style='font-size:100%')
                                              )),
                                       column(4,
                                              conditionalPanel(strong('Table of gene connectivities within module'),
                                                               condition = "input.rpe_search_type == 'By Module'",
                                                               div(DT::dataTableOutput('rpe_connect_table_mod'),style='font-size:100%')
                                              ),
                                              conditionalPanel(strong('Table of gene connectivities within module'),
                                                               condition = "input.rpe_search_type == 'By Gene'",
                                                               div(DT::dataTableOutput('rpe_connect_table_gene'),style='font-size:100%')
                                              ))),
                                     fluidRow(),
                                     fluidRow(
                                       column(8,
                                              conditionalPanel(strong('Table of enriched GO terms'),
                                                               condition = "input.rpe_search_type == 'By Module'",
                                                               div(DT::dataTableOutput('rpe_GO_table_mod'),style='font-size:100%')
                                              ),
                                              conditionalPanel(strong('Table of enriched GO terms'),
                                                               condition = "input.rpe_search_type == 'By Gene'",
                                                               div(DT::dataTableOutput('rpe_GO_table_gene'),style='font-size:100%')
                                              ))))),
                          tabPanel('RPE Network Edge Table',
                                   fluidPage(
                                     fluidRow(column(3,
                                                     selectizeInput('rpe_gene_edges',div(strong('Search for a specific gene to load its edges:')), #style="color:red"),
                                                                    choices=NULL)),
                                              column(3,
                                                     radioButtons('rpe_edge_show',strong('Connections to show:'),
                                                                  choices = c('Only Within Module',
                                                                              'Only Outside Module',
                                                                              'Both'), selected = 'Both'))
                                     ),
                                     fluidRow(DT::dataTableOutput('rpe_full_edge_table')
                                     )))),
               navbarMenu('Data',
                          tabPanel('Pan-Tissue Bulk RNA-seq Table',
                                   fluidPage(
                                     fluidRow(
                                       column(2, 
                                              selectizeInput('table_db',
                                                             strong('Dataset:'),
                                                             choices = c("Gene 2017", "Transcript 2017",
                                                                         "Gene 2019", "Transcript 2019"),
                                                             selected = 'Gene 2019')),
                                       column(2,
                                              selectizeInput('table_tissue',
                                                             strong('Tissue: '),
                                                             multiple = TRUE,
                                                             choices = NULL)),
                                       column(2,
                                              selectizeInput('table_gene',
                                                             strong('Gene: '),
                                                             multiple = TRUE,
                                                             choices = NULL)),
                                       column(6,
                                              selectizeInput('table_columns',
                                                             strong('Columns: '),
                                                             choices = NULL,
                                                             multiple = TRUE)),
                                       fluidRow(DT::dataTableOutput('table')
                                       )))),
                          tabPanel('Mouse Retina Single Cell RNA-seq Table',
                                   fluidPage(
                                     fluidRow(
                                       column(2,
                                              selectInput('SC_datatable_dataset',
                                                          strong('Dataset:'),
                                                          choices = c('Clark et al.', 'Macosko et al.')),
                                              selectizeInput('sc_datatable_tissue',
                                                             strong('Cell Type:'),
                                                             choices = NULL))),
                                     fluidRow(radioButtons('single_cell_stat_type__datatable',strong('Summary Statistic Representation:'),
                                                           choices = c('Mean Gene Expression (across all cells, split by cell type)',
                                                                       'Percentage Cells Expressing Gene (across all cells, split by cell type)', 
                                                                       'Percentage of Cell Types (across only cells expressing the gene, split by cell type)'),
                                                           selected = 'Mean Gene Expression (across all cells, split by cell type)')),
                                     fluidRow(conditionalPanel(condition = "input.single_cell_stat_type__datatable == 'Mean Gene Expression (across all cells, split by cell type)'",
                                                               div(DT::dataTableOutput('SCtable_1'), style='font-size:75%'))),
                                     fluidRow(conditionalPanel(condition = "input.single_cell_stat_type__datatable == 'Percentage Cells Expressing Gene (across all cells, split by cell type)'",
                                                               div(DT::dataTableOutput('SCtable_2'), style='font-size:75%'))),
                                     fluidRow(conditionalPanel(condition = "input.single_cell_stat_type__datatable == 'Percentage of Cell Types (across only cells expressing the gene, split by cell type)'",
                                                               div(DT::dataTableOutput('SCtable_3'), style='font-size:75%')))
                                     
                                   )),
                          tabPanel('Data Download',
                                   fluidPage(
                                     fluidRow(h3('Bulk Tissue Gene (or transcript(tx)) Expression Matrices')),
                                     fluidRow('Rows are genes, columns are samples, values are in 
                                              length scaled Transcripts Per Million (TPM) as calculated by ', 
                                              tags$a(href='https://github.com/mikelove/tximport/blob/master/R/tximport.R', 
                                                     'tximport'), '.'),
                                     withMathJax(),
                                     fluidRow(column(2, '$$X = \\frac{count\\ of\\ reads\\ mapped\\ to\\ gene * 
                                              10^{3}}{gene\\ length\\ in\\ bp}$$')),
                                     fluidRow(column(2, '$$TPM = X \\ast \\frac{1}{\\sum X} \\ast 10^{6} $$')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2017_metadata.tsv.gz',
                                                     '2017_metadata.tsv.gz')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2017_gene_TPM.tsv.gz',
                                                     '2017_gene_TPM.tsv.gz')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2017_tx_TPM.tsv.gz',
                                                     '2017_tx_TPM.tsv.gz')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2019_metadata.tsv.gz',
                                                     '2019_metadata.tsv.gz')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2019_gene_TPM.tsv.gz',
                                                     '2019_gene_TPM.tsv.gz')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2019_tx_TPM.tsv.gz',
                                                     '2019_tx_TPM.tsv.gz')),
                                     fluidRow(h3('Everything')),
                                     fluidRow('All of the data and code for this entire web application can be 
                                              retrieved by following the simple directions ', 
                                              tags$a(href='https://github.com/davemcg/eyeIntegration_app','here'), '.'),
                                     fluidRow(h3('Missing anything?')),
                                     fluidRow('If there\'s some data you want for easy download, ',
                                              tags$a(href = "mailto:mcgaugheyd@mail.nih.gov?Subject=eyeIntegration%20Comment", 
                                                     "let me know"))),  br(), br())),
               navbarMenu('Information',
                          tabPanel('Overview',
                                   fluidPage(
                                     fluidRow(column(width = 8, offset = 1, h2('eyeIntegration v1.00'))),
                                     fluidRow(column(width = 8, offset = 1, img(src='view.png'))),
                                     fluidRow(column(width = 8, offset = 1, h2('Mission'))),
                                     fluidRow(column(width = 8, offset = 1,
                                                     "The human eye has several specialized tissues which direct, capture, and pre-process information to provide vision.
                                                     RNA-seq gene expression analyses have been used extensively, for example, to profile specific eye tissues and in large consortium studies, like the GTEx project,
                                                     to study tissue-specific gene expression patterning", br(), br(), "However, there has not been an integrated study of multiple eye tissues expression patterning with other human
                                                     body tissues.", br(), br(), 
                                                     "We have collated publicly available (January 12th, 2017 and January 1st, 2019) healthy human RNA-seq datasets and a substantial subset of the GTEx project RNA-seq datasets and processed
                                                     all in a consistent bioinformatic workflow. We use this fully integrated dataset to probe the relatedness and biological processes between the cornea, retina, RPE
                                                     (choroid), and the rest of the human tissues with differential expression, clustering, and GO term enrichment tools. We also leverage our large collection of retina
                                                     and RPE (choroid) tissues to build the first human weighted gene correlation networks and use them to highlight known biological pathways and eye gene disease
                                                     enrichment.", br(), br(), 
                                                     h2("Basic Statistics"), 
                                                     fluidRow(column(6, img(src='sample_count_2017_2019.svg', align='middle', width = 600))),
                                                     tableOutput('basic_stats'), 
                                                     "We make these data, analyses, and visualizations available here with a powerful interactive web application.")),
                                     fluidRow(column(width = 8, offset = 1, h2('Attribution'))),
                                     fluidRow(column(width = 8, offset = 1, 'This project was conceived and implemented by',
                                                     tags$a(href = "mailto:mcgaugheyd@mail.nih.gov?Subject=eyeIntegration%20Comment", "David McGaughey"),
                                                     'in ', tags$a(href='https://nei.nih.gov/intramural/ogcsb','OGVFB'), '/',
                                                     tags$a(href='https://nei.nih.gov','NEI'), '/',
                                                     tags$a(href='https://www.nih.gov','NIH'), '. ',
                                                     'The retina and RPE gene networks along with their accompanying web pages were constructed by ', 
                                                     tags$a(href='mailto:john.bryan@nih.gov', 'John Bryan.'), br(), br(), 
                                                     'The 2019 automated pipeline datasets were built by ',
                                                     tags$a(href='mailto:vinay.swamy@nih.gov', 'Vinay Swamy.'), br(), br(),
                                                     
                                                     'Our analysis of the data in eyeIntegration has been published in Human Molecular Genetics. The manuscript is available ', 
                                                     tags$a(href="https://academic.oup.com/hmg/article/27/19/3325/5042913",
                                                            "here"), '. 
                                                     If you use this resource in your research we would appreciate a citation.',
                                                     br(), br(), 
                                                     'We also strongly encourage citation of the publications 
                                                     behind the datasets used in this resource. A full list can be found ',
                                                     tags$a(href="https://gitlab.com/davemcg/Human_eyeIntegration_App/blob/master/citations.md",
                                                            "here."))),
                                     
                                     fluidRow(column(width = 8, offset = 1, h2('Source Code'))),
                                     fluidRow(column(width = 8, offset = 1,
                                                     'The source code and data for this web application is available ', tags$a(href='https://gitlab.com/davemcg/Human_eyeIntegration_App', 'here.'))),
                                     fluidRow(column(width = 8, offset = 1, h2('Problems?'))),
                                     fluidRow(column(width = 8, offset = 1,
                                                     'First check the FAQ by clicking on ', strong('Information'), 'in the above header, then on ', strong('FAQs'), br(), br(), 'Other issues can be reported two ways: ',
                                                     tags$a(href = "mailto:mcgaugheyd@mail.nih.gov?Subject=eyeIntegration%20Issue", "email"), 'or ',
                                                     tags$a(href ='https://gitlab.com/davemcg/Human_eyeIntegration_App/issues', 'GitLab Issue Tracker'))),br(), br(), br(), br()
                                   )),
                          tabPanel('News',
                                   fluidPage(
                                     fluidRow(column(width = 8, offset = 1, h2('2019-01-04 | v1.00'))),
                                     fluidRow(column(width = 8, offset = 1, 'Version 1.0! We introduce a huge set of updates, including a new 2019 dataset with 224 new eye samples, four new eye tissue categories, non-protein coding quantification, heatmap visualization, custom user shortcuts, and quick gene information links. We will soon have a bioRxiv preprint describing the new automated pipeline that underlies the 2019 dataset. For users wanting to compare previous work done on eyeIntegration, we make the 2017 dataset available as a versioned choice.')),
                                     fluidRow(column(width = 8, offset = 1, h2('2018-10-12 | v0.73'))),
                                     fluidRow(column(width = 8, offset = 1, 'Now using heatmaps for SC RNA-seq data.')),
                                     fluidRow(column(width = 8, offset = 1, h2('2018-10-04 | v0.72'))),
                                     fluidRow(column(width = 8, offset = 1, 'Engineering changes to make the site more responsive.')),
                                     fluidRow(column(width = 8, offset = 1, h2('2018-09-28 | v0.70'))),
                                     fluidRow(column(width = 8, offset = 1, 'Major changes! Can now select transcript level gene expression and Clark et al. biorXiv 2018 mouse retina time series scRNA-seq data added!')),
                                     fluidRow(column(width = 8, offset = 1, h2('2018-08-21 | v0.63'))),
                                     fluidRow(column(width = 8, offset = 1, 'Updated manuscript link to Human Molecular Genetics advance print and tweaked boxplot plotting size logic')),
                                     fluidRow(column(width = 8, offset = 1, h2('2018-01-09 | v0.62'))),
                                     fluidRow(column(width = 8, offset = 1, 'Table data added for single cell data')),                                     
                                     fluidRow(column(width = 8, offset = 1, h2('2017-11-09 | v0.61'))),
                                     fluidRow(column(width = 8, offset = 1, 'More granular single cell plots and tables in the Gene Expression section')),
                                     fluidRow(column(width = 8, offset = 1, h2('2017-09-29 | v0.60'))),
                                     fluidRow(column(width = 8, offset = 1, 'Mouse retina single cell data now available in the Gene Expression and 2D Clustering sections! Also the Pan-Tissue Expression section now has metadata for each sample point on mouseover!')),
                                     fluidRow(column(width = 8, offset = 1, h2('2017-05-23 | v0.51'))),
                                     fluidRow(column(width = 8, offset = 1, 'Now the user can export data from any table')),
                                     fluidRow(column(width = 8, offset = 1, h2('2017-05-17 | v0.5'))),
                                     fluidRow(column(width = 8, offset = 1, 'SQLite used on the backend for the largest data files to reduce initialization time and memory usage. FAQ section added')),
                                     fluidRow(column(width = 8, offset = 1, h2('2017-04-13 | v0.4'))),
                                     fluidRow(column(width = 8, offset = 1, 'Added full network edge tables for retina and RPE')),
                                     fluidRow(column(width = 8, offset = 1, h2('2017-04-08 | v0.3'))),
                                     fluidRow(column(width = 8, offset = 1, 'Network visualizations changed (gene names are nodes, layouts cleaned up a little), boxplot re-plot button more visible, on load now defaults to the info | overview page'))
                                   )),
                          tabPanel('FAQs', fluidPage(
                            navlistPanel("eyeIntegration FAQs",
                                         tabPanel("How was this data generated acquired and processed?",
                                                  "For the best description, see our advance publication:", tags$a(href="https://academic.oup.com/hmg/advance-article/doi/10.1093/hmg/ddy239/5042913","https://academic.oup.com/hmg/advance-article/doi/10.1093/hmg/ddy239/5042913")),
                                         tabPanel("Data workflow",
                                                  #tags$iframe(style="height:792px, width:612px",src= "workflow.pdf"))),
                                                  img(src='workflow.svg'))),
                            navlistPanel("Pan-Tissue Expression FAQs", 
                                         tabPanel("How do I use the Pan-Tissue Expression section?",
                                                  "You can tweak the 'Genes' [1] and 'Tissues' [2] by clicking in them and starting to type (allowed values will auto-fill). You can also delete values by clicking on them and hitting the 'delete' key on your keyboard. You can tweak the display of the box plots a bit by changing the 'Number of columns' field [3]. A higher number will squeeze more plots in each row. When you are done tweaking those parameters, click the big blue '(Re)Draw Plot!' button and wait a couple of seconds.", br(), br(), 'If you mouse over a data point, you will get metadata about that particular sample.', br(), br(), img(src='pantissue_screenshot.png'), br(), br()),
                                         tabPanel("What Pan-Tissue Expression data is displayed?",
                                                  "Each gene gets its own box. The y-axis is length scaled TPM (log2 transformed). The x axis is samples, colored by tissue. The right panel [4] is a table of the absolute TPM values and the rank of the gene in the particular tissue (lower is more highly expressed).", br(), br(), img(src='pantissue_screenshot.png'), br(), br()),
                                         tabPanel("How do I interpret the 'Decile' column in the data table?",
                                                  "To give a rough sense of how highly expressed the gene is in the tissue, the decile of expression is given in [4]; 10 is the highest decile of expression and 1 is the lowest.", br(), br(), img(src='pantissue_screenshot.png'), br(), br()),
                                         tabPanel("What is that 'Fold Change' radio button?",
                                                  "This is a very simple differential expression test. When you click this radio button [1] the view changes, where the expression is being shown *relative* to the reference tissues [2]. By default it is all of the tissues. This field above the bar plots can also be edited by you. The data table on the right will then display log2 fold change, average expression, and the p-value (from a t-test) of the differential expression.", br(), br(), img(src='pantissueFC_screenshot.png'), br(), br()),
                                         tabPanel("Custom Gene / Tissue combination via the URL",
                                                  "If you use the Pan-Tissue Boxplot feature a lot, you may find it frustrating to have to change the default genes (ABCA4 and TYRP1) and tissues to something you like to use a reference. We have added the ability to use a custom url to load in the genes and tissues of your choice. Here's a complex example that loads in four genes and many tissues:", br(),br(),tags$a(href="https://eyeintegration.nei.nih.gov/?Gene=SOX2,PAX6,MITF,VSX2&Tissue=Liver,Lung,_Brain_-_Cortex_,Retina_-_Adult_Tissue,Retina_-_Stem_Cell_Line,RPE_-_Cell_Line,RPE_-_Adult_Tissue,RPE_-_Fetal_Tissue,RPE_-_Stem_Cell_Line_Cells_-_Transformed_fibroblasts_,_Cells_-_EBV-transformed_lymphocytes_,_Whole_Blood_,_Brain_-_Cerebellum_,_Pituitary_,_Brain_-_Hypothalamus_,_Skin_-_Sun_Exposed_,_Skin_-_Not_Sun_Exposed_(Suprapubic)_","https://eyeintegration.nei.nih.gov/?Gene=SOX2,PAX6,MITF,VSX2&Tissue=Liver,Lung,_Brain_-_Cortex_,Retina_-_Adult_Tissue,Retina_-_Stem_Cell_Line,RPE_-_Cell_Line,RPE_-_Adult_Tissue,RPE_-_Fetal_Tissue,RPE_-_Stem_Cell_Line_Cells_-_Transformed_fibroblasts_,_Cells_-_EBV-transformed_lymphocytes_,_Whole_Blood_,_Brain_-_Cerebellum_,_Pituitary_,_Brain_-_Hypothalamus_,_Skin_-_Sun_Exposed_,_Skin_-_Not_Sun_Exposed_(Suprapubic)_"),br(),br(),"Genes are added after the ?Gene=. Separate each gene by a , . Tissues are added after the &Tissue= . Replace the spaces in the tissue name with a _ . For non-eye tissues, place a _ at the beginning and end. You must type the tissue name exactly as shown and get the Upper and lower cases correct.")),
                            navlistPanel("Mouse Single Cell Retina Expression FAQs", 
                                         tabPanel("Where did this data come from?",
                                                  'This single-cell (~45,000) retina RNA-seq mouse P14 C57BL/6 data comes from Mackosko and McCarroll\'s field defining ', tags$a(href="https://www.ncbi.nlm.nih.gov/pubmed/26000488","paper."), 'The cluster / cell type assignments are taken from ', tags$a(href="http://mccarrolllab.com/wp-content/uploads/2015/05/retina_clusteridentities.txt","here"), br(), br()),
                                         tabPanel("How do I interpret the density plots?",
                                                  'The plot is a log2(Gene Value + 1e-6) + 20 density plot showing three views. First, distributions of gene expression, split by cell type, with the mean expression of the user selected gene(s) overlaid. Second, the The color bars show the expression of user-defined gene(s). The table on the right shows the mean expression for the selected gene(s), along with the the rank (1 is the most expressed gene in the tissue), and the decile of expression the gene falls in (10 is the most expressed decile)', br(), br(), img(src='single_cell_density_mean_gene.png'), br(), br()),
                                         tabPanel("How do I interpret the second row - Percentage of Cell Types?",
                                                  "The figure is a density plot of the percentage, for each gene, of the distribution of the twelve cell types. Each of the ~45,000 cells was assigned a cell type. So, for example we see that Pax6, 62% of the cells with detectable Pax6 expression are classified as Amacrine cells. This is much higher than most (all but 573) genes. We also see from the density plots that around 50% of genes are over 50% Rod expressed and the other 50% are not expressed in the Rods at all (bimodal distribution). We also see that Astrocytes, Fibroblasts, Microglia, Pericytes are very rarely common classifications for each gene (strong left spike).", br(), br(), img(src='single_cell_density_percentage.png'), br(), br())),
                            navlistPanel("Eye Plot FAQ",
                                         tabPanel("Where did this go?",
                                                  "Originally, this was to re-display the same data as the Pan-Tissue section for all of the eye tissues with overlays for each data. But now with the Pan-Tissue Expression section getting overlays for each datapoint, this section is superfluous.", br(), br(), img(src='eyeplot_screenshot.png'))),
                            navlistPanel("2D Clustering FAQ",
                                         tabPanel("What is Bulk RNA-seq by Tissue?",
                                                  "This shows the t-SNE tissue clustering for the bulk human eye tissues along with the GTEx data-set. Hovering the mouse over each data point will show the metadata. Changing the perplexity will demonstrate how low values artificially create sub-groups while higher value (above 30 or so) largely recapitulate tissue type. It also demonstrates that the tissue clustering is stable at higher perplexities."),
                                         tabPanel("What is Mouse Retina Single Cell...?",
                                                  "Each data point is a single cell from the ", tags$a(href="https://www.ncbi.nlm.nih.gov/pubmed/26000488","Macosko and McCarroll paper."), 'Clustering with done with the t-SNE algorithm on the 13,000 cells which had more than 900 unique genes expressed. Cluster assignments were taken from Macosko et al., same as the Mouse Single Cell Expression section.')),
                            navlistPanel("DiffExp Clustering FAQs",
                                         tabPanel("What is DiffExp showing?",
                                                  "This is short for differential expression. We have pre-calculated 55 differential expression tests. All eye tissue - origin pairs were compared to each other. We also have a synthetic human body set, made up of equal numbers of GTEx tissues (see manuscript, above, for more details). The word cloud displayed shows as many as the top 75 terms used in enriched GO terms in the selected comparison. The table data shows the actual GO terms. You can search for the comparison of your choice."),
                                         tabPanel("What are logFC, AveExp, t, P.Value, adj.Pval, and B in DiffExp?",
                                                  'These are the values taken from the limma differential expression topTable() summary table. The following has been taken from the limma manual and edited to match parameters we used (https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf):',br(),br(),'A number of summary statistics are presented by topTable() for the top genes and the selected contrast. The logFC column gives the value of the contrast. Usually this represents a log2-fold change between two or more experimental conditions although sometimes it represents a log2-expression level. The AveExpr column gives the average log2-expression level for that gene across all the arrays and channels in the experiment. Column t is the moderated t-statistic. Column P.Value is the associated p-value and adj.P.Value is the p-value adjusted for multiple testing (False Discovery Rate corrected).', br(), br(), 'The B-statistic (lods or B) is the log-odds that the gene is differentially expressed. Suppose for example that B = 1.5. The odds of differential expression is exp(1.5)=4.48, i.e, about four and a half to one. The probability that the gene is differentially expressed is 4.48/(1+4.48)=0.82, i.e., the probability is about 82% that this gene is differentially expressed. A B-statistic of zero corresponds to a 50-50 chance that the gene is differentially expressed. The B-statistic is automatically adjusted for multiple testing by assuming that 1% of the genes, or some other percentage specified by the user in the call to eBayes(), are expected to be differentially expressed. The p-values and B-statistics will normally rank genes in the same order. In fact, if the data contains no missing values or quality weights, then the order will be precisely the same.'),
                                         br()),
                            navlistPanel("GO Heatmap FAQ",
                                         tabPanel('Why is this heatmap so long?',
                                                  "This is a massive expansion of Figure 3 from our manuscript. In the manuscript we could only display the top 80 Gene Ontology (GO) terms - there is no such limit on a computer, and here we display over 2000 GO terms with a FDR-corrected p-value < 0.05. Hovering over the heatmap will show precisely the GO term and comparison displayed. Yellow is more significant, blue is less.")),
                            navlistPanel("Networks FAQ",
                                         tabPanel('What is a WGCNA network?',
                                                  'This is a weighted gene expression correlation network. The gene expression information for all retina or all RPE tissues is used to identify gene pairs whose expression is correlated with each other. All of the pair-wise correlations are assessed to build a network of interactions.'),
                                         tabPanel('How do I use the network graphic thing?',
                                                  'We imagine the most common use is to search for your gene of interest (GOI). Simply type your GOI into the search box [1]. If it is not in the network, then the name will not appear. After selecting the GOI, the network will reload to display the module the gene is in, as well as several of the most correlated partners. You can adjust the number of displayed correlated genes by changing the K-nearest genes panel [2]. Hovering over a gene name in the network will display GO terms for the gene [3].', br(), br(), img(src='network_1_screenshot.png')),
                                         tabPanel('What is the direction of the gene to gene interaction?',
                                                  'Unfortunately, we have no way of knowing this. Since the network algorithms use correlations, the gene to gene interactions have no directional information.'),
                                         tabPanel('What are the figure and tables below the network visualization?',
                                                  'The count plot [1] simply shows the number of genes in each module. A gene can only be in one module. The pair-wise gene connection strength [2] shows the strongest gene partners for the selected gene. If a module search is selected, then this table shows all gene to gene edge connection strengths (higher is stronger) in the module. The gene connectivity table [3] shows the kWithin metric for each gene in the module, which denotes how connected (and important) the gene is across the module. The GO term table [4] shows the significant GO terms for the genes in the module. This allows you to get a sense of the function of the module.', br(), br(), img(src='network_2_screenshot.png'),
                                                  br(), br(), br()),
                                         tabPanel('What is the edge table showing?',
                                                  "The edge table allows you to search for a gene and it returns all significant (edge length > 0.01) correlated genes ACROSS the entire network. Using the 'Connections to show' radio button, you can control whether only extramodular or intramodular (or both) connections are included in the table.")),
                            navlistPanel("Data Table FAQ",
                                         tabPanel('What is the Pan-Tissue Bulk RNA-seq Data table showing?',
                                                  'The data table shows, for each gene and tissue set the user selects, the most important metadata for each sample.'),
                                         tabPanel('What is the Single Cell Data table showing?',
                                                  'This data table set shows the data used to make the Mouse Single Cell Retina Expression plots in Gene Expression.'))
                          )))))))

