print('UI Start')
print(Sys.time())

library(shiny)
library(shinyjs)
library(ggplot2)
library(visNetwork)
library(colourpicker)
library(ggiraph)
library(shinythemes)
library(plotly)
# Data for PCA Visualization - created by the EiaD_build/scripts/pca_workup_for_eyeIntegration_app.Rmd script
load('./www/2023/eyeIntegration_2023_pca.Rdata')
# created by the EiaD_build/scripts/pca_workup_data_prep.R script
load('./www/2023/EiaD_pca_analysis_2023.Rdata')

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
               theme = shinytheme('cosmo'),
               selected = 'Overview',
               # Expression ---------------
               navbarMenu('Expression',
                          tabPanel('Pan-Tissue Plots',
                                   fluidPage(
                                     fluidRow(
                                       column(2,
                                              radioButtons('plot_type_gene',strong('Visualization:'),
                                                           choices = c('Box Plot', 'Heatmap', 'Table'),
                                                           selected = 'Box Plot'),
                                              conditionalPanel(condition = "input.plot_type_gene == 'Heatmap'",
                                                               checkboxGroupInput('heatmap_clustering_checkbox', strong('Clustering:'),
                                                                                  choices = list("Rows" = 1, "Columns" = 2))), br(), 
                                              selectizeInput('Database', strong('Dataset:'),
                                                             choices = c("Gene 2017", "Transcript 2017", "Gene 2019", "Transcript 2019", "DNTx v01", "Gene 2023", "Transcript 2023"), 
                                                             multiple = FALSE, selected = "Gene 2023"),
                                              selectizeInput('ID', strong('ID:'),
                                                             choices=NULL, multiple=TRUE),
                                              selectizeInput('plot_tissue_gene',strong('Tissues:'),
                                                             choices=NULL,multiple=TRUE),
                                              conditionalPanel(condition = "input.plot_type_gene != 'Heatmap' & input.Database != 'Gene 2023' & input.Database != 'Transcript 2023'",
                                                               numericInput('num_gene', strong('Number of Columns:'),
                                                                            value = 2, min = 1, max = 50)), br(), 
                                              conditionalPanel(condition = "input.plot_type_gene != 'Heatmap' & input.plot_type_gene != 'Table' & input.Database == 'Gene 2023'",
                                                               radioButtons('rotation', strong('Plot Orientation:'),
                                                                            choices = list("Samples as rows" = 1, "Samples as columns" = 2), 
                                                                            selected = 1)), br(), 
                                              conditionalPanel(condition = "input.plot_type_gene != 'Heatmap' & input.plot_type_gene != 'Table' & input.Database == 'Gene 2023'",
                                                               checkboxInput('points', 
                                                                             label = strong('Display Individual Sample Values'),
                                                                             value = FALSE)), br(),
                                              conditionalPanel(condition = "input.plot_type_gene != 'Heatmap' & input.plot_type_gene != 'Table' & input.Database == 'Transcript 2023'",
                                                               radioButtons('rotation', strong('Plot Orientation:'),
                                                                            choices = list("Samples as rows" = 1, "Samples as columns" = 2), 
                                                                            selected = 1)), br(), 
                                              conditionalPanel(condition = "input.plot_type_gene != 'Heatmap' & input.plot_type_gene != 'Table' & input.Database == 'Transcript 2023'",
                                                               checkboxInput('points', 
                                                                             label = strong('Display Individual Sample Values'),
                                                                             value = FALSE)), br(),
                                              actionButton('pan_button_gene','(Re)Draw Plot!', 
                                                           style='background-color: #3399ff; color: #ffffff'), br(), br(), 
                                              actionButton('build_pan_url','Build URL Shortcut', 
                                                           style='background-color: #808080; color: #ffffff'), br(), br(), br(), br(), br()
                                       ),
                                       column(7,
                                              conditionalPanel(condition = "input.pan_button_gene == 0", 
                                                               "Fill in the Gene and Tissue values and click the (RE)Draw Plot! button. 
                                                               It may take a few seconds for the plot and table to appear."), 
                                              conditionalPanel(condition = "input.plot_type_gene == 'Box Plot'",
                                                               girafeOutput('boxPlot_gene', height = '100%', width='100%')
                                              ),
                                              conditionalPanel(condition = "input.plot_type_gene == 'Fold Change'",
                                                               selectizeInput('Bench_gene','Select Reference Tissue(s):',
                                                                              choices = NULL ,multiple = TRUE),
                                                               plotOutput('FC_gene')
                                              ),
                                              conditionalPanel(condition = "input.plot_type_gene == 'Heatmap'",
                                                               plotOutput('bulk_tissue_heatmap')
                                              ),
                                              conditionalPanel(condition = "input.plot_type_gene == 'Table'",
                                                               div(DT::dataTableOutput('plot_table'), style='font-size:75%')
                                              )
                                       ),
                                       column(3,
                                              div(DT::dataTableOutput('gene_info'),style='font-size:75%'), br(), br(),
                                              conditionalPanel(condition = "input.plot_type_gene == 'Box Plot'",
                                                               div(DT::dataTableOutput('rankStats_gene'),style='font-size:75%'))
                                       )
                                     )
                                   ), br(), br(),
                                   fluidRow(includeHTML("www/footer.html"))
                          ),
                          # Single Cell Expression --------
                          tabPanel('Single Cell Expression',
                                   fluidPage(
                                     fluidRow(
                                       column(2,
                                              radioButtons('scplot_type_gene',strong('Visualization:'),
                                                           choices = c('Box Plot'),
                                                           selected = 'Box Plot'),
                                              selectizeInput('scGene', strong('ID:'),
                                                             choices=NULL, multiple=TRUE),
                                              selectizeInput('scmaturity', strong('Stage:'),
                                                             choices=NULL, multiple=TRUE),
                                              selectizeInput('scplot_tissue_gene',strong('Cells:'),
                                                             choices=NULL,multiple=TRUE),
                                              conditionalPanel(condition = "scplot_type_gene != 'Heatmap'",
                                                               radioButtons('sc_rotation', strong('Plot Orientation:'),
                                                                            choices = list("Cells as rows" = 1, "Cells as columns" = 2), 
                                                                            selected = 1)), 
                                              conditionalPanel(condition = "scplot_type_gene != 'Heatmap'",
                                                               checkboxInput('sc_points', 
                                                                             label = strong('Display Individual Sample Values'),
                                                                             value = FALSE)),
                                              actionButton('scpan_button_gene','(Re)Draw Plot!', 
                                                           style='background-color: #3399ff; color: #ffffff'), br(), br(), 
                                              actionButton('scbuild_pan_url','Build URL Shortcut', 
                                                           style='background-color: #808080; color: #ffffff'), br(), br(), br(), br(), br()
                                       ),
                                       column(7,
                                              conditionalPanel(condition = "input.scpan_button_gene == 0", 
                                                               "Fill in the Gene and Tissue values and click the (RE)Draw Plot! button. 
                                                               It may take a few seconds for the plot to appear."),
                                              conditionalPanel(condition = "input.scplot_type_gene == 'Box Plot'",
                                                               girafeOutput('scboxPlot_gene', height = '100%', width='100%')
                                              )
                                       )
                                     )
                                   ), br(), br(),
                                   fluidRow(includeHTML("www/footer.html"))
                          ),
                          # Pan - Eye and Body PCA ---------------
                          tabPanel('PCA Analysis Plotting',
                                   fluidPage(
                                     fluidRow(
                                       column(2,
                                              radioButtons('pca_visualization',strong('Visualization:'),
                                                           choices = c('eyeIntegration PCA Plot', 'Upload your own data'),
                                                           selected = 'eyeIntegration PCA Plot'),
                                              checkboxInput('GTEx_pca_data', 
                                                            label = 'Display GTEx Sample Data',
                                                            value = FALSE),
                                              checkboxInput('scRNA_pca_data', 
                                                            label = 'Display scRNA Sample Data',
                                                            value = FALSE),
                                              conditionalPanel(condition = "input.pca_visualization == 'eyeIntegration PCA Plot'",
                                                               checkboxInput('pc_top_genes', 
                                                                             label = 'Overlay top genes contributing to PC',
                                                                             value = FALSE)), br(),
                                              conditionalPanel(condition = "input.pca_visualization == 'Upload your own data'",
                                                               fileInput("user_samples", "Choose your CSV or TSV file for Upload",
                                                                         multiple = FALSE,
                                                                         accept = c(".csv", ".tsv", ".tsv.gz", ".csv.gz"))),
                                              conditionalPanel(condition = "input.pca_visualization == 'Upload your own data'",
                                                               selectizeInput('pca_user_plot_type', strong('How Would You Like to Display Your Data?'),
                                                                              choices = c("Faceted", "Overlayed"), 
                                                                              multiple = FALSE, selected = "Faceted")),
                                              conditionalPanel(condition = "input.pca_visualization == 'Upload your own data'",
                                                               textInput("user_given_input_project_name", label = strong("Enter Your Data Project Name:"), value = "User Generated Data"),
                                                               fluidRow(column(2, verbatimTextOutput("value")))), br(),
                                              selectizeInput('pca_component_one', strong('First PCA Component:'),
                                                             choices = c(eyeIntegration_2023_pca[[2]] %>% names() %>% head(21)), 
                                                             multiple = FALSE, selected = "PC1"),
                                              selectizeInput('pca_component_two', strong('Second PCA Component:'),
                                                             choices = c(eyeIntegration_2023_pca[[2]] %>% names() %>% head(21)), 
                                                             multiple = FALSE, selected = "PC2"),
                                              conditionalPanel(condition = "input.pca_visualization == 'eyeIntegration PCA Plot'",
                                                               actionButton('pca_button','(Re)Draw PCA Plot!', 
                                                                            style='background-color: #3399ff; color: #ffffff')),
                                              conditionalPanel(condition = "input.pca_visualization == 'Upload your own data'",
                                                               actionButton('user_generated_pca_button','(Re)Draw PCA Plot!', 
                                                                            style='background-color: #3399ff; color: #ffffff')), br(),
                                              conditionalPanel(condition = "input.pca_visualization == 'eyeIntegration PCA Plot'",
                                                               downloadButton('PCA_eyeIntegration_data','Download Plot Data',
                                                                              style='background-color: #3399ff; color: #ffffff')),
                                              conditionalPanel(condition = "input.pca_visualization == 'Upload your own data'",
                                                               downloadButton('PCA_ei_user_combined_data','Download Plot Data',
                                                                              style='background-color: #3399ff; color: #ffffff'))
                                       ),
                                       column(10,
                                              conditionalPanel(condition = "input.pca_button == 0 & input.pca_visualization == 'eyeIntegration PCA Plot'", 
                                                               "Select a first and second PCA component then click the (RE)Draw PCA Plot! button. 
                                                               It may take a few seconds for the plot to appear."),
                                              conditionalPanel(condition = "input.user_generated_pca_button == 0 & input.pca_visualization == 'Upload your own data'",
                                                               "Select a first and second PCA component then click the (RE)Draw PCA Plot! button. 
                                                               It may take a few seconds for the plot to appear."),
                                              conditionalPanel(condition = "input.pca_button != 0 & input.pca_visualization == 'eyeIntegration PCA Plot'",
                                                               plotlyOutput("eyeIntegration_pca_plot", height = "1080", width = "100%")
                                              ),
                                              conditionalPanel(condition = "input.user_generated_pca_button != 0 & input.pca_visualization == 'Upload your own data'",
                                                               plotlyOutput("user_pca_plot", height = "1080", width = "100%")
                                              )
                                       )
                                     )
                                   ), br(), br(),
                                   fluidRow(includeHTML("www/footer.html"))
                          ),
                          tabPanel('BP Level Expression', 
                                   
                                   fluidPage(
                                     fluidRow(column(width = 8, offset = 1, h2('UCSC Track for Base Pair Level Expression Coverage'))),
                                     fluidRow(column(width = 8, offset = 1, 'All links are external')),
                                     fluidRow(column(width = 8, offset = 1, img(src='ucsc_tracks_screenshot.png',width = '900px'))),
                                     fluidRow(column(width = 8, offset = 1, includeHTML("www/ucsc_tracks.html")))
                                   ))
               ),
               # navbarMenu("Find a Friend",
               #            tabPanel("Gene - Gene Euclidean Distance",
               #                     fluidPage(
               #                       fluidRow(strong("Euclidean distance of gene expression profiles across GTEx and Eye Tissues."),"In other words, the 100 genes with the most similar expression pattern across the EiaD dataset are returned. Co-expressed genes often have similar function and this is a quick way to look for genes with similar function to your favorite gene."),
               #                       fluidRow(
               #                         column(2,
               #                                selectizeInput('FaF_ID', strong('ID:'),
               #                                               choices = NULL, multiple = F)),
               #                         column(4,
               #                                div(DT::dataTableOutput('FaF_euc_dist'), style='font-size:75%'))
               #                       )
               #                     )
               #            )),
               navbarMenu('Data',
                          tabPanel('Pan-Tissue Bulk RNA-seq Table',
                                   fluidPage(
                                     fluidRow(
                                       column(2,
                                              selectizeInput('table_db',
                                                             strong('Dataset:'),
                                                             choices = c("Gene 2023", "Transcript 2023",
                                                                         "Gene 2019", "Transcript 2019",
                                                                         "Gene 2017", "Transcript 2017",
                                                                         "DNTx v01"),
                                                             selected = 'Gene 2023')),
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
                                       column(10, fluidRow(div(DT::dataTableOutput('table'), style='font-size:75%'))
                                       )))),
                          # Data Download ---------------
                          tabPanel('Data Download',
                                   fluidPage(
                                     br(),
                                     fluidRow(h2('As of 2023-04-21, we are do no yet have public facing download links for the new 2023 data. We will rectify this soon. Please email if you need
                                                 more immediate access to this data.')),
                                     br(),
                                     fluidRow(h3('References')),
                                     fluidRow(tags$a(href='ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.transcripts.fa.gz', 'gencode.v29.transcripts.fa.gz')),
                                     fluidRow(tags$a(href='ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz', 'gencode.v29.annotation.gtf.gz')),
                                     fluidRow(h3('Bulk Tissue Gene (or transcript(tx)) Raw Count Matrices')),
                                     
                                     fluidRow('Rows are genes, columns are samples, values are
                           raw counts as calculated by salmon.'),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2019_12_gene_counts_04.csv.gz',
                                                     '2019_12_gene_counts_04.csv.gz')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2019_12_transcript_counts_04.csv.gz',
                                                     '2019_12_transcript_counts_04.csv.gz')),
                                     fluidRow(h3('Bulk Tissue Gene (or transcript(tx)) Expression Matrices')),
                                     fluidRow('Rows are genes, columns are samples, values are in
                           length scaled Transcripts Per Million (TPM) as calculated by ',
                                              tags$a(href='https://github.com/mikelove/tximport/blob/master/R/tximport.R',
                                                     'tximport'), '.'),
                                     withMathJax(),
                                     fluidRow(column(2, '$$X = \\frac{count\\ of\\ reads\\ mapped\\ to\\ gene *
                           10^{3}}{gene\\ length\\ in\\ bp}$$')),
                                     fluidRow(column(2, '$$TPM = X \\ast \\frac{1}{\\sum X} \\ast 10^{6} $$')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2017_metadata_04.tsv.gz',
                                                     '2017_metadata_04.tsv.gz')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2017_gene_TPM_04.tsv.gz',
                                                     '2017_gene_TPM_04.tsv.gz')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2017_tx_TPM_04.tsv.gz',
                                                     '2017_tx_TPM_04.tsv.gz')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2019_metadata_04.tsv.gz',
                                                     '2019_metadata_04.tsv.gz')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2019_gene_TPM_04.tsv.gz',
                                                     '2019_gene_TPM_04.tsv.gz')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2019_tx_TPM_04.tsv.gz',
                                                     '2019_tx_TPM_04.tsv.gz')),
                                     fluidRow(h3(tags$em('De novo '), 'transcriptome data')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2019_DNTx_tx_TPM_00.tsv.gz',
                                                     '2019_DNTx_gene_TPM_01.tsv.gz')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2019_DNTx_tx_TPM_00.tsv.gz',
                                                     '2019_DNTx_tx_TPM_01.tsv.gz')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/DNTx_v00.fa.gz',
                                                     'DNTx_v01.fa.gz')),
                                     fluidRow(tags$a(href='https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/DNTx_v00.gtf.gz',
                                                     'DNTx_v01.gtf.gz')),
                                     fluidRow(h3('Everything')),
                                     fluidRow('All of the data and code for this entire web application can be
                           retrieved by following the simple directions ',
                                              tags$a(href='https://github.com/davemcg/eyeIntegration_app','here'), '.'),
                                     fluidRow(h3('Missing anything?')),
                                     fluidRow('If there\'s some data you want for easy download, ',
                                              tags$a(href = "mailto:mcgaugheyd@mail.nih.gov?Subject=eyeIntegration%20Comment",
                                                     "let me know"))),  br(), br())),
               # Information ---------------
               navbarMenu('Information',
                          tabPanel('Overview',
                                   fluidPage(
                                     fluidRow(column(width = 8, offset = 1, h2('eyeIntegration v2.00'))),
                                     fluidRow(column(width = 8, offset = 1, img(src='2023_eyeIntegration_Overview_drawIO.drawio.svg', width = 300))),
                                     fluidRow(column(width = 8, offset = 1, h2('Mission'))),
                                     fluidRow(column(width = 8, offset = 1,
                                                     "The human eye has several specialized tissues which direct, capture, and pre-process information to provide vision.
                                  RNA-seq gene expression analyses have been used extensively, for example, to profile specific eye tissues and in large consortium studies, like the GTEx project,
                                  to study tissue-specific gene expression patterning.", br(), br(), "However, there has not been an integrated study of multiple eye tissues expression patterning with other human
                                  body tissues.", br(), br(),
                                                     "We have collated publicly available (deposited between January 1st, 2019 and November 1st, 2023) healthy human RNA-seq datasets and a substantial subset of the GTEx project RNA-seq datasets and processed
                                  all of these samples in a consistent bioinformatic workflow. We use this fully integrated dataset to build informative visualizations, a novel PCA tool, and UCSC genome browser to provide the ophthalmic community with a powerful and quick means to formulate and test hypotheses on human 
                                  gene and transcript expression.", br(), br(),
                                                     
                                                     # h2("Basic Statistics"),
                                                     # fluidRow(column(6, img(src='sample_count_2023.svg', align='middle', width = 1000))),
                                                     # tableOutput('basic_stats'),
                                                     "We make these data, analyses, and visualizations available here with a powerful interactive web application.")),
                                     fluidRow(column(width = 8, offset = 1, h2('Attribution'))),
                                     fluidRow(column(width = 8, offset = 1, 'This project was conceived and implemented by',
                                                     tags$a(href = "mailto:mcgaugheyd@mail.nih.gov?Subject=eyeIntegration%20Comment", "David McGaughey"),
                                                     'in ', tags$a(href='https://nei.nih.gov/intramural/ogcsb','OGVFB'), '/',
                                                     tags$a(href='https://nei.nih.gov','NEI'), '/',
                                                     tags$a(href='https://www.nih.gov','NIH'), ', ',
                                                     ' in 2017. The retina and RPE gene networks along with their accompanying web pages were constructed by ',
                                                     tags$a(href='mailto:john.bryan@nih.gov', 'John Bryan'), '.', br(), br(),
                                                     'The 2019 automated pipeline datasets were built by ',
                                                     tags$a(href='mailto:vinay.swamy@nih.gov', 'Vinay Swamy'), '.', br(), br(),
                                                     'The 2023 update with the PCA analysis tool, its accompanying web page, and the ', tags$a(href="https://genome.ucsc.edu/s/parikhpp/Tissue%20Level%20BigWig%20Data", "tissue"),
                                                     '/', tags$a(href="https://genome.ucsc.edu/s/parikhpp/Sample%20Level%20BigWig%20Data", "sample"), 'level UCSC genome browser were built by ',
                                                     tags$a(href='mailto:pparikh1020@gmail.com', 'Prashit Parikh'),'. The new metadata curation schema along with the new samples is thanks to collaborative
                                                            work with ', tags$a(href = "https://jasonmiller.lab.medicine.umich.edu/", "Jason Miller"), ' and ', 
                                                     tags$a(href = 'https://prasov.lab.medicine.umich.edu', 'Lev Prasov'), '.', br(), br(),
                                                     
                                                     'Our analysis of the 2017 data in eyeIntegration has been published in Human Molecular Genetics. The manuscript is available ',
                                                     tags$a(href="https://academic.oup.com/hmg/article/27/19/3325/5042913",
                                                            "here"), '.
                                  If you use this resource in your research we would appreciate a citation.',
                                                     br(), br(),
                                                     'We also strongly encourage citation of the publications
                                  behind the datasets used in this resource. A full list can be found ',
                                                     tags$a(href="https://github.com/davemcg/eyeIntegration_app/blob/master/inst/citations.md",
                                                            "here"), ".")),
                                     
                                     fluidRow(column(width = 8, offset = 1, h2('Source Code'))),
                                     fluidRow(column(width = 8, offset = 1,
                                                     'The source code and data for this web application are available ', tags$a(href='https://gitlab.com/davemcg/Human_eyeIntegration_App', 'here'), '.')),
                                     fluidRow(column(width = 8, offset = 1, h2('Problems?'))),
                                     fluidRow(column(width = 8, offset = 1,
                                                     'First check the FAQ by clicking on ', strong('Information'), 'in the above header, then on ', strong('FAQs'), '.', br(), br(), 'Other issues can be reported two ways: ',
                                                     tags$a(href = "mailto:mcgaugheyd@mail.nih.gov?Subject=eyeIntegration%20Issue", "Email"), 'or ',
                                                     tags$a(href ='https://github.com/davemcg/eyeintegration_app/issues', 'GitHub Issue Tracker'))),br(), br(),
                                     fluidRow(includeHTML("www/footer.html")),
                                     br(), br()
                                   )),
                          # analysis -----------
                          tabPanel("Analysis and Extension", 
                                   fluidRow(column(width = 8, offset = 1, h2('Advanced Analysis'))),
                                   fluidRow(column(width = 8, offset = 1, "All links are external. Here we present a tutorial on how to use the data in EiaD 
                                   and recount3 to do custom differential testing. We also provide some brief guidance on how to use your own private data 
                                                   to run custom diff testing.")),
                                   fluidRow(column(width = 8, offset = 1, includeHTML("www/analyses.html")))),
                          # News ---------------
                          tabPanel('News',
                                   fluidPage(
                                     fluidRow(column(width = 8, offset = 1, h2('2023-05-22 | v2.01'))),
                                     fluidRow(column(width = 8, offset = 1, 'Swap out the recount3-based quant to a \"vanilla\" salmon quant. Has little effect on the data, but enables more straightforward outside comparison and brings back transcript level quant.')),
                                     fluidRow(column(width = 8, offset = 1, h2('2023-04-19 | v2.00'))),
                                     fluidRow(column(width = 8, offset = 1, 'Version 2.0! We introduce another huge set of updates, including a new 2023 dataset with 287 new eye samples, three new tissue categories, cell type level expression data, bulk RNA-seq expression boxplots to better express our new metadata, a new PCA tool with user-inputted data compatibility, and a UCSC genome browser for visualization of base-pair level expression counts. Click', tags$a(href="https://genome.ucsc.edu/s/parikhpp/Tissue%20Level%20BigWig%20Data", "here"), 'to view the tissue-level genome browser and', tags$a(href="https://genome.ucsc.edu/s/parikhpp/Sample%20Level%20BigWig%20Data", "here"), 'for the sample-level genome browser.')),
                                     fluidRow(column(width = 8, offset = 1, h2('2020-02-14 | v1.05'))),
                                     fluidRow(column(width = 8, offset = 1, 'Updated DNTx to v01. Removed v00 as we have made SUBSTANTIAL improvements to the precision and reliability of the results. We do not recommend v00 be used.')),
                                     fluidRow(column(width = 8, offset = 1, h2('2020-01-31 | v1.04'))),
                                     fluidRow(column(width = 8, offset = 1, 'Add raw counts to data download, as there were a handful of requests from users. Also fix data repo link from gitlab to github.')),
                                     fluidRow(column(width = 8, offset = 1, h2('2019-06-09 | v1.04'))),
                                     fluidRow(column(width = 8, offset = 1, 'Big update, which addresses (I hope) the comments from the reviewers of IOVS. More GTEx samples added per tissue. New GTEx tissues added (bladder, bone marrow, cervix uteri, fallopian tube, ovary, prostrate, testis, uterus, and vagina). Ratnapriya et al. AMD (MGS 1 == normal, 2,3,4 are increasing severity of AMD) retina dataset added. Modified lengthScaledTPM scores to adjust for tissue design and use mapping rate as covariate with limma\'s batchEffects() function. The differential expression test now uses mapping rate as covariate. On the UI side, now using fixed (consistent) color scheme for tissue in the box plots. Update summary stats on loading page with new tissues, numbers. Also showing GTEx count differences from 2017 to 2019.')),
                                     fluidRow(column(width = 8, offset = 1, h2('2019-04-26 | v1.03'))),
                                     fluidRow(column(width = 8, offset = 1, 'Added prototype ocular ', tags$em('de novo') ,' transcriptomes Vinay Swamy has built as a database option the pan-tissue visualizations and the data tables. We also make the ', tags$em('de novo') ,' gene models (GTF) and sequence (fasta) available for download in the Data -> Data Download section. Again, this is version 00 prototype data and ', tags$b('will'), ' change in the future. Depending on how much the project develops, we will expand eyeIntegration to show more information on the', tags$em('de novo') ,'transcript models or move parts of this project out into a new web site.')),
                                     fluidRow(column(width = 8, offset = 1, h2('2019-03-06 | v1.02'))),
                                     fluidRow(column(width = 8, offset = 1, 'Tweaked Retina Stem Cell / Organoid samples labelling. Added temporal heatmap for retina fetal and organoid time points. Modified bulk RNA-seq heatmap to use ComplexHeatmap, which handles long column names better. Optional row and clustering for the heatmaps.')),
                                     fluidRow(column(width = 8, offset = 1, h2('2019-01-16 | v1.01'))),
                                     fluidRow(column(width = 8, offset = 1, 'Fixed some tissue mislabels in EiaD 2019, removed unused sub-tissue, added eye sub-tissue vs human body tissue differential tests and GO term enrichments. Updated the workflow and 2017 to 2019 tables on the main loading page.')),
                                     fluidRow(column(width = 8, offset = 1, h2('2019-01-16 | v1.00'))),
                                     fluidRow(column(width = 8, offset = 1, 'Version 1.0! We introduce a huge set of updates, including a new 2019 dataset with 224 new eye samples, four new eye tissue categories, non-protein coding quantification, heatmap visualization, custom user shortcuts, quick gene information links, and easy data downloads. We will soon have a bioRxiv preprint describing the new automated pipeline that underlies the 2019 dataset. For users wanting to compare previous work done on eyeIntegration, the 2017 dataset is available as a versioned option.')),
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
                          # FAQ ----------
                          tabPanel('FAQs', fluidPage(
                            navlistPanel("eyeIntegration FAQs",
                                         tabPanel("How was the original (2017) data generated, acquired and processed?",
                                                  "See our publication:",
                                                  tags$a(href="https://academic.oup.com/hmg/advance-article/doi/10.1093/hmg/ddy239/5042913",
                                                         "https://academic.oup.com/hmg/advance-article/doi/10.1093/hmg/ddy239/5042913")),
                                         tabPanel("2017 Data Workflow",
                                                  img(src='2017_workflow.svg')),
                                         tabPanel("How was the 2019 data generated, acquired and processed?",
                                                  "We will soon have a pre-print describing the new 2019 automated data analysis workflow. The code-base for the build can be found ", tags$a(href="https://github.com/davemcg/eyeIntegration_data_build", "here.")),
                                         tabPanel("2019 Data Workflow",
                                                  img(src='2019_workflow.svg', width = 900)),
                                         tabPanel("How was the 2023 data generated, acquired and processed?",
                                                  "We will soon have a pre-print describing the new 2023 automated data analysis workflow. The code-base for the build can be found ", tags$a(href="https://github.com/davemcg/EiaD_build/tree/recount3", "here.")),
                                         tabPanel("2023 Data Workflow",
                                                  "We will upload a finalized 2023 data workflow soon.")),
                            navlistPanel("Pan-Tissue Expression FAQs",
                                         # How to use the boxplot split by iteration of EiaD
                                         tabPanel("How do I use the Pan-Tissue Expression 'Box Plot' section for data after 2019?",
                                                  "After selecting the 'Box Plot' radio button [1], you can choose a 2023 dataset [2]. 
                                                  Then, you can tweak which 'IDs' [3] and 'Tissues' [4] you would like to display by 
                                                  clicking in the respective boxes and starting to type (allowed values will auto-fill). 
                                                  You can also delete values by clicking on them and hitting the 'delete' key on your keyboard. 
                                                  The display of the box plot can be changed depending on whether you want your samples to be 
                                                  displayed as rows or columns [5]. Furthermore, the 'Display Individual Sample Values' checkbox [6] 
                                                  will enable you to hover your mouse over a data point and show the metadata for that particular sample [8]. 
                                                  When you are done tweaking these parameters, you can click the big blue '(Re)Draw Plot!' button [7] 
                                                  and wait a few seconds for the plot to appear.", br(), br(), 
                                                  img(src='pantissue_boxplot_2023.png', width = 900), br(), br()),
                                         tabPanel("How do I use the Pan-Tissue Expression ‘Box Plot' section for data from 2019 and prior?",
                                                  "After selecting the ‘Box Plot' radio button [1], you can choose a 2017, 2019, or DNTx dataset [2]. 
                                                  Then, you can tweak which ‘IDs' [3] and ‘Tissues' [4] you would like to display by clicking in the 
                                                  respective boxes and starting to type (allowed values will auto-fill). You can also delete values by 
                                                  clicking on them and hitting the ‘delete' key on your keyboard. You can change the display of the box 
                                                  plots by selecting a different value for the 'Number of Columns' field [5]. A lower number will squeeze 
                                                  more plots in each column. When you are done tweaking these parameters, you can click the big blue 
                                                  '(Re)Draw Plot!' button [6] and wait a few seconds for the plot to appear.", br(), br(), 
                                                  'If you mouse over a data point, you will get metadata about that particular sample [7].',
                                                  br(), br(), img(src='pantissue_boxplot_pre_2023.png', width = 900), br(), br()),
                                         # What the boxplot shows split by iteration of EiaD
                                         tabPanel("What data is displayed in the Pan-Tissue 'Box Plot' for data after 2019?",
                                                  "Each gene and tissue combination is given its own box. Depending on how the plot is oriented, 
                                                  one axis is log1p transformed z-counts, and the other axis contains the samples, colored by tissue. 
                                                  The right panel contains tables with external links to gene info [9], as well as the zCount values and 
                                                  rank of each gene in the chosen tissues (lower is more highly expressed).",
                                                  br(), br(), img(src='pantissue_boxplot_2023.png', width = 900), br(), br()),
                                         tabPanel("What Pan-Tissue Expression data is displayed in the ‘Box Plot' for data from 2019 and prior?",
                                                  "Each gene gets its own box. The y-axis is length scaled TPM (log2 transformed), and the x-axis is samples, 
                                                  colored by tissue. The right panel contains tables with external links to gene info [8], 
                                                  as well as the absolute TPM values and rank of each gene in the chosen tissues (lower is more highly expressed).", 
                                                  br(), br(), img(src='pantissue_boxplot_pre_2023.png', width = 900), br(), br()),
                                         # Other features which are uniform for data before and after 2019
                                         tabPanel("How do I use the Pan-Tissue Expression 'Heatmap' radio button?",
                                                  "After selecting the ‘Heatmap’ radio button [1], you can choose a 2017, 2019, or DNTx dataset [3]. 
                                                  Then, you can tweak which ‘IDs’ [4] and ‘Tissues’ [5] you would like to display by clicking in the respective boxes 
                                                  and starting to type (allowed values will auto-fill). You can also delete values by clicking on them and hitting the 
                                                  ‘delete’ key on your keyboard. The display of the heatmap can be changed depending on whether you want to cluster 
                                                  your samples by rows or columns by clicking the appropriate checkboxes [2]. When you are done tweaking these parameters, 
                                                  you can click the big blue '(Re)Draw Plot!' button [6] and wait a few seconds for the plot to appear", 
                                                  br(), br(), 
                                                  "The heatmap is an efficient way to display the expression of many genes and tissues. More yellow indicates higher expression, 
                                                  and further information about each chosen gene can be found by following the external links in the table to the right [7].", br(), br(),
                                                  img(src='pantissue_heatmap_pre_2023.png', width = 900), br(), br()),
                                         tabPanel("Why does the Heatmap for the 2023 datasets look different?",
                                                  "The 2023 heatmap operates at the tissue level, which results in more samples being present.", 
                                                  br(), br(), img(src='pantissue_heatmap_2023.png', width = 900), br(), br()),
                                         tabPanel("What is the 'Table' radio button?",
                                                  "This produces a table containing metadata for the gene and tissue combo selected by the user.", br(), br()),
                                          tabPanel("How do I use the Single-Cell Expression section?",
                                                  "First you select the genes [1], stages of development [2], and cell types [3] you would like to view by clicking in their respective boxes 
                                                  and starting to type (allowed values will auto-fill). You can delete values by clicking on them and hitting the 'delete' key on your keyboard. 
                                                  Next, select whether you would like to view your cell types as rows or columns in the 'Plot Orientation' section [4]. 
                                                  Furthermore, the ‘Display Individual Sample Values’ checkbox [5] will enable you to hover your mouse over a data point and show the metadata 
                                                  for that particular sample [7]. When you are done tweaking those parameters, click the big blue '(Re)Draw Plot!' button [6] and wait a few seconds.", 
                                                  br(), br(), img(src='singlecell_boxplot_2023.png', width = 900), br(), br()),
                                         tabPanel("What data is displayed in the Single-Cell Expression section?", 
                                                  "Each gene and developmental stage combination gets its own box. These features are further faceted between the back and front of the eye. 
                                                  Depending on how the plot is oriented, one axis is length scaled count value (log2 transformed), and the other axis contains the cell types. 
                                                  If hover data is toggled on [5], then each point is colored by an independent study and the size of the point is a log2 scaled percentage of cells 
                                                  that have detected expression of the gene.", br(), br(), img(src='singlecell_boxplot_2023.png', width = 900), br(), br()),
                                         tabPanel("How do I use the PCA Analysis Plotting tool using eyeIntegration’s built-in database?",
                                                  "This tool produces a plot of principal components from PCA (principal component analysis) conducted on our eyeIntegration 2.0 database. 
                                                  To begin, the user must select the 'eyeIntegration PCA Plot' radio button [1]. Then, you can select and unselect various checkboxes to 
                                                  include or exclude GTEx and single-cell RNA-seq data, as well as the visualization for which genes contribute most to your chosen principal 
                                                  components of interest [2]. Once you have decided which data you would like to visualize, you can select two principal components to plot [3]. 
                                                  When you are done tweaking those parameters, click the big blue '(Re)Draw Plot!' button [4] and wait a few seconds for the plot to appear. 
                                                  Since this plot was built using the ggplotly R package, you can hover your mouse over a point to see that sample’s metadata [6], 
                                                  and adjust the plot window using various scaling parameters provided on the top right of the plot [7]. Tissue types can be included and excluded 
                                                  by clicking directly on the name of the tissue within the legend on the right side of the plot [8]. 
                                                  Finally, if you would like to export the data used to make this plot, you can click the big blue 'Download Plot Data' button and download a CSV 
                                                  containing the data [5].", br(), br(), img(src='ei_pca_visualization.png', width = 900), br(), br()),
                                         tabPanel("How do I use my own data within the PCA Analysis Plotting tool?",
                                                  "This tool produces a plot of principal components from PCA (principal component analysis) conducted on our eyeIntegration 2.0 database and 
                                                  projects user-inputted data onto this PCA space. To begin, the user must select the 'Upload your own data' radio button [1]. 
                                                  Then, you can select and unselect checkboxes to include or exclude GTEx and single-cell RNA-seq data from eyeIntegration [2]. 
                                                  Once you have decided which samples you would like to include in visualization, you can click the ‘Browse…’ button and search your computer 
                                                  for a dataset to upload for projection [3]. This data should be in one of the following formats: csv, tsv, csv.gz, or tsv.gz. 
                                                  In this tool, the projected data can be either faceted against, or overlayed onto the existing eyeIntegration database. 
                                                  This can be selected from the dropdown menu under, ‘How Would You Like to Display Your Data?’ [4]. 
                                                  Finally, you can type in a name for your dataset to be included within the visualization [5], then select which principal components to plot [6]. 
                                                  When you are done tweaking those parameters, click the big blue '(Re)Draw Plot!' button [7] and wait a few seconds for the plot to appear. 
                                                  Since this plot was built using the ggplotly R package, you can hover your mouse over a point to see that sample’s metadata [9], 
                                                  and adjust the plot window using various scaling parameters provided on the top right of the plot [10]. Tissue types can be included and excluded 
                                                  by clicking directly on the name of the tissue within the legend on the right side of the plot [11]. 
                                                  Finally, if you would like to export the data used to make this plot, you can click the big blue 'Download Plot Data' button and download a CSV 
                                                  containing the data [8]. This data will include a scaled version of the original user data to plot alongside the eyeIntegration 2023 dataset.", 
                                                  br(), br(), strong("Faceted View"), br(), br(), img(src='user_faceted_pca_visualization.png', width = 900), br(), br(), 
                                                  strong("Overlayed View"), br(), br(), img(src='user_overlayed_pca_visualization.png', width = 900), br(), br())),
                            
                            navlistPanel("Data Table FAQ",
                                         tabPanel('What is the Pan-Tissue Bulk RNA-seq Data table showing?',
                                                  'The data table shows, for each gene and tissue set the user selects, the most important metadata for each sample.')),
                            
                            navlistPanel("Deprecated Features FAQ",
                                         tabPanel('What is the Deprecated Section?',
                                                  'These are all features from the 2019 iteration of eyeIntegration which have either been replaced with new features or deprecated due to limited use.'),
                                         
                                         "Differential Expression:",
                                         tabPanel("What is Differential Expression?",
                                                  "This is short for differential expression. We have pre-calculated 55+ differential expression tests. All eye tissue - origin pairs were compared to each other. We also have a synthetic human body set, made up of equal numbers of GTEx tissues (see manuscript, above, for more details). The word cloud displayed shows as many as the top 75 terms used in enriched GO terms in the selected comparison. The table data shows the actual GO terms. You can search for the comparison of your choice."),
                                         tabPanel("What are logFC, AveExp, t, P.Value, adj.Pval, and B in Differential Expression?",
                                                  'These are the values taken from the limma differential expression topTable() summary table. The following has been taken from the limma manual and edited to match parameters we used (https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf):',br(),br(),'A number of summary statistics are presented by topTable() for the top genes and the selected contrast. The logFC column gives the value of the contrast. Usually this represents a log2-fold change between two or more experimental conditions although sometimes it represents a log2-expression level. The AveExpr column gives the average log2-expression level for that gene across all the arrays and channels in the experiment. Column t is the moderated t-statistic. Column P.Value is the associated p-value and adj.P.Value is the p-value adjusted for multiple testing (False Discovery Rate corrected).', br(), br(), 'The B-statistic (lods or B) is the log-odds that the gene is differentially expressed. Suppose for example that B = 1.5. The odds of differential expression is exp(1.5)=4.48, i.e, about four and a half to one. The probability that the gene is differentially expressed is 4.48/(1+4.48)=0.82, i.e., the probability is about 82% that this gene is differentially expressed. A B-statistic of zero corresponds to a 50-50 chance that the gene is differentially expressed. The B-statistic is automatically adjusted for multiple testing by assuming that 1% of the genes, or some other percentage specified by the user in the call to eBayes(), are expected to be differentially expressed. The p-values and B-statistics will normally rank genes in the same order. In fact, if the data contains no missing values or quality weights, then the order will be precisely the same.'),
                                         
                                         "Mouse Single Cell Data:",
                                         tabPanel("Mouse Single Cell Retina Expression Datasets?",
                                                  'The Macosko data is a single-cell (~45,000) retina RNA-seq mouse P14 C57BL/6 dataset  from Mackosko and McCarroll\'s field defining ', tags$a(href="https://www.ncbi.nlm.nih.gov/pubmed/26000488","paper."), 'The cluster / cell type assignments are taken from ', tags$a(href="http://mccarrolllab.com/wp-content/uploads/2015/05/retina_clusteridentities.txt","here."), 'The Clark data is a 100,000 cell plus mouse retina RNA-seq time series dataset. Their pre-publication manuscript is on ', tags$a(href="https://www.ncbi.nlm.nih.gov/pubmed/26000488","bioRxiv. "), 'Data was pulled from ', tags$a(href="https://github.com/gofflab/developing_mouse_retina_scRNASeq", 'here.'), br(), br()),
                                         tabPanel("Mouse Single Cell Retina Expression Heatmaps / Tables?",
                                                  "To efficiently display a huge amount of information, expression across many individual cells is averaged by cell type, (if available) age, and gene. You can select the Macosko or Clark dataset [1], then one gene [2] to plot. The gene expression is displayed as a heatmap, with each row being a retina cell type (derived by the respective authors) and each column [4] is a time point, arranged from youngest to oldest. More yellow is higher expression [5].", br(), br(), img(src='sc_heatmap_top.png', width = 900), br(), img(src='sc_heatmap_bottom.png', width = 900), br(), br()),
                                         tabPanel("Mouse Single Cell Retina Expression Overlay?",
                                                  "You can add the rank of expression (or rank of percentage of cells with detectable expression of selected gene) with this radio", br(), br(), img(src='sc_heatmap_overlay.png', width = 900), br(), br()),
                                         tabPanel('What is the Mouse Single Cell Data table showing?',
                                                  'This data table set shows the data used to make the Mouse Single Cell Retina Expression plots in Gene Expression.'),
                                         
                                         "2D Clustering FAQ:",
                                         tabPanel("What is Bulk RNA-seq by Tissue?",
                                                  "This shows the t-SNE tissue clustering for the bulk human eye tissues along with the GTEx data-set. Hovering the mouse over each data point will show the metadata. Changing the perplexity will demonstrate how low values artificially create sub-groups while higher value (above 30 or so) largely recapitulate tissue type. It also demonstrates that the tissue clustering is stable at higher perplexities.", br(), br(), img(src='bulk_tsne.png', width = 900), br(), br()),
                                         tabPanel("What is Mouse Retina Single Cell...?",
                                                  "Each data point is a single cell from the ", tags$a(href="https://www.biorxiv.org/content/early/2018/07/30/378950","Macosko and McCarroll"), 'or the ', tags$a(href="https://www.ncbi.nlm.nih.gov/pubmed/26000488","Clark and Blackshaw"), '[2]. Dimensionality reduction with done with the t-SNE (Macosko) or UMAP (Clark) algorithm. Cluster assignments were taken from the respective papers. While we did the t-SNE on the Macosko data, the Clark authors provided the UMAP coordinates. The Clark dataset was generated across multiple time-points during development and thus, you can select time points of interest [4]. Only one gene can selected at a time [3], as it is very computationally expensive to plot many points. Points (cells) expressing the gene of interest are plotted in darker color [arrow]. Hovering over each point [6] will show the metadata for the cell.', br(), br(), img(src='sc_2d.png', width = 900), br(), br()),
                                         
                                         "Networks FAQ:",
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
                                                  "The edge table allows you to search for a gene and it returns all significant (edge length > 0.01) correlated genes ACROSS the entire network. Using the 'Connections to show' radio button, you can control whether only extramodular or intramodular (or both) connections are included in the table."))))),
               
               # Deprecated -----
               navbarMenu('Deprecated',
                          # Differential --------
                          tabPanel('Differential Expression',
                                   fluidPage(
                                     fluidRow(selectizeInput('diff_database', strong('Dataset:'),
                                                             choices = c('Gene 2017', 'Gene 2019', 'Transcript 2019', 'DNTx v01'),
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
                                     conditionalPanel(condition = "input.diff_database == 'Gene 2017' || input.diff_database == 'Gene 2019'",
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
                                     )), br(), br(),
                                   fluidRow(includeHTML("www/footer.html"))),
                          ## Retina Development Time Series -----
                          tabPanel('Retina Development Time Series',
                                   fluidPage(
                                     fluidRow(br()),
                                     fluidRow(strong('Heatmap of developing retina. Rows are genes, columns are age (in days) of retina or organoid.')),
                                     fluidRow(br()),
                                     fluidRow(selectizeInput('temporal_retina_heatmap_table', strong('Table:'),
                                                             choices = c('Gene 2019','Transcript 2019'))),
                                     fluidRow(selectizeInput('temporal_retina_heatmap_ID', strong('ID:'),
                                                             choices = NULL, multiple = T)),
                                     fluidRow(radioButtons('temporal_retina_heatmap_viz', strong('Type:'),
                                                           choices = c('Split by type','Merged'))),
                                     fluidRow(checkboxGroupInput('temporal_retina_heatmap_clustering', strong('Clustering:'),
                                                                 selected = 1,
                                                                 choices = list("Rows" = 1))), br(),
                                     actionButton('pan_button_temporal_heatmap','(Re)Draw Heatmap!',
                                                  style='background-color: #3399ff; color: #ffffff'), br(), br(),
                                     fluidRow(br()),
                                     fluidRow(
                                       column(10, strong('Gene Expression Levels across Developing Fetal and Organoid Retina (days)'),
                                              plotOutput('temporal_retina_heatmap', height = '2000px')))
                                   ), br(), br(),
                                   fluidRow(includeHTML("www/footer.html"))
                          ),
                          ## Mouse Single Cell Retina Expression --------
                          tabPanel('Mouse Single Cell Retina Expression',
                                   fluidPage(
                                     fluidRow(br()),
                                     fluidRow(strong('Mouse retina Single Cell Gene Expression Statistics by User Selected Gene(s)')),
                                     fluidRow('Please be patient with this section as calculating the density plots for all genes in across all cell types can take several seconds.'), 
                                     fluidRow(br()),
                                     fluidRow(actionButton('SC_density_pan_button','(Re)Draw Plot and Table!', style='background-color: #3399ff; color: #ffffff'), br(), br()),
                                     fluidRow(selectizeInput('SC_dataset', strong('Dataset:'),
                                                             choices = c('Clark et al. [pub. labels]', 'Clark et al. [bioRxiv labels]', 'Macosko et al.'))),
                                     fluidRow(selectizeInput('mGene', strong('Mouse Genes:'),
                                                             choices=NULL, multiple=FALSE)),
                                     fluidRow(radioButtons('single_cell_stat_type',strong('Summary Statistic Representation:'),
                                                           choices = c('Decile of Mean Gene Expression (across all cells, split by cell type)',
                                                                       'Percentage Cells Expressing Gene (across all cells, split by cell type and age (if available))'), 
                                                           selected = 'Decile of Mean Gene Expression (across all cells, split by cell type)')),
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
                          ),
                          ## Mouse Single Cell Retina RNA-seq Table ---------
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
                                     column(10,
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
                                            
                                     ))),
                          ## 2D Clustering ---------
                          tabPanel('2D Clustering: Bulk RNA-Seq by Tissue',
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
                                              girafeOutput('tsne', height = '1200px', width = '100%')))
                                   ), br(), br(),
                                   fluidRow(includeHTML("www/footer.html"))
                          ),
                          tabPanel('2D Clustering: Mouse Retina Single Cell RNA-Seq',
                                   fluidPage(
                                     fluidRow(strong('Visualization of single cell mouse retina sample expression patterning (tSNE or UMAP based)')),
                                     fluidRow(column(10, 'Each point is a individual cell with dimensionality reduction and cell labelling done by the publishing scientists. Similar cell types cluster together. Darker opacity are cells which express detectable levels of the user selected gene. The user can also select the minimum level of expression of the gene. Most gene\'s expression level ranges from 0 to 5 per cell. ', br(),br())),
                                     fluidRow(
                                       column(2,
                                              selectizeInput('SC_dataset_tsne', strong('Dataset:'),
                                                             choices = c('Clark et al.', 'Macosko et al.')),
                                              selectizeInput('mGene_tsne', strong('Mouse Gene:'),
                                                             choices=NULL, multiple=FALSE),
                                              selectizeInput('age_tsne', strong('Up to 4 Time Points:'),
                                                             choices=NULL, multiple=TRUE, options = list(maxItems = 4)),
                                              numericInput('min_single_cell_gene_count', strong('Gene Count Greater Than:'), value=0, min=0, max=15),
                                              actionButton('SC_tsne_pan_button','(Re)Draw!', style='background-color: #3399ff; color: #ffffff'), br(), br()),
                                       column(9,
                                              girafeOutput('single_cell_tsne_plot', height = '800px'))
                                     ))),
                          ## Eye Networks -------
                          tabPanel('Eye Networks: Retina Network',
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
                                              )))),
                                   br(), br(),br(),
                                   fluidRow(includeHTML("www/footer.html"))),
                          tabPanel('Eye Networks: Retina Network Edge Table',
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
                          tabPanel('Eye Networks: RPE Network',
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
                          tabPanel('Eye Networks: RPE Network Edge Table',
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
                                     )))
               ))),
  tags$style(HTML("
footer{
  background: #000000;
    margin: 0 auto;
  height:  70px;
  padding-top: 10px;
  padding-bottom: 10px;
  width: 100%;
  text-align: center;
}
footer ul{
  margin: 0px;
  padding:0;
}
footer li {
  /*  float: left;  */
    list-style-type: none;
  line-height: 4px;
  height: 11px;
  /* border-right: 1px solid #354052; */
  padding: 0px 10px;
  display:inline;
}
footer li a{
  text-decoration: none;
  color: #FFFFFF;
    font-size: 12px;
}

.footer {
  text-align: center;
}
.footer img {
  position: absolute;
  left: 30px;
  height: 50px;
}"))))