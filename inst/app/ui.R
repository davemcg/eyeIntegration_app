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
                                                             choices = c("Gene 2017", "Transcript 2017", "Gene 2019", "Transcript 2019", "DNTx v01", "Gene 2022"), 
                                                             multiple = FALSE, selected = "Gene 2022"),
                                              selectizeInput('ID', strong('ID:'),
                                                             choices=NULL, multiple=TRUE),
                                              selectizeInput('plot_tissue_gene',strong('Tissues:'),
                                                             choices=NULL,multiple=TRUE),
                                              conditionalPanel(condition = "input.plot_type_gene != 'Heatmap' & input.Database != 'Gene 2022'",
                                                               numericInput('num_gene', strong('Number of columns:'),
                                                                            value = 2, min = 1, max = 50)), br(), 
                                              conditionalPanel(condition = "input.plot_type_gene != 'Heatmap' & input.Database == 'Gene 2022'",
                                                               radioButtons('rotation', strong('Plot Orientation:'),
                                                                            choices = list("Samples as rows" = 1, "Samples as columns" = 2), 
                                                                            selected = 1)), br(), 
                                              conditionalPanel(condition = "input.plot_type_gene != 'Heatmap' & input.Database == 'Gene 2022'",
                                                               checkboxInput('points', 
                                                                            label = strong('Display individual sample values'),
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
                                              div(DT::dataTableOutput('gene_info'),style='font-size:75%')
                                       )
                                     ),
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
                                              selectizeInput('scmaturity', strong('Maturity Stage:'),
                                                             choices=NULL, multiple=TRUE),
                                              selectizeInput('scplot_tissue_gene',strong('Tissues:'),
                                                             choices=NULL,multiple=TRUE),
                                              conditionalPanel(condition = "scplot_type_gene != 'Heatmap'",
                                                               numericInput('scnum_gene', strong('Number of columns:'),
                                                                            value = 2, min = 1, max = 8)), br(), 
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
                                                               girafeOutput('scboxPlot_gene', height = '1000px', width='100%')
                                              ),
                                       )
                                     )
                                   ), br(), br(),
                                   fluidRow(includeHTML("www/footer.html"))
                          ),
                          # Eye PCA ---------------
                          tabPanel('PCA Analysis Plotting',
                                   fluidPage(
                                     fluidRow(
                                       column(2,
                                              radioButtons('pca_data_format',strong('Visualization:'),
                                                           choices = c('PCA Plot'),
                                                           selected = 'PCA Plot'),
                                              selectizeInput('pca_component_one', strong('First PCA Component:'),
                                                             choices = c(names(eye_pca_data) %>% head(10)), 
                                                             multiple = FALSE, selected = "PC1"),
                                              selectizeInput('pca_component_two', strong('Second PCA Component:'),
                                                             choices = c(names(eye_pca_data) %>% head(10)), 
                                                             multiple = FALSE, selected = "PC2"),
                                              actionButton('pca_button','(Re)Draw PCA Plot!', 
                                                           style='background-color: #3399ff; color: #ffffff'),
                                       ),
                                       column(10,
                                              conditionalPanel(condition = "input.pca_button == 0", 
                                                               "Select a first and second PCA component then click the (RE)Draw PCA Plot! button. 
                                                               It may take a few seconds for the plot to appear."),
                                              conditionalPanel(condition = "input.pca_data_format == 'PCA Plot'",
                                                               plotlyOutput("eye_pca_plot", height = "1080", width = "100%"),
                                              ),
                                       ),
                                     ),
                                   ), br(), br(),
                                   fluidRow(includeHTML("www/footer.html"))
                          ),
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
                                                             choices = c("Gene 2022",
                                                                         "Gene 2019", "Transcript 2019",
                                                                        "Gene 2017", "Transcript 2017",
                                                                         "DNTx v01"),
                                                             selected = 'Gene 2022')),
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
                          # Data Download ---------------
                          tabPanel('Data Download',
                                   fluidPage(
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
                                     fluidRow(column(width = 8, offset = 1, img(src='simplified_workflow.svg', width = 300))),
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
                                                     fluidRow(column(6, img(src='sample_count_2017_2019.svg', align='middle', width = 1000))),
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
                                                     tags$a(href="https://github.com/davemcg/eyeIntegration_app/blob/master/inst/citations.md",
                                                            "here."))),
                                     
                                     fluidRow(column(width = 8, offset = 1, h2('Source Code'))),
                                     fluidRow(column(width = 8, offset = 1,
                                                     'The source code and data for this web application is available ', tags$a(href='https://gitlab.com/davemcg/Human_eyeIntegration_App', 'here.'))),
                                     fluidRow(column(width = 8, offset = 1, h2('Problems?'))),
                                     fluidRow(column(width = 8, offset = 1,
                                                     'First check the FAQ by clicking on ', strong('Information'), 'in the above header, then on ', strong('FAQs'), br(), br(), 'Other issues can be reported two ways: ',
                                                     tags$a(href = "mailto:mcgaugheyd@mail.nih.gov?Subject=eyeIntegration%20Issue", "email"), 'or ',
                                                     tags$a(href ='https://github.com/davemcg/eyeintegration_app/issues', 'GitHub Issue Tracker'))),br(), br(),
                                     fluidRow(includeHTML("www/footer.html")),
                                     br(), br()
                                   )),
                          # News ---------------
                          tabPanel('News',
                                   fluidPage(
                                     fluidRow(column(width = 8, offset = 1, h2('2022-11-03 | v2.00'))),
                                     fluidRow(column(width = 8, offset = 1, 'EVERYTHING CHANGE PANIC')),
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
                          tabPanel('FAQs', fluidPage(
                            navlistPanel("eyeIntegration FAQs",
                                         tabPanel("How was the original (2017) data generated acquired and processed?",
                                                  "See our publication:",
                                                  tags$a(href="https://academic.oup.com/hmg/advance-article/doi/10.1093/hmg/ddy239/5042913",
                                                         "https://academic.oup.com/hmg/advance-article/doi/10.1093/hmg/ddy239/5042913")),
                                         tabPanel("2017 Data Workflow",
                                                  img(src='2017_workflow.svg')),
                                         tabPanel("How was the 2019 data generated and processed?",
                                                  "We will soon have a pre-print describing the new 2019 automated data analysis workflow. The code-base for the build can be found ", tags$a(href="https://github.com/davemcg/eyeIntegration_data_build", "here.")),
                                         tabPanel("2019 Data Workflow",
                                                  img(src='2019_workflow.svg', width = 900))),
                            navlistPanel("Pan-Tissue Expression FAQs",
                                         tabPanel("How do I use the Pan-Tissue Expression section?",
                                                  "First you pick the dataset (2017 or 2019) [2]. Then you can tweak the 'Genes' [3] and 'Tissues' [4] by clicking in them and starting to type (allowed values will auto-fill). You can also delete values by clicking on them and hitting the 'delete' key on your keyboard. You can tweak the display of the box plots a bit by changing the 'Number of columns' field [5]. A higher number will squeeze more plots in each column. When you are done tweaking those parameters, click the big blue '(Re)Draw Plot!' button [6] and wait a few seconds.", br(), br(), 'If you mouse over a data point [8], you will get metadata about that particular sample.', br(), br(), img(src='pantissue_screenshot.png', width = 900), br(), br()),
                                         tabPanel("What Pan-Tissue Expression data is displayed?",
                                                  "Each gene gets its own box. The y-axis is length scaled TPM (log2 transformed). The x axis is samples, colored by tissue. The right panels [9, 10] are tables with external links to gene info [9] and the absolute TPM values and the rank of the gene in the particular tissue (lower is more highly expressed) [10].", br(), br(), img(src='pantissue_screenshot.png', width = 900), br(), br()),
                                         tabPanel("How do I interpret the 'Decile' column in the data table?",
                                                  "To give a rough sense of how highly expressed the gene is in the tissue, the decile of expression is given in [10]; 10 is the highest decile of expression and 1 is the lowest.", br(), br(), img(src='pantissue_screenshot.png', width = 900), br(), br()),
                                         tabPanel("What is that 'Fold Change' radio button?",
                                                  "This is a very simple differential expression test. When you click this radio button [1] the view changes, where the expression is being shown *relative* to the reference tissues [2]. By default it is all of the tissues; in the screenshot Retina - Organoid and Retina - Stem Cell Line are the baseline expression samples. The data table on the right [3] will then display log2 fold change, average expression, and the p-value (from a t-test) of the differential expression.", br(), br(), img(src='pantissueFC_screenshot.png', width = 900), br(), br()),
                                         tabPanel("What is the 'Heatmap' radio button?",
                                                  "This produces a 2D visualization, with each gene as a row and each tissue as a column. More yellow is more expressed. It is a efficient way to display the expression of many genes and tissues.", br(), br(), img(src='pantissue_heatmap.png', width = 900), br(), br()),
                                         tabPanel("Custom Gene / Tissue combination via the URL",
                                                  "If you use the Pan-Tissue Boxplot feature a lot, you may find it frustrating to have to input in your favorite genes and tissues. We have added the ability to use a custom url to load in the genes and tissues of your choice. Previously you had to build this link youself - but now there's a handy button [7] you can click that will re-create the parameters. One downside is that the web app is continually using the URL dataset, which makes it impossible for you to change it. You can simply reload the web page with the custom bits.", br(), br(), img(src='pantissue_screenshot.png', width = 900), br(), br()),
                                         tabPanel("What is Differential?",
                                                  "This is short for differential expression. We have pre-calculated 55+ differential expression tests. All eye tissue - origin pairs were compared to each other. We also have a synthetic human body set, made up of equal numbers of GTEx tissues (see manuscript, above, for more details). The word cloud displayed shows as many as the top 75 terms used in enriched GO terms in the selected comparison. The table data shows the actual GO terms. You can search for the comparison of your choice."),
                                         tabPanel("What are logFC, AveExp, t, P.Value, adj.Pval, and B in Differential?",
                                                  'These are the values taken from the limma differential expression topTable() summary table. The following has been taken from the limma manual and edited to match parameters we used (https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf):',br(),br(),'A number of summary statistics are presented by topTable() for the top genes and the selected contrast. The logFC column gives the value of the contrast. Usually this represents a log2-fold change between two or more experimental conditions although sometimes it represents a log2-expression level. The AveExpr column gives the average log2-expression level for that gene across all the arrays and channels in the experiment. Column t is the moderated t-statistic. Column P.Value is the associated p-value and adj.P.Value is the p-value adjusted for multiple testing (False Discovery Rate corrected).', br(), br(), 'The B-statistic (lods or B) is the log-odds that the gene is differentially expressed. Suppose for example that B = 1.5. The odds of differential expression is exp(1.5)=4.48, i.e, about four and a half to one. The probability that the gene is differentially expressed is 4.48/(1+4.48)=0.82, i.e., the probability is about 82% that this gene is differentially expressed. A B-statistic of zero corresponds to a 50-50 chance that the gene is differentially expressed. The B-statistic is automatically adjusted for multiple testing by assuming that 1% of the genes, or some other percentage specified by the user in the call to eBayes(), are expected to be differentially expressed. The p-values and B-statistics will normally rank genes in the same order. In fact, if the data contains no missing values or quality weights, then the order will be precisely the same.'),
                                         tabPanel("Mouse Single Cell Retina Expression Datasets?",
                                                  'The Macosko data is a single-cell (~45,000) retina RNA-seq mouse P14 C57BL/6 dataset  from Mackosko and McCarroll\'s field defining ', tags$a(href="https://www.ncbi.nlm.nih.gov/pubmed/26000488","paper."), 'The cluster / cell type assignments are taken from ', tags$a(href="http://mccarrolllab.com/wp-content/uploads/2015/05/retina_clusteridentities.txt","here."), 'The Clark data is a 100,000 cell plus mouse retina RNA-seq time series dataset. Their pre-publication manuscript is on ', tags$a(href="https://www.ncbi.nlm.nih.gov/pubmed/26000488","bioRxiv. "), 'Data was pulled from ', tags$a(href="https://github.com/gofflab/developing_mouse_retina_scRNASeq", 'here.'), br(), br()),
                                         tabPanel("Mouse Single Cell Retina Expression Heatmaps / Tables?",
                                                  "To efficiently display a huge amount of information, expression across many individual cells is averaged by cell type, (if available) age, and gene. You can select the Macosko or Clark dataset [1], then one gene [2] to plot. The gene expression is displayed as a heatmap, with each row being a retina cell type (derived by the respective authors) and each column [4] is a time point, arranged from youngest to oldest. More yellow is higher expression [5].", br(), br(), img(src='sc_heatmap_top.png', width = 900), br(), img(src='sc_heatmap_bottom.png', width = 900), br(), br()),
                                         tabPanel("Mouse Single Cell Retina Expression Overlay?",
                                                  "You can add the rank of expression (or rank of percentage of cells with detectable expression of selected gene) with this radio", br(), br(), img(src='sc_heatmap_overlay.png', width = 900), br(), br())),
                            navlistPanel("2D Clustering FAQ",
                                         tabPanel("What is Bulk RNA-seq by Tissue?",
                                                  "This shows the t-SNE tissue clustering for the bulk human eye tissues along with the GTEx data-set. Hovering the mouse over each data point will show the metadata. Changing the perplexity will demonstrate how low values artificially create sub-groups while higher value (above 30 or so) largely recapitulate tissue type. It also demonstrates that the tissue clustering is stable at higher perplexities.", br(), br(), img(src='bulk_tsne.png', width = 900), br(), br()),
                                         tabPanel("What is Mouse Retina Single Cell...?",
                                                  "Each data point is a single cell from the ", tags$a(href="https://www.biorxiv.org/content/early/2018/07/30/378950","Macosko and McCarroll"), 'or the ', tags$a(href="https://www.ncbi.nlm.nih.gov/pubmed/26000488","Clark and Blackshaw"), '[2]. Dimensionality reduction with done with the t-SNE (Macosko) or UMAP (Clark) algorithm. Cluster assignments were taken from the respective papers. While we did the t-SNE on the Macosko data, the Clark authors provided the UMAP coordinates. The Clark dataset was generated across multiple time-points during development and thus, you can select time points of interest [4]. Only one gene can selected at a time [3], as it is very computationally expensive to plot many points. Points (cells) expressing the gene of interest are plotted in darker color [arrow]. Hovering over each point [6] will show the metadata for the cell.', br(), br(), img(src='sc_2d.png', width = 900), br(), br())),
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
                                                  'This data table set shows the data used to make the Mouse Single Cell Retina Expression plots in Gene Expression.'))))),
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
                                                             choices = c('Gene 2019','Transcript 2019', 'Gene 2022'))),
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
                                              ggiraphOutput('tsne', height = '1200px', width = '100%')))
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