#' Table Maker
#' 
#' Extracts tables for user download to put on biowulf2
#' 
#' @param app_location The path to the working App (so script can find sqlite 
#' databases)
#' @param save_dir Where to save the tsv files
#' @param version_append What to append to file 
#' 
#' @examples
#' \dontrun{extract_data_for_download(app_location = system.file('app', 
#' package = 'eyeIntegrationApp'), version_append = 01) }
#' @importFrom magrittr "%>%"

extract_data_for_download <- function(app_location, 
         save_dir = '~/Desktop',
         version_append){
  gene_pool_2017 <- pool::dbPool(drv = RSQLite::SQLite(), dbname = paste0(app_location, "/www/2017/*sqlite") %>% Sys.glob(), idleTimeout = 3600000)
  gene_pool_2019 <- pool::dbPool(drv = RSQLite::SQLite(), dbname = paste0(app_location, "/www/2019/*sqlite") %>% Sys.glob(), idleTimeout = 3600000)
  
  metadata_2017 <- gene_pool_2017 %>% dplyr::tbl('metadata') %>% tibble::as_tibble() %>% dplyr::rowwise() %>% 
    dplyr::mutate(sample_accession = gsub('E-MTAB-','E.MTAB.',sample_accession),
           Sub_Tissue = gsub('_',' - ', Sub_Tissue)) 
  readr::write_tsv(metadata_2017, path = paste0(save_dir, '2017_metadata_', version_append, '.tsv.gz'))
  gene_pool_2017 %>% dplyr::tbl('lsTPM_gene') %>% tibble::as_tibble() %>% 
    dplyr::filter(sample_accession %in% metadata_2017$sample_accession) %>% tidyr::spread(sample_accession, value) %>% 
    readr::write_tsv(., path = paste0(save_dir, '2017_gene_TPM_', version_append, '.tsv.gz'))
  gene_pool_2017 %>% dplyr::tbl('lsTPM_TX') %>% tibble::as_tibble() %>% 
    dplyr::filter(sample_accession %in% metadata_2017$sample_accession) %>% tidyr::spread(sample_accession, value) %>% 
    readr::write_tsv(., path = paste0(save_dir, '2017_tx_TPM_', version_append, '.tsv.gz'))
  
  metadata_2019 <- gene_pool_2019 %>% dplyr::tbl('metadata') %>% tibble::as_tibble() %>% dplyr::rowwise() %>% 
    dplyr::mutate(sample_accession = gsub('E-MTAB-','E.MTAB.',sample_accession),
           Sub_Tissue = gsub('_',' - ', Sub_Tissue)) 
  readr::write_tsv(metadata_2019, path = paste0(save_dir, '2019_metadata_', version_append, '.tsv.gz'))
  gene_pool_2019 %>% dplyr::tbl('lsTPM_gene') %>% tibble::as_tibble() %>% 
    dplyr::filter(sample_accession %in% metadata_2019$sample_accession) %>% tidyr::spread(sample_accession, value) %>% 
    readr::write_tsv(., path = paste0(save_dir, '2019_gene_TPM_', version_append, '.tsv.gz'))
  gene_pool_2019 %>% dplyr::tbl('lsTPM_TX') %>% tibble::as_tibble() %>% 
    dplyr::filter(sample_accession %in% metadata_2019$sample_accession) %>% tidyr::spread(sample_accession, value) %>% 
    readr::write_tsv(., path = paste0(save_dir, '2019_tx_TPM_', version_append, '.tsv.gz'))
}
