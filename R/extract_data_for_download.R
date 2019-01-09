library(tidyverse)
library(pool)
library(RSQLite)
gene_pool_2017 <- dbPool(drv = SQLite(), dbname = "./www/2017/eyeIntegration_human_2017_01.sqlite", idleTimeout = 3600000)
gene_pool_2019 <- dbPool(drv = SQLite(), dbname = "./www/2019/eyeIntegration_human_expression_2019_01.sqlite", idleTimeout = 3600000)

metadata_2017 <- gene_pool_2017 %>% tbl('metadata') %>% as.tibble() %>% rowwise() %>% 
  mutate(sample_accession = gsub('E-MTAB-','E.MTAB.',sample_accession),
         Sub_Tissue = gsub('_',' - ', Sub_Tissue)) 
write_tsv(metadata_2017, path = '~/Desktop/2017_metadata.tsv.gz')
gene_pool_2017 %>% tbl('lsTPM_gene') %>% as.tibble() %>% 
  filter(sample_accession %in% metadata_2017$sample_accession) %>% spread(sample_accession, value) %>% 
  write_tsv(., path = '~/Desktop/2017_gene_TPM.tsv.gz')
gene_pool_2017 %>% tbl('lsTPM_TX') %>% as.tibble() %>% 
  filter(sample_accession %in% metadata_2017$sample_accession) %>% spread(sample_accession, value) %>% 
  write_tsv(., path = '~/Desktop/2017_tx_TPM.tsv.gz')

metadata_2019 <- gene_pool_2019 %>% tbl('metadata') %>% as.tibble() %>% rowwise() %>% 
  mutate(sample_accession = gsub('E-MTAB-','E.MTAB.',sample_accession),
         Sub_Tissue = gsub('_',' - ', Sub_Tissue)) 
write_tsv(metadata_2019, path = '~/Desktop/2019_metadata.tsv.gz')
gene_pool_2019 %>% tbl('lsTPM_gene') %>% as.tibble() %>% 
  filter(sample_accession %in% metadata_2019$sample_accession) %>% spread(sample_accession, value) %>% 
  write_tsv(., path = '~/Desktop/2019_gene_TPM.tsv.gz')
gene_pool_2019 %>% tbl('lsTPM_TX') %>% as.tibble() %>% 
  filter(sample_accession %in% metadata_2019$sample_accession) %>% spread(sample_accession, value) %>% 
  write_tsv(., path = '~/Desktop/2019_tx_TPM.tsv.gz')
