library(tidyverse)
#https://eyeintegration.nei.nih.gov -> Data -> Data Download for links
metadata <- read_tsv('https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2019_metadata_04.tsv.gz')
gene_tpm <- read_tsv('https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2019_gene_TPM_04.tsv.gz')
metadata %>% sample_n(10)
gene_tpm[1:10,c(1,100:110)]

gene_tpm_long <- gene_tpm %>% 
  pivot_longer(cols=2:ncol(gene_tpm), names_to = 'sample_accession', values_to = 'TPM')
gene_tpm_long %>% sample_n(10)

gene_tpm_long <- left_join(gene_tpm_long, 
                           metadata %>% select(sample_accession, Tissue, Sub_Tissue), 
                           by = 'sample_accession')
gene_tpm_long %>% sample_n(10)


gene_tpm_long$Sub_Tissue %>% table()

gene_tpm_long %>% 
  filter(ID %in% c('ADCY10', 'RAPGEF3','RAPGEF4', 'RPE65'), 
         Sub_Tissue %in% c('ESC', 'RPE - Fetal Tissue','RPE - Stem Cell Line', 'RPE - Adult Tissue','ESC - Stem Cell Line')) %>% 
  ggplot(aes(x=Sub_Tissue,y=log2(TPM+1), color = Tissue)) + 
  geom_boxplot(width = 0.5) +
  ggbeeswarm::geom_quasirandom(alpha = 0.3) +
  facet_wrap(~ID) +
  coord_flip() + 
  cowplot::theme_cowplot() +
  xlab('')

gene_tpm_long %>% 
  filter(ID %in% c('ADCY10', 'RAPGEF3','RAPGEF4', 'RPE65'), 
         Sub_Tissue %in% c('ESC', 'RPE - Fetal Tissue','RPE - Stem Cell Line', 'RPE - Adult Tissue','ESC - Stem Cell Line')) %>% 
  write_csv(., file = '~/Desktop/rpe_adcy10_rapgef_rpe65.csv')
