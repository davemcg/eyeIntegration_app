#' Build bar plot graphic
#' 
#' Builds bar plot on main page 
#' 
#' @examples
#' \dontrun{source('R/sample_count_2017_2019.R')} 
#' 
library(tidyverse)
library(pool)
library(RSQLite)
library(colorspace)

sample_count_2017_2019 <- function(){
  app_location <- '/Volumes/ARC168/eyeIntegration_app'
  gene_pool_2019 <- dbPool(drv = SQLite(), dbname = paste0(app_location, "/www/2019/EiaD_human_expression_2019_04.sqlite"))
  gene_pool_2017 <- dbPool(drv = SQLite(), dbname = paste0(app_location, "/www/2017/eyeIntegration_human_2017_01.sqlite"))
  core_tight_2017 <- gene_pool_2017 %>% tbl('metadata') %>% as_tibble()
  core_tight_2019 <- gene_pool_2019 %>% tbl('metadata') %>% as_tibble()
  core_both <- bind_rows(core_tight_2019 %>% mutate(Version = '2019 EiaD') %>% filter(Kept == 'Kept', Tissue != 'Choroid Plexus'), core_tight_2017 %>% mutate(Version = '2017') %>% mutate(Sub_Tissue = gsub('_',' - ', Sub_Tissue))) %>% unique()
  
  
  # fix tissue <-> color
  meta <- 'core_tight_2019'
  tissue_col <- scale_fill_manual(values = setNames(c(pals::glasbey(n = 32), pals::kelly(n = get(meta) %>% pull(Tissue) %>% unique() %>% length() - 32)) %>% colorspace::lighten(0.3), get(meta) %>% pull(Tissue) %>% unique() %>% sort()))
  
  core_both %>% mutate(GTEx = case_when(study_accession == 'SRP012682' ~ 'GTEx', TRUE ~ 'Eye')) %>% 
    select(sample_accession, Tissue, Sub_Tissue, GTEx, Version) %>% unique() %>% 
    #mutate(Sub_Tissue = case_when(sample_accession == 'SRS1955479' ~ 'Retina - Adult Tissue', TRUE ~ Sub_Tissue)) %>%
    group_by(Tissue, Sub_Tissue, Version, GTEx) %>% summarise(Count=n()) %>%
    ungroup() %>% 
    mutate(Dataset = Version, Sub_Tissue = trimws(Sub_Tissue)) %>% 
    ggplot(aes(x=Sub_Tissue, y=Count, fill=Tissue, alpha = Dataset)) + geom_bar(stat = 'identity', position = 'dodge') + theme_minimal() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.3)) + scale_colour_manual(values = c('grey','black')) + xlab('') + scale_alpha_manual(values = c(0.5,1)) + 
    facet_grid(~GTEx, scales = 'free_x', space = 'free_x')+
    theme(text = element_text(family = 'Lucida Console', size = 10)) +
    scale_fill_discrete_sequential(palette = 'viridis')  +
    tissue_col
  
  
  ggsave(filename = paste0(app_location, '/www/sample_count_2017_2019.svg'),dpi = 'retina', height = 9, width=12)
}

