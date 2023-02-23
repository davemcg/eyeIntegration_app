#' Build bar plot graphic
#' 
#' Builds bar plot on main page 
#' 
#' @examples
#' \dontrun{source('R/sample_count_2017_2019_2022.R')}
#' 
library(tidyverse)
library(pool)
library(RSQLite)
library(colorspace)

#sample_count_2022 <- function(){
  app_location <- '/Users/parikhpp/git/eyeIntegration_app_v2_pp/inst/app'
  gene_pool_2022 <- dbPool(drv = SQLite(), dbname = paste0(app_location, "/www/2022/eyeIntegration_2022_human.sqlite"))
  core_tight_2022 <- gene_pool_2022 %>% tbl('metadata') %>% as_tibble()
  core_tight_2022 <- core_tight_2022 %>% mutate(Sub_Tissue = case_when(is.na(Sub_Tissue) ~ '', TRUE ~ Sub_Tissue), 
                                                Source = case_when(is.na(Source) ~ '', TRUE ~ Source), 
                                                Age = case_when(is.na(Age) ~ '', TRUE ~ Age),
                                                Perturbation = case_when(is.na(Perturbation) ~ '', TRUE ~ Perturbation),
                                                Tissue = case_when(is.na(Tissue) ~ '', TRUE ~ Tissue))
  
  # fix tissue <-> color
  meta <- 'core_tight_2022'
  tissue_col <- scale_fill_manual(values = setNames(c(pals::glasbey(n = 32), 
                                                      pals::kelly(n = get(meta) %>% pull(Tissue) %>% unique() %>% length() - 32 + 1)[-1]) %>% 
                                                      colorspace::lighten(0.3), get(meta) %>% pull(Tissue) %>% unique() %>% sort()))
  
  core_tight_2022 %>% 
    arrange(Tissue) %>% 
    mutate(GTEx = case_when(study_accession == 'SRP012682' ~ 'GTEx', TRUE ~ 'Eye')) %>% 
    select(run_accession, Tissue, Sub_Tissue, Age, Source, Perturbation, GTEx) %>% unique() %>%
    mutate(expanded_name = paste0(Tissue, sep = " | ", Sub_Tissue, sep = " | ", 
                                  Age, sep = " | ", Perturbation, sep = " | ", 
                                  Source)) %>%
    group_by(Tissue, Sub_Tissue, Age, Source, Perturbation, expanded_name, GTEx) %>% 
    count(name="Count") %>% 
    ungroup() %>% 
    ggplot(aes(x= expanded_name, y= Count, fill=Tissue)) + geom_bar(stat = 'identity', position = 'dodge') + 
    theme_minimal() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.3)) + 
    xlab('') + 
    facet_grid(~GTEx, scales = 'free_x', space = 'free_x')+
    theme(text = element_text(family = 'Arial', size = 10)) +
    scale_fill_discrete_sequential(palette = 'viridis')  +
    tissue_col
    
  ggsave(filename = paste0(app_location, '/www/sample_count_2022.svg'),dpi = 'retina', height = 10, width=15)
#}
