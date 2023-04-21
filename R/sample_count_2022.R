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
library(svglite)

sample_count_2022 <- function(){
  app_location <- '/Users/mcgaugheyd//git/eyeIntegration_app/inst/app'
  #gene_pool_2017 <- dbPool(drv = SQLite(), dbname = paste0(app_location, "/www/2017/eyeIntegration_human_2017_01.sqlite"), idleTimeout = 3600000)
  #gene_pool_2019 <- dbPool(drv = SQLite(), dbname = paste0(app_location, "/www/2019/EiaD_human_expression_2019_04.sqlite"), idleTimeout = 3600000)
  gene_pool_2022 <- dbPool(drv = SQLite(), dbname = paste0(app_location, "/www/2022/eyeIntegration_2022_human.sqlite"))
  
  core_tight_2017 <- gene_pool_2017 %>% tbl('metadata') %>% as_tibble()
  core_tight_2019 <- gene_pool_2019 %>% tbl('metadata') %>% as_tibble()
  core_tight_2022 <- gene_pool_2022 %>% tbl('metadata') %>% as_tibble()
  core_tight_2022 <- core_tight_2022 %>% mutate(Sub_Tissue = case_when(is.na(Sub_Tissue) ~ '', TRUE ~ Sub_Tissue), 
                                                Source = case_when(is.na(Source) ~ '', TRUE ~ Source), 
                                                Age = case_when(is.na(Age) ~ '', TRUE ~ Age),
                                                Perturbation = case_when(is.na(Perturbation) ~ '', TRUE ~ Perturbation),
                                                Tissue = case_when(is.na(Tissue) ~ '', TRUE ~ Tissue))
  
  # # fix tissue <-> color
  # meta <- 'core_tight_2022'
  # tissue_col <- scale_fill_manual(values = setNames(c(pals::glasbey(n = 32), 
  #                                                     pals::kelly(n = get(meta) %>% pull(Tissue) %>% unique() %>% length() - 32 + 1)[-1]) %>% 
  #                                                     colorspace::lighten(0.3), get(meta) %>% pull(Tissue) %>% unique() %>% sort()))
  
  # Use global tissue values to match remainder of eyeIntegration app
  tissues <- c(core_tight_2017$Tissue, core_tight_2019$Tissue, core_tight_2022$Tissue)%>% unique() %>% sort() 
  tissue_fill <- scale_fill_manual(values = setNames(c(pals::polychrome()[3:36],  pals::kelly()[c(3:7,10:21)])[1:length(tissues)], tissues %>% sort()))
  
  a <- core_tight_2022 %>% 
    arrange(Tissue) %>% 
    mutate(GTEx = case_when(study_accession == 'SRP012682' ~ 'GTEx', TRUE ~ 'Eye')) %>% 
    select(run_accession, Tissue, Sub_Tissue, Age, Source, Perturbation, GTEx) %>% unique() %>%
    mutate(Sub_Tissue = case_when(is.na(Sub_Tissue) ~ '', TRUE ~ Sub_Tissue), 
           Source = case_when(is.na(Source) ~ '', TRUE ~ Source), 
           Age = case_when(is.na(Age) ~ '', TRUE ~ Age),
           Perturbation = case_when(is.na(Perturbation) ~ '', TRUE ~ Perturbation)) %>% 
    mutate(Sub_Tissue = glue::glue("<span style='color:#000000FF'>{Sub_Tissue}</span>"),
           Source = glue::glue("<span style='color:#1E46A2FF'>{Source}</span>"),
           Age = glue::glue("<span style='color:#FB323BFF'>{Age}</span>"),
           Perturbation = glue::glue("<span style='color:#85660D'>{Perturbation}</span>")
    ) %>% 
    mutate(expanded_name = paste0(Tissue, sep = " | ", Sub_Tissue, sep = " | ", 
                                  Age, sep = " | ", Perturbation, sep = " | ", 
                                  Source)) %>%
    group_by(GTEx, Tissue, Sub_Tissue, Age, Source, Perturbation, expanded_name) %>% 
    filter(GTEx == 'Eye') %>% 
    count(name="Count") %>% 
    ungroup() %>% 
    #mutate(Perturbation = case_when(grepl('MGS', Source_details) ~ Source_details)) %>% 
    ggplot(data=.,aes(x=interaction(Source, Sub_Tissue, Age, Perturbation, sep = ' | '),y=Count, 
                      fill = Tissue)) +
    #geom_violin(alpha=0.5, scale = 'width') +
    geom_bar(stat = 'identity', position = 'dodge') + 
    cowplot::theme_cowplot(font_size = 15) + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.2)) +
    ylab("Count") +
    theme(strip.background = element_rect(fill = 'black'),
          strip.text = element_text(color = 'white'),
          panel.background = element_rect(fill = 'gray90'),
          plot.margin=grid::unit(c(0,0,0,0.1), "cm"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size= unit(0.2, "cm"),
          legend.spacing = unit(0.2, "cm"))  +
    tissue_fill + 
    coord_flip() + 
    facet_grid(rows = vars(Tissue),  scales = 'free_y', space = 'free') +
    theme(strip.text.y.right = element_text(angle = 0)) +
    theme(
      axis.text.y = element_markdown(),
      axis.title.y = element_markdown()) +
    labs(x = "<span style='color:#1E46A2FF'>Source</span> | 
       <span style='color:#000000FF'>Sub Tissue</span> |
       <span style='color:#FB323BFF'>Age</span> |
       <span style='color:#85660D'>Perturbation</span>") +
    theme(legend.position = "none")
  
  b <- core_tight_2022 %>% 
    arrange(Tissue) %>% 
    mutate(GTEx = case_when(study_accession == 'SRP012682' ~ 'GTEx', TRUE ~ 'Eye')) %>% 
    select(run_accession, Tissue, Sub_Tissue, Age, Source, Perturbation, GTEx) %>% unique() %>%
    mutate(Sub_Tissue = case_when(is.na(Sub_Tissue) ~ '', TRUE ~ Sub_Tissue), 
           Source = case_when(is.na(Source) ~ '', TRUE ~ Source), 
           Age = case_when(is.na(Age) ~ '', TRUE ~ Age),
           Perturbation = case_when(is.na(Perturbation) ~ '', TRUE ~ Perturbation)) %>% 
    mutate(Sub_Tissue = glue::glue("<span style='color:#000000FF'>{Sub_Tissue}</span>"),
           Source = glue::glue("<span style='color:#1E46A2FF'>{Source}</span>"),
           Age = glue::glue("<span style='color:#FB323BFF'>{Age}</span>"),
           Perturbation = glue::glue("<span style='color:#85660D'>{Perturbation}</span>")
    ) %>% 
    mutate(expanded_name = paste0(Tissue, sep = " | ", Sub_Tissue, sep = " | ", 
                                  Age, sep = " | ", Perturbation, sep = " | ", 
                                  Source)) %>%
    group_by(GTEx, Tissue, Sub_Tissue, Age, Source, Perturbation, expanded_name) %>% 
    filter(GTEx == 'GTEx') %>% 
    count(name="Count") %>% 
    ungroup() %>% 
    #mutate(Perturbation = case_when(grepl('MGS', Source_details) ~ Source_details)) %>% 
    ggplot(data=.,aes(x=interaction(Source, Sub_Tissue, Age, Perturbation, sep = ' | '),y=Count, 
                      fill = Tissue)) +
    #geom_violin(alpha=0.5, scale = 'width') +
    geom_bar(stat = 'identity', position = 'dodge') + 
    cowplot::theme_cowplot(font_size = 15) + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.2)) +
    ylab("Count") +
    theme(strip.background = element_rect(fill = 'black'),
          strip.text = element_text(color = 'white'),
          panel.background = element_rect(fill = 'gray90'),
          plot.margin=grid::unit(c(0,0,0,0.1), "cm"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size= unit(0.2, "cm"),
          legend.spacing = unit(0.2, "cm"))  +
    tissue_fill + 
    coord_flip() + 
    facet_grid(rows = vars(Tissue),  scales = 'free_y', space = 'free') +
    theme(strip.text.y.right = element_text(angle = 0)) +
    theme(
      axis.text.y = element_markdown(),
      axis.title.y = element_markdown()) +
    labs(x = "<span style='color:#1E46A2FF'>Source</span> | 
       <span style='color:#000000FF'>Sub Tissue</span> |
       <span style='color:#FB323BFF'>Age</span> |
       <span style='color:#85660D'>Perturbation</span>") +
    theme(legend.position = "none") + xlab('')
  
  cowplot::plot_grid(plotlist = list(a,b),ncol =2)
  ggsave(filename = paste0(app_location, '/www/sample_count_2022.svg'),dpi = 'retina', height = 15, width=15)
}
