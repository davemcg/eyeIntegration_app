
#' normalize_data
#'
#' Applies *sample* level normalization
#' 
#' This function is used by both metamoRph::metamoRph and metamoRph::run_pca
#'
#' @param feature_by_sample features (row) by sample (column) matrix 
#' @param log1p Default is TRUE
#' @return  a list object containing the prcomp "center" and "scale"
#' values for use in the metamoRph function
#' @export
normalize_data <- function(feature_by_sample, sample_scale = 'cpm', log1p = TRUE){
  new_scale <- feature_by_sample # in case no applicable sample scaling is selected
  
  
  message("Sample CPM scaling")
  new_scale <-  ((feature_by_sample) / 
                   colSums(feature_by_sample)) * 1e6
  
  if (log1p){
    message("log1p scaling")
    new_scale <- log1p(new_scale)
  } else {
    message("no log1p scaling")
  }
  
  
  return(new_scale)
}



#' metamoRph
#'
#' Takes in a count matrix (where genes (features) are rows and samples are
#' columns) as well as a named vector with the eigenvalues (see [metamoRph::run_pca()])
#' and pulls the gene (feature) information from the rotation vector and cuts down
#' the new_counts matrix to match the rotation vector gene (feature) names. Any
#' genes (features) missing from the input new_counts matrix will be replaced with zeros.
#'
#' The function will scale the new_counts matrix in the same manner as [metamoRph::run_pca()]
#' and matrix multiply by the rotation vector. The output is equivalent
#' to the prcomp "$x" matrix.
#'
#' @param new_counts raw gene count matrix (where genes are rows and samples are
#' columns)
#' @param rotation matrix where the row names are genes and the col names are the
#' principal components. If you used metamoRph::run_pca() then this would be
#' in the output$PCA$rotation slot.
#' @param center_scale list object where the $center slot has the center values
#' and the $scale slot has the scale value for the "scale" function. If you do not
#' give a value here, then feature/gene scaling WILL NOT HAPPEN.
#' @param sample_scale Options include `cpm`, `seurat`, `zscale`,  and `none`. This MUST match
#' the normalization used in the data for `run_pca` (or where your rotation came from). 
#' @param log1p Default is TRUE. Again, must match the normalization used in the data for 
#' `run_pca` (or where your rotation came from).
#' @param normalization Default is TRUE, if set to FALSE will override `sample_scale` and `log1p` 
#' and not do any sample scaling
#' @param feature_scale Default is FALSE, if TRUE will apply feature (gene) scaling to the input (query)
#' data (if center_scale is left empty). VERY DANGEROUS as this will likely put the query data in a different scale than your
#' reference data. Use with deliberate intent.
#' @return A matrix with the transformated eigenvalue matrix which should be equivalent
#' to the original rotation matrix's eigenvalue/pattern matrix (The $x slot from
#' the output of prcomp)
#' @import dplyr
#' @import tidyr
#' @importFrom Matrix t colSums
#' @importFrom tibble enframe
#' @importFrom stats prcomp
#' @export
metamoRph <- function(new_counts,
                      rotation,
                      center_scale = NULL,
                      normalization = TRUE,
                      feature_scale = FALSE,
                      log1p = TRUE){
  value <- ensgene <- name <- NULL # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  
  if (is.null(center_scale)){
    warning(
      "No center and scale values provided. Projection may be erroneous 
(unless you already ensured your input data
is scaled the same as the rotation/PCA data input data")
  }
  if (normalization){
    new_scale <- normalize_data(new_counts, 
                                log1p = log1p)
  } else {
    message("Normalization skipped")
    new_scale <- new_counts 
  }
  
  # extract gene names from new data and make upper case
  row_genes <- row.names(new_scale) %>% toupper()
  row.names(new_scale) <- row_genes
  # this bit is for the mcgaughey (sc)EiaD data which uses a gene naming
  # scheme that pastes together the common gene name with the ENSEMBL id
  suppressWarnings(feature_id_table <- row.names(rotation) %>%
                     toupper() %>%
                     enframe() %>%
                     separate(value, c("feature_id", "ensgene"), sep = " \\(") %>%
                     mutate(ensgene = gsub(")", "", ensgene)) %>% dplyr::select(-name))
  overlap_with_ID <- row_genes[row_genes %in% feature_id_table$feature_id]
  overlap_with_ens <- row_genes[row_genes %in% feature_id_table$ensgene]
  
  if ((length(overlap_with_ID) < 2) &
      (length(overlap_with_ens) < 2)){
    stop("Failure to align gene (feature) names, please check your
    input matrix rownames against the names of your rotation vector")
  } else {
    # select column ID type to use
    if (length(overlap_with_ID) >
        length(overlap_with_ens)){
      column_val <- 1
    } else {
      column_val <- 2
    }
  }
  
  # Only genes that match the rotation/eigenvalues genes used
  feature_universe <- feature_id_table %>% pull(column_val) %>% make.unique()
  feature_universe[feature_universe == ''] = 'X'
  cutdown <- new_scale[feature_universe[feature_universe %in% row.names(new_scale)], , drop = FALSE]
  message(paste0("Aligned ", max(nrow(cutdown)), " features/genes."))
  # add in missing genes
  # by building a matrix composed of the missing genes and the samples
  # then rbind it with the existing data
  not_included <- feature_universe[!feature_universe %in% intersect(row.names(cutdown), feature_universe)]
  data <- matrix(nrow=length(not_included), ncol=ncol(cutdown))
  row.names(data) <- not_included
  colnames(data) <- colnames(cutdown)
  cutdown <- rbind(cutdown, data)
  
  # Match order of feature_ids from reference PCA
  ## ensure the input data is in the
  ## same order as the rotation matrix
  cutdown <- cutdown[feature_universe, , drop = FALSE]
  
  if (!missing(center_scale)){
    message("Applying feature scaling from reference data")
    scaled_cutdown <- scale(Matrix::t(cutdown), center_scale$center,
                            center_scale$scale)
  } else if (feature_scale) {
    message("Running feature scaling on query data")
    scaled_cutdown <- Matrix::t(scale(cutdown))
  } else {
    scaled_cutdown <- Matrix::t(cutdown)
  }
  # Replace NAs with 0
  scaled_cutdown[is.na(scaled_cutdown)] <- 0
  
  # project new data onto the PCA space
  ## matrix multiply the expression of the scaled new data against
  ## the supplied eigenvector
  projected_PC_space <- scaled_cutdown %*% rotation %>% as.data.frame()
  
  return((projected_PC_space))
}