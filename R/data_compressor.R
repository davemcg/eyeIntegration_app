#' Data Compressor
#' 
#' When given www/ location, tar gzips it for upload to biowulf2
#' 
#' @param www_location The path which contains the www/ App folder
#' @param tar_name path/name_of_tar_file.tar.gz
#' @param include_DB Should the SQLite databases be included?
#' 
#' @examples
#' \dontrun{tar_maker(www_location = paste0(system.file('app', 
#' package = 'eyeIntegrationApp')), tar_name = 'eyeIntegration_v100_01.tar.gz') }

data_compressor <- function(www_location, tar_name = '~/Desktop/eyeIntegration_v100_02.tar.gz', include_DB = T){
  setwd(www_location)
  if (include_DB){
    system(paste('tar --exclude=".*" -cvf - www/ | pigz >', tar_name))
  } else {
    system(paste('tar --exclude=".*" --exclude="*sqlite" -cvf - www/ | pigz >', tar_name)) 
  }
}
