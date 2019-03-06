#' Download eyeIntegration data
#' 
#' Downloads the data underlying eyeIntegration.
#' 
#' @param file The file name on biowulf2
#' @param destdir The destination for the eyeIntegration tar.gz file. 
#' @param destfile The file name for the eyeIntegration tar.gz file
#' @param method Download method
#' @param delete_tar Delete the tar file after data extraction?
#' @param verbose Turn on verbosity for untar()
#' 
#' @examples
#' \dontrun{get_eyeIntegration_datasets()} 
#' @export
get_eyeIntegration_datasets <-
  function (file = "eyeIntegration_v102_03_heavyweight.tar.gz",
            destdir = system.file('app', package = "eyeIntegrationApp"), 
            destfile = "eyeIntegration_data.tar.gz",
            method,
            delete_tar = TRUE,
            verbose = TRUE)
  {
    # try to make app destination directory
    if (!dir.exists(destdir)){
      dir.create(destdir)
    }
    if (missing(method))
      method <- ifelse(!is.null(getOption("download.file.method")),
                       getOption("download.file.method"), "auto")
    localfile <- file.path(destdir, destfile)   
    options(warn=-1)
    
    url = paste0('https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/', file)
    message("Huge download starting!\n")
    utils::download.file(url, destfile = localfile, mode = "wb", method=method)
    message("Decompressing data...\n")
    utils::untar(localfile, exdir = destdir, verbose = verbose)
    # copy server.R and ui.R to custom location, if used
    if ( destdir != system.file('app', package = "eyeIntegrationApp")) {
      file.copy(from = system.file('app/server.R', package = "eyeIntegrationApp"), 
                to = destdir)
      file.copy(from = system.file('app/ui.R', package = "eyeIntegrationApp"), 
                to = destdir)
    }
    message("Done!\n")
    if (isTRUE(delete_tar)){
      file.remove(localfile)
    }
  }
