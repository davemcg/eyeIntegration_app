#' Run eyeIntegration App
#' 
#' Wrapper that starts the eyeIntegration Shiny App. You must run 
#' \code{\link{get_eyeIntegration_datasets}} first. 
#' 
#' @examples
#' \dontrun{run_eyeIntegration()} 
#' @export
run_eyeIntegration <- function() { shiny::runApp('./inst/app/') }
