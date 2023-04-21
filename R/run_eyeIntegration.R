#' Run eyeIntegration App
#' 
#' Wrapper that starts the eyeIntegration Shiny App. You must run 
#' \code{\link{get_eyeIntegration_datasets}} first. 
#' 
#' @param app_path The path to the folder holding the server.R, ui.R, and www/ 
#' folder.
#' @examples
#' \dontrun{run_eyeIntegration()} 
#' @export
run_eyeIntegration <- function(app_path = system.file('app', package = "eyeIntegrationApp")) { shiny::runApp(app_path) }
