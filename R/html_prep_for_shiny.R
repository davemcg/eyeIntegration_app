#' Trim down html
#' 
#' Removes html outside <body> and </body>
#' 
#' @param html_in Input html file
#' @param html_out Output html file path
#' 
#' @examples
#' \dontrun{get_eyeIntegration_datasets()} 
#' @import xml2
#' @import rvest

# takes knit html file as input and strips out html outside <head>
args <- commandArgs(trailingOnly = TRUE)
html_prep_for_shiny <- function(html_in, html_out){
  write_html(html_node(read_html(html_in), "body"), 
             file = html_out)
}
