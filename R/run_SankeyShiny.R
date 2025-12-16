#' Run the Sankey Shiny App
#'
#' Launches the interactive Shiny application for Sankey combination plots.
#' @return A shiny app object
#' @export
#' @import shiny
#' @import ggplot2
#' @import dplyr
#' @importFrom ggsankey geom_sankey
run_SankeyShiny <- function() {
  # 寻找包内部 inst/app 目录下的应用
  app_dir <- system.file("app/SankeyShiny.R", package = "LittleTools")
  
  if (app_dir == "") {
    stop("Could not find example directory. Try re-installing `SankeyShiny`.", call. = FALSE)
  }
  
  # 启动应用
  shiny::runApp(app_dir, display.mode = "normal")
}




