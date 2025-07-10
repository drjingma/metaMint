#' Visualize microbiome-metabolite network results using ComplexHeatmap
#'
#' @param result The output from analyze_microbiome_metabolite_network()
#' @param cluster_rows Logical. Whether to cluster rows (default: TRUE)
#' @param cluster_columns Logical. Whether to cluster columns (default: TRUE)
#' @param show_values Logical. Whether to show matrix values in the heatmap (default: FALSE)
#' @param value_cex Size of the text values if shown (default: 0.5)
#' @param sig_threshold_for_viz What to use for white color in case method is "corr_test" (default: 0.05)
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap
#'
#' @return A ComplexHeatmap object
#' 
#' @examples
#' result = readRDS("data/bisbm_test.rds")
#' visualize_network_heatmap(result, cluster_rows = T, cluster_columns = T,
#'                           row_names_gp=grid::gpar(fontsize=7),
#'                           column_names_gp=grid::gpar(fontsize=7))
#'                           
#' @export
visualize_network_heatmap <- function(result, 
                                      cluster_rows = TRUE,
                                      cluster_columns = TRUE,
                                      show_values = FALSE,
                                      value_cex = 0.5,
                                      sig_threshold_for_viz = NULL, 
                                      ...) {
  method <- strsplit(names(result)[1], "\\.")[[1]][2]
  
  if (method != "corr_test") {
    if (!is.null(sig_threshold_for_viz)) {
      warning("sig_threshold_for_viz is only used when method = 'corr_test'; the provided value will be ignored.")
    }
  } else {
    if (is.null(sig_threshold_for_viz)) {
      sig_threshold_for_viz = 0.05
      message(sprintf("`sig_threshold_for_viz` not provided, using default: '%s'.", sig_threshold_for_viz))
    }
  }
  
  
  if (method %in% c("corr", "pcorr")) {
    mat <- as.matrix(result[[1]])
    title <- paste("Heatmap of", ifelse(method == "corr", "Correlation", "Partial Correlation"))
    col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  } else if (method == "corr_test") {
    mat <- as.matrix(result[[1]])
    title <- "Heatmap of FDR-adjusted p-values"
    col_fun <- circlize::colorRamp2(c(0, 0.05, 1), c("red", "white", "blue"))
  }
  
  ComplexHeatmap::Heatmap(mat,
                          name = switch(method,
                                        "corr" = "Corr",
                                        "pcorr" = "Partial\nCorr",
                                        "corr_test" = "FDR"),
                          col = col_fun,
                          cluster_rows = cluster_rows,
                          cluster_columns = cluster_columns,
                          show_row_names = TRUE,
                          show_column_names = TRUE,
                          cell_fun = if (show_values) {
                            function(j, i, x, y, width, height, fill) {
                              grid::grid.text(sprintf("%.3f", mat[i, j]), x, y, 
                                              gp = grid::gpar(fontsize = 10 * value_cex))
                            }
                          } else NULL,
                          column_title = title,
                          ...)
}
