#' Visualize microbiome-metabolite network results using ComplexHeatmap
#'
#' @param result The output from analyze_microbiome_metabolite_network()
#' @param cluster_rows Logical. Whether to cluster rows (default: TRUE)
#' @param cluster_columns Logical. Whether to cluster columns (default: TRUE)
#' @param show_values Logical. Whether to show matrix values in the heatmap (default: FALSE)
#' @param value_cex Size of the text values if shown (default: 0.5)
#' @param sig_threshold_for_viz What to use for white color in case method is "corr_test" (default: 0.05)
#' @param save_pdf Either FALSE (default) or a file name (e.g., "plot.pdf") to save heatmap
#' @param width Width of heatmap in inches if saved (default: 7)
#' @param height Height of heatmap in inches if saved (default: 7)
#' @param save_path Path to save to if saved (default: script directory)
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap
#'
#' @return A ComplexHeatmap object
#' 
#' @examples
#' result = readRDS("data/corr_test.rds")
#' ht = visualize_network_heatmap(result, cluster_rows = T, cluster_columns = T,
#'                                row_names_gp=grid::gpar(fontsize=7),
#'                                column_names_gp=grid::gpar(fontsize=7))
#' ht                   
#' result = readRDS("data/corr.rds")
ht = visualize_network_heatmap(result, cluster_rows = T, cluster_columns = T,
                               row_names_gp=grid::gpar(fontsize=10),
                               column_names_gp=grid::gpar(fontsize=10),
                               save_pdf = "plot.pdf", width=12, height=12)
#' @export
visualize_network_heatmap <- function(result, 
                                      cluster_rows = TRUE,
                                      cluster_columns = TRUE,
                                      show_values = FALSE,
                                      value_cex = 0.5,
                                      sig_threshold_for_viz = NULL,
                                      save_pdf = FALSE, 
                                      width = 7, 
                                      height = 7,
                                      save_path = NULL,
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
  
  ht = ComplexHeatmap::Heatmap(mat,
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
                               width = grid::unit(width, "in"), 
                               height = grid::unit(height, "in"),
                               ...)
      
  if (is.character(save_pdf)) {
    if (is.null(save_path)) {
      if (is.null(sys.frame(1)$ofile)) {
        # if running interactively
        save_path <- getwd()
      } else {
        save_path <- dirname(sys.frame(1)$ofile)
      }
    }
    full_path <- file.path(save_path, save_pdf)
    grDevices::pdf(full_path, width=width+5, height=height+5)
    ComplexHeatmap::draw(ht)
    grDevices::dev.off()
    message(sprintf("Saved heatmap to '%s'", full_path))
    invisible(ht)
  } else {
    return(ht)
  }
}

