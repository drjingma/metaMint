#' Visualize microbiome-metabolite network results using ComplexHeatmap
#'
#' This function creates a heatmap of microbiome-metabolite associations using
#' results from \code{\link{analyze_microbiome_metabolite_network}}. It supports
#' correlation, partial correlation, and thresholded network methods (e.g., FDR-based)
#' and uses the \pkg{ComplexHeatmap} package for rich annotation support.
#' 
#' Depending on the input method, different defaults and color scales are used:
#' \itemize{
#'   \item \code{"corr"}: visualizes correlation matrix with values from -1 to 1.
#'   \item \code{"pcorr"}: visualizes partial correlation matrix from -1 to 1.
#'   \item \code{"corr_test"}: binary adjacency matrix thresholded by FDR.
#' }
#' 
#' @param result The result object from \code{analyze_microbiome_metabolite_network()}.
#' @param adjacency_threshold Numeric; the threshold used to determine edges in the network.
#' Interpretation depends on the method:
#' \itemize{
#'   \item For \code{"corr"} or \code{"pcorr"}: values with absolute correlation below this threshold are set to 0 (no edge).
#'   \item For \code{"corr_test"}: p-values below this threshold are considered significant (i.e., edge present).
#' }
#' Default is \code{1e-6}.
#' @param show_values Logical; whether to overlay numeric values on heatmap cells (default: \code{FALSE}).
#' @param value_cex Numeric; scaling factor for value label size if \code{show_values = TRUE} (default: \code{0.5}).
#' @param save_pdf Logical or character; if \code{FALSE} (default), shows plot interactively. If a file name (e.g., \code{"plot.pdf"}), saves heatmap to file.
#' @param width Numeric; width of heatmap in inches when saving to file (default: \code{7}).
#' @param height Numeric; height of heatmap in inches when saving to file (default: \code{7}).
#' @param save_path Character; directory to save heatmap if saving (default: current working directory).
#' @param rownames_size Numeric; font size for row names (default: \code{20}).
#' @param colnames_size Numeric; font size for column names (default: \code{20}).
#' @param ... Additional arguments passed to \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap}} object (invisible if saved to PDF).
#'  
#' @export
visualize_network_heatmap <- function(result, 
                                      adjacency_threshold = 10^-6, 
                                      show_values = FALSE,
                                      value_cex = 0.5,
                                      save_pdf = FALSE, 
                                      width = 7, 
                                      height = 7,
                                      save_path = ".",
                                      rownames_size = 20,
                                      colnames_size = 20, 
                                      ...) {
  method <- strsplit(names(result)[1], "\\.")[[1]][2]
  
  if (method %in% c("corr", "pcorr")) {
    cluster_rows = TRUE
    cluster_columns = TRUE
    mat <- as.matrix(result[[1]])
    mat[abs(mat) < adjacency_threshold] <- 0
    title <- paste("Heatmap of", ifelse(method == "corr", "Correlation", "Partial Correlation"))
    col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  } else if (method == "corr_test") {
    cluster_rows = FALSE
    cluster_columns = FALSE
    mat <- as.matrix(result[[1]])
    edge_idx = mat < adjacency_threshold
    mat[edge_idx] <- 1
    mat[!edge_idx] <- 0
    title <- paste0("Heatmap of inferred network at FDR=", adjacency_threshold*100, "%")
    col_fun <- c("0" = "gray", "1" = "blue")
  }
  
  heatmap_args <- list(
    mat,
    name = switch(method,
                  "corr" = "Corr",
                  "pcorr" = "Partial\nCorr",
                  "corr_test" = "FDR"),
    col = col_fun,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_title = title,
    width = grid::unit(width, "in"), 
    height = grid::unit(height, "in"),
    row_names_gp = grid::gpar(fontsize = rownames_size),
    column_names_gp = grid::gpar(fontsize = colnames_size),
    ...
  )
  
  # Add cell_fun if needed
  if (show_values) {
    heatmap_args$cell_fun <- function(j, i, x, y, width, height, fill) {
      grid::grid.text(sprintf("%.3f", mat[i, j]), x, y,
                      gp = grid::gpar(fontsize = 10 * value_cex))
    }
  }
  
  # Add top_annotation only if it's corr_test
  if (method == "corr_test" && !is.null(result$metabolite_clusters)) {
    metabolite_cluster_levels <- sort(unique(result$metabolite_clusters))
    metabolite_cluster_colors <- setNames(RColorBrewer::brewer.pal(length(metabolite_cluster_levels), 
                                                                   "Set1"), 
                                          metabolite_cluster_levels)
    
    heatmap_args$top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      metabolite_clusters = as.factor(result$metabolite_clusters),
      col = list(metabolite_clusters = metabolite_cluster_colors),
      annotation_legend_param = list(
        metabolite_clusters = list(title = "Metabolite Cluster")
      ),
      show_annotation_name=F
    )
    metabolite_order <- order(result$metabolite_clusters)
    heatmap_args$column_order <- metabolite_order
    
    microbiome_cluster_levels <- sort(unique(result$microbiome_clusters))
    microbiome_cluster_colors <- setNames(RColorBrewer::brewer.pal(length(microbiome_cluster_levels), 
                                                                   "Set2"), 
                                          microbiome_cluster_levels)
    
    heatmap_args$left_annotation <- ComplexHeatmap::rowAnnotation(
      microbiome_clusters = as.factor(result$microbiome_clusters),
      col = list(microbiome_clusters = microbiome_cluster_colors),
      annotation_legend_param = list(
        microbiome_clusters = list(title = "Microbiome Cluster")
      ),
      show_annotation_name=F
    )
    microbiome_order <- order(result$microbiome_clusters)
    heatmap_args$row_order <- microbiome_order
    
    heatmap_args$heatmap_legend_param <- list(
      at = c("0", "1"),
      labels = c("Absent", "Present"),
      title = "Edge"
    )
  }
  ht <- do.call(ComplexHeatmap::Heatmap, heatmap_args)
  if (is.character(save_pdf) || save_pdf) {
    if (!is.character(save_pdf)) {
      save_pdf = "plot.pdf"
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

