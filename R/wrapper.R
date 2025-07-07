#' Microbiome-Metabolite Network Analysis Wrapper
#' 
#' @param X A p x n matrix of microbiome abundances (p features, n samples)
#' @param Y A q x n matrix of metabolite abundances (q features, n samples)
#' @param method Character string specifying the method to use.
#' One of: 
#'   \describe{
#'     \item{"corr"}{Marginal correlation}
#'     \item{"pcorr"}{Partial correlation (e.g., graphical lasso)}
#'     \item{"bisbm"}{Bipartite SBM on marginal associations}
#'   }
#' @param zscore_method Character string specifying z-score normalization for marginal SBM method
#' One of: 
#'   \describe{
#'     \item{"pearson+fisher"}{Uses Pearson correlation, followed by Fisher transform}
#'     \item{"spearman+fisher"}{Uses Spearman correlation, followed by Fisher transform}
#'     \item{"CL"}{Uses the method from Cai and Liu}
#'   }
#' @param lambda A sequence of decreasing positive numbers for regularization (used in pcorr method)
#' @param alpha Significance level for graph inference (default: 0.05)
#' @param sbm_model Character string specifying the model type for SBM method 
#' (default: \code{"Gauss"}). For other model options, see the documentation 
#' for \code{fitNSBM} in the \pkg{noisySBM} package.
#' @param sbm_params List of parameters for SBM method (Q1, Q2, explor)
#' @param init_params List of initialization parameters for SBM method
#' @param nb_cores Number of cores for parallel processing (default: 1)
#' 
#' @return A named list with method-dependent contents:
#' 
#' \strong{If \code{method = "corr"}}:
#' \describe{
#'   \item{correlation_matrix}{Estimated Pearson correlation matrix (p+q) x (p+q)}
#'   \item{adjacency_matrix}{Binary matrix indicating significant correlations after FDR correction}
#'   \item{qvalues}{FDR-adjusted p-values (same dimensions as correlation matrix)}
#' }
#' 
#' \strong{If \code{method = "pcorr"}}:
#' \describe{
#'   \item{correlation_matrix}{Estimated covariance matrix (p+q) x (p+q)}
#'   \item{partial_correlation_matrix}{Estimated partial correlation matrix, computed from the precision matrix}
#'   \item{adjacency_matrix}{Binary matrix indicating stable edges across regularization path}
#'   \item{lambda_path}{Sequence of lambda values used for graphical lasso}
#'   \item{best_lambda}{Lambda selected via stability selection}
#'   \item{parameters}{List containing lambda selection details (best index and path)}
#' }
#' 
#' \strong{If \code{method = "bisbm"}}:
#' \describe{
#'   \item{adjacency_matrix}{Binary matrix indicating inferred edges between microbiome and metabolite features}
#'   \item{qvalues}{Matrix of FDR-adjusted p-values for each cross-feature test}
#'   \item{clustering}{List with clustering assignments for rows (microbiome) and columns (metabolites)}
#'   \item{parameters}{Estimated SBM parameters (e.g., block means)}
#'   \item{model_selection}{List with selected model dimensions (Q1, Q2) and ICL score}
#'   \item{sbm_result}{Raw SBM fit object (all candidate models)}
#'   \item{zscore_method}{Z-score normalization method used ('CL', 'pearson+fisher', 'spearman+fisher')}
#' }
#' 
#' \strong{Common fields}:
#' \describe{
#'   \item{method}{Character string indicating the method used}
#'   \item{data_dimensions}{List with \code{p}, \code{q}, and \code{n} (numbers of microbiome features, metabolite features, and samples)}
#' }
#' 
#' @export
analyze_microbiome_metabolite_network <- function(X, Y, 
                                                  method,
                                                  lambda = NULL,
                                                  alpha = 0.05,
                                                  zscore_method = NULL,
                                                  sbm_model = NULL,
                                                  sbm_params = NULL,
                                                  init_params = NULL,
                                                  nb_cores = 1) {
  
  # Input validation and default values
  valid_methods <- c("corr", "pcorr", "bisbm")
  
  if (!method %in% valid_methods) {
    stop(sprintf("Invalid method '%s'. Must be one of: %s", method, paste(valid_methods, collapse = ", ")))
  }
  
  # Check matrix dimensions
  if (ncol(X) != ncol(Y)) {
    stop("X and Y must have the same number of samples (columns)")
  }
  
  if (method != "pcorr") {
    if (!is.null(lambda)) {
      warning("lambda is only used when method = 'pcorr'; the provided value will be ignored.")
    }
  }
  
  if (method %in% c("corr", "pcorr")) {
    if (!is.null(zscore_method)) {
      warning("zscore_method is only used when method = 'bisbm'; the provided value will be ignored.")
    }
    if (!is.null(sbm_model)) {
      warning("sbm_model is only used when method = 'bisbm'; the provided value will be ignored.")
    }
    if (!is.null(sbm_params)) {
      warning("sbm_params is only used when method = 'bisbm'; the provided value will be ignored.")
    }
    if (!is.null(init_params)) {
      warning("init_params is only used when method = 'bisbm'; the provided value will be ignored.")
    }
    
    if (method == "pcorr") {
      # Set default lambda if not provided
      if (is.null(lambda)) {
        message("`lambda` not provided, using default sequence: ",
                "seq(0.1, 3, 0.1) * sqrt(log(p) / n)")
        lambda <- seq(0.1, 3, 0.1) * sqrt(log(nrow(X) + nrow(Y)) / ncol(X))
      }
    }
  }
  
  if (method == "bisbm") {
    if (is.null(zscore_method)) {
      zscore_method <- "pearson+fisher"
      message("`zscore_method` not provided, using default: '", zscore_method, "'.")
    } else {
      valid_zscore_methods = c("pearson+fisher", "spearman+fisher", "CL")
      if (!zscore_method %in% valid_zscore_methods) {
        stop(sprintf("Invalid zscore_method '%s'. Must be one of: %s", method, paste(valid_zscore_methods, collapse = ", ")))
      }
    }
    
    if (is.null(sbm_model)) {
      sbm_method <- "Gauss"
      message("`sbm_method` not provided, using default: '", sbm_method, "'.")
    } else {
      valid_sbm_methods = c('Gauss','Gauss0','Gauss01','GaussEqVar',
                            'Gauss0EqVar','Gauss0Var1','Gauss02distr',
                            'Gauss2distr','GaussAffil', 'Exp','ExpGamma','Gamma')
      if (!sbm_method %in% valid_sbm_methods) {
        stop(sprintf("Invalid sbm_method '%s'. Must be one of: %s", sbm_method, paste(valid_sbm_methods, collapse = ", ")))
      }
    }
    
    if (!is.null(sbm_params)) {
      required_fields <- c("Q1", "Q2", "explor")
      missing_fields <- setdiff(required_fields, names(sbm_params))
      if (length(missing_fields) > 0) {
        stop("`sbm_params` must contain: Q1, Q2, explor. Missing: ", paste(missing_fields, collapse = ", "))
      }
      if (!is.numeric(sbm_params$Q1) || !is.numeric(sbm_params$Q2)) {
        stop("`sbm_params$Q1` and `sbm_params$Q2` must be numeric vectors.")
      }
      if (!is.numeric(sbm_params$explor) || length(sbm_params$explor) != 1) {
        stop("`sbm_params$explor` must be a single numeric value.")
      }
    } else {
      sbm_params <- list(Q1 = 1:5, Q2 = 1:5, explor = 1.5)
      message("`sbm_params` not provided, using default: Q1 = 1:5, Q2 = 1:5, explor = 1.5")
    }
    
    if (!is.null(init_params)) {
      required_fields <- c("nbOfbeta", "nbOfPointsPerbeta", "maxNbOfPasses", "minNbOfPasses")
      missing_fields <- setdiff(required_fields, names(init_params))
      if (length(missing_fields) > 0) {
        stop("`init_params` must contain: ", paste(required_fields, collapse = ", "), 
             ". Missing: ", paste(missing_fields, collapse = ", "))
      }
      if (!is.numeric(init_params$nbOfbeta) || init_params$nbOfbeta < 1) {
        stop("`init_params$nbOfbeta` must be a numeric value >= 1.")
      }
      if (!is.numeric(init_params$nbOfPointsPerbeta) || init_params$nbOfPointsPerbeta < 1) {
        stop("`init_params$nbOfPointsPerbeta` must be a numeric value >= 1.")
      }
      if (!is.numeric(init_params$maxNbOfPasses) || init_params$maxNbOfPasses < init_params$minNbOfPasses) {
        stop("`init_params$maxNbOfPasses` must be a numeric value >= minNbOfPasses.")
      }
      if (!is.numeric(init_params$minNbOfPasses) || init_params$minNbOfPasses < 1) {
        stop("`init_params$minNbOfPasses` must be a numeric value >= 1.")
      }
    } else {
      init_params <- list(nbOfbeta=NULL, nbOfPointsPerbeta=NULL,
                          maxNbOfPasses=NULL, minNbOfPasses=1)
      message("`init_params` not provided, using default: nbOfbeta=NULL, nbOfPointsPerbeta=NULL, maxNbOfPasses=NULL, minNbOfPasses=1")
    }
  }
  
  # if (verbose) {
  cat("Starting microbiome-metabolite network analysis with method:", method, "\n")
  # }
  
  p <- nrow(X)  # number of microbiome features
  q <- nrow(Y)  # number of metabolite features
  n <- ncol(X)  # number of samples
  
  #if (verbose) {
  cat("Data dimensions: p =", p, ", q =", q, ", n =", n, "\n")
  #}
  
  # Initialize result list
  result <- list(
    method = method,
    data_dimensions = list(p = p, q = q, n = n)
  )
  
  # Method-specific analysis
  if (method %in% c("corr", "pcorr")) {
    # Correlation methods
    result <- run_correlation_method(X, Y, method, lambda, alpha, result)
    
  } else if (method == "bisbm") {
    # SBM method
    result <- run_sbm_method(X, Y, method, zscore_method, alpha, sbm_model, 
                             sbm_params, init_params, nb_cores, result)
  }
  
  #if (verbose) {
  cat("Analysis completed successfully.\n")
  #}
  
  return(result)
}

#' Internal function for correlation methods
run_correlation_method <- function(X, Y, method, lambda, alpha, result) {
  
  # Combine X and Y matrices (samples as rows for cggm functions)
  combined_data <- t(rbind(X, Y))  # n x (p+q) matrix
  
  if (method == "corr") {
    #if (verbose) 
    cat("Computing correlation matrix using cggm.corr...\n")
    
    # Estimate correlation matrix
    correlation_matrix <- cggm.corr(combined_data)
    result$correlation_matrix <- correlation_matrix
    
    # Create adjacency matrix based on correlation threshold
    # Using Fisher's z-transform for significance testing
    n_samples <- nrow(combined_data)
    z_scores <- 0.5 * log((1 + correlation_matrix) / (1 - correlation_matrix))
    p_values <- pnorm(z_scores * sqrt(n_samples - 3), lower.tail = FALSE)
    
    # Adjust p-values using FDR
    q_values <- p.adjust(p_values, method = "fdr")
    result$qvalues <- matrix(q_values, nrow = nrow(correlation_matrix))
    
    # Create adjacency matrix
    adjacency_matrix <- matrix(0, nrow = nrow(correlation_matrix), ncol = ncol(correlation_matrix))
    adjacency_matrix[q_values < alpha] <- 1
    diag(adjacency_matrix) <- 0  # Remove self-loops
    
    result$adjacency_matrix <- adjacency_matrix
    
  } else if (method == "pcorr") {
    cat("Computing partial correlation matrix using cggm.pcorr...\n")
    
    # Estimate partial correlation matrices
    pcorr_result <- cggm.pcorr(combined_data, lambda = lambda, method = "glasso")
    
    # Select best lambda using cggm.stars
    best_lambda <- cggm.stars(pcorr_result)
    best_lambda_idx <- best_lambda$opt.index
      
    result$correlation_matrix <- pcorr_result$cov
    result$partial_correlation_matrix <- solve(pcorr_result$icov[[best_lambda_idx]])
    result$adjacency_matrix <- pcorr_result$path[[best_lambda_idx]]
    result$lambda_path <- pcorr_result$lambda
    result$best_lambda <- pcorr_result$lambda[best_lambda_idx]
    result$parameters$lambda_selection <- list(
      best_index = best_lambda_idx,
      lambda_values = pcorr_result$lambda
    )
  }
  
  return(result)
}

#' Internal function for marginal SBM method
run_sbm_method <- function(X, Y, method, zscore_method, alpha, sbm_model, 
                           sbm_params, init_params, nb_cores, result) {
  
  #if (verbose) 
  cat("Running SBM analysis\n")
  
  # Run SBM analysis
  sbm_result <- runSBM(X, Y, method, zscore_method, alpha, sbm_model, 
                       sbm_params, init_params, nb_cores)
  
  # Store results
  result$sbm_result <- sbm_result$sbm_fit
  result$adjacency_matrix <- sbm_result$adjacency_matrix
  result$qvalues <- sbm_result$qvalues
  result$clustering <- sbm_result$clustering
  result$parameters <- sbm_result$parameters
  result$model_selection <- sbm_result$model_selection
  result$zscore_method <- zscore_method
  
  return(result)
}

#' Internal function to compute CL test statistic (Cai & Liu method)
#'
#' Computes test statistics T_ij and associated p-values for testing cov(X_i, Y_j) = 0
#' using method from Cai & Liu (2016), see T-statistic page 238.
#'
#' @param X A p x n matrix (microbiome features x samples)
#' @param Y A q x n matrix (metabolite features x samples)
#' @return A list with test_stat (T_ij), p_values, and q_values
compute_CL_test_statistic <- function(X, Y) {
  p <- nrow(X)
  q <- nrow(Y)
  n <- ncol(X)
  
  X_centered <- X - rowMeans(X)
  Y_centered <- Y - rowMeans(Y)
  
  T_mat <- matrix(0, nrow = p, ncol = q)
  pval_mat <- matrix(1, nrow = p, ncol = q)
  
  for (i in 1:p) {
    for (j in 1:q) {
      x <- X_centered[i, ]
      y <- Y_centered[j, ]
      sigma_hat <- mean(x * y)  # sample covariance
      
      theta_hat <- mean((x * y - sigma_hat)^2)
      
      T_ij <- sqrt(n) * sigma_hat / sqrt(theta_hat)
      p_val <- 2 * pnorm(abs(T_ij), lower.tail = FALSE)
      
      T_mat[i, j] <- T_ij
      pval_mat[i, j] <- p_val
    }
  }
  
  qval_mat <- matrix(p.adjust(as.vector(pval_mat), method = "fdr"), nrow = p, ncol = q)
  
  return(list(
    test_stat = T_mat,
    p_values = pval_mat,
    q_values = qval_mat
  ))
}

#' Internal function to run SBM analysis
runSBM <- function(X, Y, method, zscore_method, alpha, sbm_model, sbm_params, 
                   init_params, nb_cores) {
  
  if (zscore_method == "CL") {
    #if (verbose) 
    cat("Computing CL test statistics...\n")
    cl_result <- compute_CL_test_statistic(X, Y)
    data_matrix <- cl_result$test_stat 
  } else {
    p <- nrow(X)
    q <- nrow(Y)
    n <- ncol(X)
    
    data_matrix <- matrix(0, nrow = p, ncol = q)  # p x q matrix for cross-correlations
    
    # Compute cross-correlations between X and Y
    for (i in 1:p) {
      for (j in 1:q) {
        if (zscore_method == "spearman+fisher") {
          corr_val <- cor(X[i, ], Y[j, ], method = "spearman")
        } else if (zscore_method == "pearson+fisher") {
          corr_val <- cor(X[i, ], Y[j, ], method = "pearson")
        } else {
          stop("Unsupported zscore_method: ", zscore_method)
        }
        
        # Fisher transformation
        z_val <- 0.5 * log((1 + corr_val) / (1 - corr_val))
        data_matrix[i, j] <- z_val
      }
    }
  }
  
  #if (verbose) 
  cat("Running fitNobiSBM...\n")
  bisbm_fit <- fitNobiSBM(data_matrix,
                          model = sbm_model,
                          exclusive_Rows_Cols = FALSE,
                          sbmSize = sbm_params,
                          initParam = init_params,
                          nbCores = nb_cores)
  
  # Select best model
  best_idx <- which.max(sapply(bisbm_fit, function(x) x$sbmParam$ICL))
  best_solution <- bisbm_fit[[best_idx]]
  
  # Perform graph inference
  graph_result <- graphInferenceNobiSBM(data_matrix,
                                        best_solution$clustering_row,
                                        best_solution$clustering_col,
                                        best_solution$theta,
                                        alpha = alpha)
  
  result <- list(
    sbm_fit = bisbm_fit,
    adjacency_matrix = graph_result$A,
    qvalues = graph_result$qvalues,
    clustering = list(
      rows = best_solution$clustering_row,
      cols = best_solution$clustering_col
    ),
    parameters = best_solution$theta,
    model_selection = list(
      best_Q1 = best_solution$sbmParam$Q1,
      best_Q2 = best_solution$sbmParam$Q2,
      ICL = best_solution$sbmParam$ICL,
      best_index = best_idx
    )
  )
  
  return(result)
}

#' Print method for network analysis results
#' @export
print.microbiome_network_result <- function(x, ...) {
  cat("Microbiome-Metabolite Network Analysis Results\n")
  cat("==============================================\n")
  cat("Method:", x$method, "\n")
  cat("Data dimensions: p =", x$data_dimensions$p, 
      ", q =", x$data_dimensions$q, 
      ", n =", x$data_dimensions$n, "\n")
  
  if (!is.null(x$adjacency_matrix)) {
    cat("Number of inferred edges:", sum(x$adjacency_matrix), "\n")
  }
  
  if (!is.null(x$model_selection)) {
    if (x$method == "sbm") {
      cat("Best number of blocks (Q):", x$model_selection$best_Q, "\n")
    } else if (x$method == "bisbm") {
      cat("Best number of blocks: Q1 =", x$model_selection$best_Q1, 
          ", Q2 =", x$model_selection$best_Q2, "\n")
    }
    cat("ICL:", x$model_selection$ICL, "\n")
  }
}
