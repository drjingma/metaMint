#' Microbiome-Metabolite Network Analysis Wrapper
#' 
#' @param X A p x n matrix of microbiome abundances (p features, n samples)
#' @param Y A q x n matrix of metabolite abundances (q features, n samples)
#' @param method Character string specifying the method to use.
#' One of: 
#'   \describe{
#'     \item{"corr"}{Marginal correlation}
#'     \item{"pcorr"}{Partial correlation (e.g., graphical lasso)}
#'     \item{"corr_test"}{Marginal correlation after accounting for multiple comparisons}
#'   }
#' @param zscore_method Character string specifying z-score normalization for corr_test method
#' One of: 
#'   \describe{
#'     \item{"pearson+fisher"}{Uses Pearson correlation, followed by Fisher transform}
#'     \item{"spearman+fisher"}{Uses Spearman correlation, followed by Fisher transform}
#'     \item{"CL"}{Uses the method from Cai and Liu}
#'   }
#' @param lambda A sequence of decreasing positive numbers for regularization (used in pcorr method)
#' @param alpha Significance level for graph inference (default: 0.05)
#' @param sbm_model Character string specifying the model type for corr_test method 
#' One of: 
#'   \describe{
#'     \item{"Gauss0"}{the mean of the null distribution is set to 0}
#'     \item{"Gauss01"}{the null distribution is set to N(0,1)}
#'     \item{"Gauss02distr"}{the mean of the null distribution is set to 0 and the alternative distribution is a single Gaussian distribution, i.e. the block memberships of the nodes do not influence on the alternative distribution}
#'   }
#' @param sbm_params List of parameters for corr_test method (Q1, Q2, explor)
#' @param nb_cores Number of cores for parallel processing for corr_test method
#' 
#' @return A named list with method-dependent result
#' 
#' \strong{If \code{method = "corr"}}:
#' \describe{
#'   \item{result}{Estimated correlation matrix (p+q) x (p+q)}
#' }
#' 
#' \strong{If \code{method = "pcorr"}}:
#' \describe{
#'   \item{result}{Estimated partial correlation matrix (p+q) x (p+q)}
#' }
#' 
#' \strong{If \code{method = "corr_test"}}:
#' \describe{
#'   \item{qvalues}{(p x q) matrix of FDR-adjusted p-values for each cross-feature test}
#' }
#' 
#' @export
analyze_microbiome_metabolite_network <- function(X, Y, 
                                                  method,
                                                  lambda = NULL,
                                                  alpha = NULL,
                                                  zscore_method = NULL,
                                                  sbm_model = NULL,
                                                  sbm_params = NULL,
                                                  nb_cores = NULL) {
  
  # Input validation and default values
  valid_methods <- c("corr", "pcorr", "corr_test")
  
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
      warning("zscore_method is only used when method = 'corr_test'; the provided value will be ignored.")
    }
    if (!is.null(sbm_model)) {
      warning("sbm_model is only used when method = 'corr_test'; the provided value will be ignored.")
    }
    if (!is.null(sbm_params)) {
      warning("sbm_params is only used when method = 'corr_test'; the provided value will be ignored.")
    }
    if (!is.null(nb_cores)) {
      warning("nb_cores is only used when method = 'corr_test'; the provided value will be ignored.")
    }
    if (!is.null(alpha)) {
      warning("alpha is only used when method = 'corr_test'; the provided value will be ignored.")
    }
    

    if (method == "pcorr") {
      # Set default lambda if not provided
      if (is.null(lambda)) {
        message("`lambda` not provided, using default sequence: seq(0.1, 3, 0.1) * sqrt(log(p) / n)")
        lambda <- seq(0.1, 3, 0.1) * sqrt(log(nrow(X) + nrow(Y)) / ncol(X))
      }
    }
  }
  
  if (method == "corr_test") {
    if (is.null(zscore_method)) {
      zscore_method <- "CL"
      message(sprintf("`zscore_method` not provided, using default: '%s'.", zscore_method))
    } else {
      valid_zscore_methods = c("pearson+fisher", "spearman+fisher", "CL")
      if (!zscore_method %in% valid_zscore_methods) {
        stop(sprintf("Invalid zscore_method '%s'. Must be one of: %s", method, paste(valid_zscore_methods, collapse = ", ")))
      }
    }
    
    if (is.null(nb_cores)) {
      nb_cores = parallel::detectCores()
      message(sprintf("`nb_cores` not provided, using default: '%s'.", nb_cores))
    }
    
    if (is.null(sbm_model)) {
      sbm_model <- "Gauss0"
      message(sprintf("`sbm_model` not provided, using default: '%s'.", sbm_model))
    } else {
      valid_sbm_models = c('Gauss0','Gauss01','Gauss02distr')
      if (!sbm_model %in% valid_sbm_models) {
        stop(sprintf("Invalid sbm_method '%s'. Must be one of: %s", sbm_model, 
                     paste(valid_sbm_models, collapse = ", ")))
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
  }

  cat("Starting microbiome-metabolite network analysis with method:", method, "\n")

  p <- nrow(X)  # number of microbiome features
  q <- nrow(Y)  # number of metabolite features
  n <- ncol(X)  # number of samples
  
  cat("Data dimensions: p =", p, ", q =", q, ", n =", n, "\n")
  
  # Method-specific analysis
  if (method %in% c("corr", "pcorr")) {
    # Correlation methods
    result <- run_correlation_method(X, Y, method, lambda)
  } else if (method == "corr_test") {
    # corr_test method
    result <- run_corr_test_method(X, Y, zscore_method, alpha, sbm_model, 
                                   sbm_params, nb_cores)
  }
  
  cat("Analysis completed successfully.\n")

  return(result)
}

#' Internal function for correlation methods
run_correlation_method <- function(X, Y, method, lambda) {
  
  # Combine X and Y matrices (samples as rows for cggm functions)
  combined_data <- t(rbind(X, Y))  # n x (p+q) matrix
  
  if (method == "corr") {
    cat("Computing correlation matrix using cggm.corr...\n")
    
    # Estimate correlation matrix
    correlation_matrix <- cggm.corr(combined_data)
    result = list(result=correlation_matrix)
  } else if (method == "pcorr") {
    cat("Computing partial correlation matrix using cggm.pcorr...\n")
    
    # Estimate partial correlation matrices
    pcorr_result <- cggm.pcorr(combined_data, lambda = lambda, method = "glasso")
    
    # Select best lambda using cggm.stars
    best_lambda <- cggm.stars(pcorr_result)
    best_lambda_idx <- best_lambda$opt.index
      
    partial_correlation_matrix <- solve(pcorr_result$icov[[best_lambda_idx]])
    result = list(result=partial_correlation_matrix)
  }
  
  return(result)
}

#' Internal function for marginal corr_test method
run_corr_test_method <- function(X, Y, zscore_method, alpha, sbm_model, 
                                 sbm_params, nb_cores) {
  
  cat("Running SBM analysis\n")
  
  # Run SBM analysis
  sbm_result <- runSBM(X, Y, zscore_method, alpha, sbm_model, sbm_params, 
                       nb_cores)
  return(sbm_result)
}

#' Internal function to compute CL test statistic (Cai & Liu method)
#'
#' Computes test statistics T_ij for testing cov(X_i, Y_j) = 0
#' using method from Cai & Liu (2016), see T-statistic page 238.
#'
#' @param X A p x n matrix (microbiome features x samples)
#' @param Y A q x n matrix (metabolite features x samples)
#' @return test statistic matrix 
compute_CL_test_statistic <- function(X, Y) {
  p <- nrow(X)
  q <- nrow(Y)
  n <- ncol(X)
  
  X_centered <- X - rowMeans(X)
  Y_centered <- Y - rowMeans(Y)
  
  T_mat <- matrix(0, nrow = p, ncol = q)

  for (i in 1:p) {
    for (j in 1:q) {
      x <- X_centered[i, ]
      y <- Y_centered[j, ]
      sigma_hat <- mean(x * y)  # sample covariance
      
      theta_hat <- mean((x * y - sigma_hat)^2)
      
      T_ij <- sqrt(n) * sigma_hat / sqrt(theta_hat)
      T_mat[i, j] <- T_ij
    }
  }
  
  return(T_mat)
}

#' Internal function to run SBM analysis
runSBM <- function(X, Y, zscore_method, alpha, sbm_model, sbm_params, nb_cores) {
  if (zscore_method == "CL") {
    cat("Computing CL test statistics...\n")
    data_matrix <- compute_CL_test_statistic(X, Y)
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
  
  cat("Running fitNobiSBM...\n")
  bisbm_fit <- fitNobiSBM(data_matrix,
                          model = sbm_model,
                          exclusive_Rows_Cols = FALSE,
                          sbmSize = sbm_params,
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
  
  result <- list(result = graph_result$qvalues)
  return(result)
}