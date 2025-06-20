#' Microbiome-Metabolite Network Analysis Wrapper
#' 
#' @param X A p x n matrix of microbiome abundances (p features, n samples)
#' @param Y A q x n matrix of metabolite abundances (q features, n samples)
#' @param method Character string specifying the method to use.
#' One of: 
#'   \describe{
#'     \item{"corr"}{Marginal correlation}
#'     \item{"pcorr"}{Partial correlation (e.g., graphical lasso)}
#'     \item{"sbm"}{Stochastic Block Model on marginal associations}
#'     \item{"bisbm"}{Bipartite SBM on marginal associations}
#'   }
#' @param zscore_method Character string specifying z-score normalization for marginal SBM methods
#' One of: 
#'   \describe{
#'     \item{"pearson+fisher"}{Uses Pearson correlation, followed by Fisher transform}
#'     \item{"CL"}{Uses the method from Cai and Liu}
#'   }
#' @param lambda A sequence of decreasing positive numbers for regularization (used in pcorr method)
#' @param alpha Significance level for graph inference (default: 0.1)
#' @param sbm_model Character string specifying the model type for SBM methods 
#' (default: \code{"Gauss01"}). For other model options, see the documentation 
#' for \code{fitNSBM} in the \pkg{noisySBM} package.
#' @param sbm_params List of parameters for SBM methods (Q1, Q2, explor)
#' @param init_params List of initialization parameters for SBM methods
#' @param nb_cores Number of cores for parallel processing (default: 1)
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#' 
#' @return A list containing:
#' \describe{
#'   \item{method}{The method used}
#'   \item{correlation_matrix}{Estimated correlation/partial correlation matrix (for correlation methods)}
#'   \item{adjacency_matrix}{Inferred adjacency matrix}
#'   \item{clustering}{Node clustering results (for SBM methods)}
#'   \item{parameters}{Estimated model parameters}
#'   \item{qvalues}{Q-values for edge significance (where applicable)}
#'   \item{model_selection}{Model selection results (for SBM methods)}
#' }
#' 
#' @export
analyze_microbiome_metabolite_network <- function(X, Y, 
                                                  method,
                                                  lambda = NULL,
                                                  alpha = 0.1,
                                                  zscore_method = NULL,
                                                  sbm_model = NULL,
                                                  sbm_params = NULL,
                                                  init_params = NULL,
                                                  nb_cores = 1,
                                                  verbose = TRUE) {
  
  # Input validation and default values
  method <- match.arg(method, c("corr", "pcorr", "sbm", "bisbm"))
  
  if (method != "pcorr") {
    if (!is.null(lambda)) {
      warning("lambda is only used when method = 'pcorr'; the provided value will be ignored.")
    }
  }
  
  if (method %in% c("corr", "pcorr")) {
    if (!is.null(zscore_method)) {
      warning("zscore_method is only used when method = 'sbm' or 'bisbm'; the provided value will be ignored.")
    }
    if (!is.null(sbm_model)) {
      warning("sbm_model is only used when method = 'sbm' or 'bisbm'; the provided value will be ignored.")
    }
    if (!is.null(sbm_params)) {
      warning("sbm_params is only used when method = 'sbm' or 'bisbm'; the provided value will be ignored.")
    }
    if (!is.null(init_params)) {
      warning("init_params is only used when method = 'sbm' or 'bisbm'; the provided value will be ignored.")
    }
  }
  
  if (method %in% c("sbm", "bisbm")) {
    if (is.null(zscore_method)) {
      # default
      zscore_method <- "pearson+fisher"
      message("`zscore_method` not provided, using default: '", zscore_method, "'.")
    } else {
      zscore_method <- match.arg(zscore_method, c("pearson+fisher", "CL"))
    }
  }

  if (verbose) {
    cat("Starting microbiome-metabolite network analysis with method:", method, "\n")
  }
  
  # Check matrix dimensions
  if (ncol(X) != ncol(Y)) {
    stop("X and Y must have the same number of samples (columns)")
  }
  
  p <- nrow(X)  # number of microbiome features
  q <- nrow(Y)  # number of metabolite features
  n <- ncol(X)  # number of samples
  
  if (verbose) {
    cat("Data dimensions: p =", p, ", q =", q, ", n =", n, "\n")
  }
  
  # Initialize result list
  result <- list(
    method = method,
    data_dimensions = list(p = p, q = q, n = n),
    parameters = list()
  )
  
  # Method-specific analysis
  if (method %in% c("corr", "pcorr")) {
    # Correlation methods
    result <- run_correlation_method(X, Y, method, lambda, alpha, verbose, result)
    
  } else if (method %in% c("sbm", "bisbm")) {
    # SBM methods
    result <- run_sbm_method(X, Y, method, zscore_method, alpha, sbm_model, 
                             sbm_params, init_params, nb_cores, verbose, result)
  }
  
  if (verbose) {
    cat("Analysis completed successfully.\n")
  }
  
  return(result)
}

#' Internal function for correlation methods
run_correlation_method <- function(X, Y, method, lambda, alpha, verbose, result) {
  
  # Combine X and Y matrices (samples as rows for cggm functions)
  combined_data <- t(rbind(X, Y))  # n x (p+q) matrix
  
  if (method == "corr") {
    if (verbose) cat("Computing correlation matrix using cggm.corr...\n")
    
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
    if (verbose) cat("Computing partial correlation matrix using cggm.pcorr...\n")
    
    # Set default lambda if not provided
    if (is.null(lambda)) {
      lambda <- seq(0.5, 0.01, length.out = 20)
    }
    
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

#' Internal function for marginal SBM methods
run_sbm_method <- function(X, Y, method, zscore_method, alpha, sbm_model, 
                           sbm_params, init_params, nb_cores, verbose, result) {
  
  if (verbose) cat("Running SBM analysis with method:", method, "\n")
  
  # Run SBM analysis
  sbm_result <- runSBM(X, Y, method, zscore_method, alpha, sbm_model, 
                       sbm_params, init_params, nb_cores, verbose)
  
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

#' Internal function to compute CL test statistic
compute_CL_test_statistic <- function(X, Y) {
  # X: p x n1 matrix
  # Y: p x n2 matrix
  p <- nrow(X)
  n1 <- ncol(X)
  n2 <- ncol(Y)
  
  # Standardize X and Y
  X_centered <- scale(t(X), center = TRUE, scale = TRUE) # n1 x p
  Y_centered <- scale(t(Y), center = TRUE, scale = TRUE) # n2 x p
  
  # Compute sample correlations
  rho_X <- cor(X_centered)
  rho_Y <- cor(Y_centered)
  
  T_mat <- matrix(0, p, p)
  
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      rho1 <- rho_X[i, j]
      rho2 <- rho_Y[i, j]
      
      # Compute theta_hat for X
      x_i <- X_centered[, i]
      x_j <- X_centered[, j]
      theta1 <- mean((2 * x_i * x_j - rho1 * x_i^2 - rho1 * x_j^2)^2)
      
      # Compute theta_hat for Y
      y_i <- Y_centered[, i]
      y_j <- Y_centered[, j]
      theta2 <- mean((2 * y_i * y_j - rho2 * y_i^2 - rho2 * y_j^2)^2)
      
      # Compute T'
      T_val <- 2 * (rho1 - rho2) / sqrt(theta1 / n1 + theta2 / n2)
      T_mat[i, j] <- T_val
      T_mat[j, i] <- T_val
    }
  }
  return(T_mat)
}

#' Internal function to run SBM analysis
runSBM <- function(X, Y, method, zscore_method, alpha, sbm_model, sbm_params, init_params, nb_cores, verbose) {
  combined_data <- t(rbind(X, Y))  # dimension: n x (p + q)
  n_features <- ncol(combined_data)
  data_matrix <- matrix(0, n_features, n_features)
  
  # Fill with pairwise relationships between features
  for (i in 1:(n_features - 1)) {
    for (j in (i + 1):n_features) {
      # Compute correlation between features i and j (across samples)
      if (zscore_method == "spearman+fisher") {
        corr_val <- cor(combined_data[, i], combined_data[, j], method = "spearman")
      } else if (zscore_method == "pearson+fisher") {
        corr_val <- cor(combined_data[, i], combined_data[, j], method = "pearson")
      }
      
      # Fisher transformation
      z_val <- 0.5 * log((1 + corr_val) / (1 - corr_val))
      
      data_matrix[i, j] <- data_matrix[j, i] <- z_val
    }
  }
  
  if (method == "sbm") { # (p+q) x (p+q)
    if (verbose) cat("Running fitNSBM...\n")
    sbm_fit <- fitNSBM(data_matrix, 
                       model = sbm_model,
                       sbmSize = list(Qmin = min(sbm_params$Q1), 
                                      Qmax = max(sbm_params$Q1), 
                                      explor = sbm_params$explor),
                       initParam = init_params)
    
    # Perform graph inference
    best_idx <- which.max(sapply(sbm_fit, function(x) x$sbmParam$ICL))
    best_solution <- sbm_fit[[best_idx]]
    
    graph_result <- noisySBM::graphInference(data_matrix,
                                             best_solution$clustering,
                                             best_solution$theta, 
                                             alpha = alpha)
    
    result <- list(
      sbm_fit = sbm_fit,
      adjacency_matrix = graph_result$A,
      qvalues = graph_result$qvalues,
      clustering = list(nodes = best_solution$clustering),
      parameters = best_solution$theta,
      model_selection = list(
        best_Q = best_solution$sbmParam$Q,
        ICL = best_solution$sbmParam$ICL,
        best_index = best_idx
      )
    )
    
  } else if (method == "bisbm") { # (p x q)
    # For biSBM, use bipartite structure
    if (verbose) cat("Running fitNobiSBM...\n")
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
  }
  
  return(result)
}


# # Apply z-score normalization
# if (verbose) cat("Applying z-score normalization method:", zscore_method, "\n")
# 
# if (zscore_method == "pearson+fisher") {
#   # Pearson correlation + Fisher z-transform
#   X_norm <- apply(X, 1, scale)  # standardize each feature
#   Y_norm <- apply(Y, 1, scale)
#   combined_norm <- rbind(t(X_norm), t(Y_norm))
# } else if (zscore_method == "CL") {
#   # combined_norm <- rbind(X_clr, Y_clr)
# }
# 

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

# Example usage function
#' @export
example_usage <- function() {
  cat("Example usage:\n")
  cat("# Generate example data\n")
  cat("set.seed(123)\n")
  cat("n <- 50  # samples\n")
  cat("p <- 20  # microbiome features\n") 
  cat("q <- 15  # metabolite features\n")
  cat("X <- matrix(abs(rnorm(p*n)), nrow=p, ncol=n)\n")
  cat("Y <- matrix(abs(rnorm(q*n)), nrow=q, ncol=n)\n")
  cat("\n")
  cat("# Run analysis\n")
  cat("result <- analyze_microbiome_metabolite_network(\n")
  cat("  X = X, Y = Y,\n")
  cat("  method = 'bisbm',\n")
  cat("  zscore_method = 'CL',\n")
  cat("  alpha = 0.1\n")
  cat(")\n")
  cat("\n")
  cat("# View results\n")
  cat("print(result)\n")
}

