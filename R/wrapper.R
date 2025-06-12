#' Microbiome-Metabolite Network Analysis Wrapper
#' 
#' @param X A p x n matrix of microbiome abundances (p features, n samples)
#' @param Y A q x n matrix of metabolites (q features, n samples)
#' @param method Character string specifying the method: "conditional", "pc_conditional", "sbm_marginal", or "bisbm_marginal"
#' @param zscore_method Character string specifying z-score normalization for marginal methods: "pearson+fisher" or "CLR"
#' @param lambda A sequence of decreasing positive numbers for regularization (used in pc_conditional method)
#' @param alpha Significance level for graph inference (default: 0.1)
#' @param model Model type for SBM methods (default: "Gauss01")
#' @param sbm_params List of parameters for SBM methods (Q1, Q2, explor)
#' @param init_params List of initialization parameters for SBM methods
#' @param nb_cores Number of cores for parallel processing (default: 1)
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#' 
#' @return A list containing:
#' \describe{
#'   \item{method}{The method used}
#'   \item{correlation_matrix}{Estimated correlation/partial correlation matrix (for conditional methods)}
#'   \item{adjacency_matrix}{Inferred adjacency matrix}
#'   \item{clustering}{Node clustering results (for SBM methods)}
#'   \item{parameters}{Estimated model parameters}
#'   \item{qvalues}{Q-values for edge significance (where applicable)}
#'   \item{model_selection}{Model selection results (for SBM methods)}
#' }
#' 
#' @export
analyze_microbiome_metabolite_network <- function(X, Y, 
                                                  method = c("conditional", "pc_conditional", "sbm_marginal", "bisbm_marginal"),
                                                  zscore_method = c("pearson+fisher", "CLR"),
                                                  lambda = NULL,
                                                  alpha = 0.1,
                                                  model = "Gauss01",
                                                  sbm_params = list(Q1 = 1:5, Q2 = 1:5, explor = 1.5),
                                                  init_params = list(nbOfbeta = 1, nbOfPointsPerbeta = NULL,
                                                                     maxNbOfPasses = 2, minNbOfPasses = 1),
                                                  nb_cores = 1,
                                                  verbose = TRUE) {
  
  # Input validation
  method <- match.arg(method)
  zscore_method <- match.arg(zscore_method)
  
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
  if (method %in% c("conditional", "pc_conditional")) {
    # Conditional methods
    result <- analyze_conditional_methods(X, Y, method, lambda, alpha, verbose, result)
    
  } else if (method %in% c("sbm_marginal", "bisbm_marginal")) {
    # Marginal SBM methods
    result <- analyze_marginal_sbm_methods(X, Y, method, zscore_method, alpha, model, 
                                           sbm_params, init_params, nb_cores, verbose, result)
  }
  
  if (verbose) {
    cat("Analysis completed successfully.\n")
  }
  
  return(result)
}

#' Internal function for conditional methods
analyze_conditional_methods <- function(X, Y, method, lambda, alpha, verbose, result) {
  
  # Combine X and Y matrices (samples as rows for cggm functions)
  combined_data <- t(rbind(X, Y))  # n x (p+q) matrix
  
  if (method == "conditional") {
    if (verbose) cat("Computing correlation matrix using cggm.corr...\n")
    
    # Estimate correlation matrix
    correlation_matrix <- cggm.corr(combined_data)
    result$correlation_matrix <- correlation_matrix
    
    # Create adjacency matrix based on correlation threshold
    # Using Fisher's z-transform for significance testing
    n_samples <- nrow(combined_data)
    z_scores <- 0.5 * log((1 + abs(correlation_matrix)) / (1 - abs(correlation_matrix)))
    p_values <- 2 * pnorm(abs(z_scores) * sqrt(n_samples - 3), lower.tail = FALSE)
    
    # Adjust p-values using FDR
    q_values <- p.adjust(p_values, method = "fdr")
    result$qvalues <- matrix(q_values, nrow = nrow(correlation_matrix))
    
    # Create adjacency matrix
    adjacency_matrix <- matrix(0, nrow = nrow(correlation_matrix), ncol = ncol(correlation_matrix))
    adjacency_matrix[q_values < alpha] <- 1
    diag(adjacency_matrix) <- 0  # Remove self-loops
    
    result$adjacency_matrix <- adjacency_matrix
    
  } else if (method == "pc_conditional") {
    if (verbose) cat("Computing partial correlation matrix using cggm.pcorr...\n")
    
    # Set default lambda if not provided
    if (is.null(lambda)) {
      lambda <- seq(0.5, 0.01, length.out = 20)
    }
    
    # Estimate partial correlation matrices
    pcorr_result <- cggm.pcorr(combined_data, lambda = lambda, method = "glasso")
    
    # Select best lambda using cross-validation or BIC (simplified approach)
    best_lambda_idx <- select_best_lambda(pcorr_result, combined_data)
    
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
analyze_marginal_sbm_methods <- function(X, Y, method, zscore_method, alpha, model, 
                                         sbm_params, init_params, nb_cores, verbose, result) {
  
  if (verbose) cat("Running SBM analysis with method:", method, "\n")
  
  # Run SBM analysis
  sbm_result <- runSBM(X, Y, method, zscore_method, alpha, model, 
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

#' Internal function to run SBM analysis
runSBM <- function(X, Y, method, zscore_method, alpha, model, sbm_params, init_params, nb_cores, verbose) {
  
  # Apply z-score normalization
  if (verbose) cat("Applying z-score normalization method:", zscore_method, "\n")
  
  if (zscore_method == "pearson+fisher") {
    # Pearson correlation + Fisher z-transform
    X_norm <- apply(X, 1, scale)  # standardize each feature
    Y_norm <- apply(Y, 1, scale)
    combined_norm <- rbind(t(X_norm), t(Y_norm))
  } else if (zscore_method == "CLR") {
    # Centered log-ratio transformation
    X_clr <- apply(X, 2, function(x) log(x / exp(mean(log(x[x > 0])))))
    Y_clr <- apply(Y, 2, function(x) log(x / exp(mean(log(x[x > 0])))))
    combined_norm <- rbind(X_clr, Y_clr)
  }
  
  # Create data matrix for SBM
  if (method == "sbm_marginal") {
    # For SBM, create full symmetric matrix
    n_total <- nrow(combined_norm)
    data_matrix <- matrix(0, n_total, n_total)
    
    # Fill with pairwise relationships
    for (i in 1:(n_total-1)) {
      for (j in (i+1):n_total) {
        # Compute correlation or distance metric
        corr_val <- cor(combined_norm[i, ], combined_norm[j, ])
        data_matrix[i, j] <- data_matrix[j, i] <- corr_val
      }
    }
    
    if (verbose) cat("Running fitNSBM...\n")
    sbm_fit <- fitNSBM(data_matrix, 
                       model = model,
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
    
  } else if (method == "bisbm_marginal") {
    # For biSBM, use bipartite structure
    # Create bipartite data matrix (microbiome x metabolite)
    data_matrix <- t(X) %*% Y / ncol(X)  # Simple cross-correlation matrix
    
    if (verbose) cat("Running fitNobiSBM...\n")
    bisbm_fit <- fitNobiSBM(data_matrix,
                            model = model,
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

#' Helper function to select best lambda for partial correlation
select_best_lambda <- function(pcorr_result, data) {
  # Simple BIC-based selection (can be improved with cross-validation)
  n <- nrow(data)
  p <- ncol(data)
  
  bic_values <- sapply(1:length(pcorr_result$lambda), function(i) {
    adj_matrix <- pcorr_result$path[[i]]
    num_edges <- sum(adj_matrix) / 2  # undirected graph
    
    # Compute log-likelihood (simplified)
    precision_matrix <- pcorr_result$icov[[i]]
    log_det <- determinant(precision_matrix, logarithm = TRUE)$modulus
    
    if (is.finite(log_det)) {
      log_likelihood <- n/2 * (log_det - sum(diag(cov(data) %*% precision_matrix)))
      bic <- -2 * log_likelihood + log(n) * num_edges
    } else {
      bic <- Inf
    }
    
    return(bic)
  })
  
  return(which.min(bic_values))
}

#' Placeholder function for graphInferenceNobiSBM (to be implemented based on your package)
graphInferenceNobiSBM <- function(dataMatrix, clustering_row, clustering_col, theta, alpha) {
  # This function should implement graph inference for bipartite SBM
  # Based on the structure shown in your examples
  
  # Placeholder implementation - replace with actual function
  warning("graphInferenceNobiSBM is a placeholder function. Please implement based on your package.")
  
  n1 <- nrow(dataMatrix)
  n2 <- ncol(dataMatrix)
  
  # Simple thresholding as placeholder
  p_values <- 2 * pnorm(abs(as.vector(dataMatrix)), lower.tail = FALSE)
  q_values <- p.adjust(p_values, method = "fdr")
  
  adjacency_matrix <- matrix(0, n1, n2)
  adjacency_matrix[q_values < alpha] <- 1
  
  return(list(
    A = adjacency_matrix,
    qvalues = matrix(q_values, nrow = n1, ncol = n2)
  ))
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
    if (x$method == "sbm_marginal") {
      cat("Best number of blocks (Q):", x$model_selection$best_Q, "\n")
    } else if (x$method == "bisbm_marginal") {
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
  cat("  method = 'bisbm_marginal',\n")
  cat("  zscore_method = 'CLR',\n")
  cat("  alpha = 0.1\n")
  cat(")\n")
  cat("\n")
  cat("# View results\n")
  cat("print(result)\n")
}

