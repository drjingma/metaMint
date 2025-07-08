# # Create adjacency matrix based on correlation threshold
# adjacency_matrix <- matrix(0, nrow = nrow(correlation_matrix), ncol = ncol(correlation_matrix))
# adjacency_matrix[correlation_matrix < corr_threshold] <- 1
# diag(adjacency_matrix) <- 0  # Remove self-loops
# 
# result$adjacency_matrix <- adjacency_matrix
