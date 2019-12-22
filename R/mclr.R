#' Function to transform microbial count data
#' @param Z A sample by species count matrix
#' @param base The base of logarithm transformation. The default value is \code{e}.
#' @param eps A small constant added to adjust negative values after centered log-ratio transformation. See details.
#' @description This function transforms microbial count data using the modified centered log-ratio (mclr) transformation.
#' @details
#' The centered log-ratio (clr) transformation is commonly used when working with microbiome data that is subject to total sum constraint.
#' In order to apply the clr transformation, one needs to have a sequence of counts that are strictly positive. When there are zeros,
#' a small pseudocount is necessary to ensure all values are strictly positive. However, this may unfairly bias rare variables. The modified clr (mclr)
#' addresses this limitation by applying the clr to positive counts only and then shifting all values to be positive
#' with a small constant \code{eps}. The small constant \code{eps} is defined as the absolute value of the smallest values
#' after clr transformation to the positive counts. For more details, see Yoon et al. (2019).
#' @return A data matrix of the same size as \code{Z} after the modified centered log-ratio transformation
#' @references
#' Yoon, Grace, Irina Gaynanova, and Christian L. Müller. "Microbial networks in SPRING-Semi-parametric rank-based correlation and partial correlation estimation for quantitative microbiome data." Frontiers in Genetics 10 (2019).
#' @examples
#' set.seed(1)
#' x <- c(rep(0,10),rpois(10,100))
#' y <- c(rep(0,5),rpois(15,10))
#' z <- rbind(x,y)
#' z_mclr <- mclr(z)
#' # compare with clr
#' z_clr <- apply(z,1,compositions::clr)
#' 
#' @export
mclr <- function(Z, base=exp(1), eps=0.1){
  index <- which(Z>0 & !is.na(Z))
  a <- apply(Z,1,function(x){
    y <- x[x > 0 & !is.na(x)]
    y <- log(y / gm_mean_pos(y), base)
    x[x > 0 & !is.na(x)] <- y
    x
  })
  a <- t(a)
  Z[index] <- a[index] + abs(min(a)) + eps

  return(Z)
}

#' @param x A vector of sequencing counts from one observation
#' @return The geometric mean for positive components only
gm_mean_pos <- function(x){
  y <- x[x > 0 & !is.na(x)]
  if (length(y)==0){
    stop("The input vector is all zero!")
  } else {
    exp(sum(log(y)) / length(y))
  }
}



