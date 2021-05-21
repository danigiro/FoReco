#' @title Shrinkage of the covariance matrix
#'
#' @description Shrinkage of the covariance matrix according to \enc{Schäfer}{Schafer} and Strimmer (2005).
#'
#' @param x residual matrix
#' @param minT this param allows to calculate the covariance matrix according
#' to the original hts formulation (\code{TRUE}) or according to the standard
#' approach (\code{FALSE}).
#'
#' @author
#' This function is a modified version of the \code{shrink_estim}() hidden function of \pkg{hts}.
#'
#' @return A list with two objects: the first (\code{$scov}) is the shrunk covariance matrix
#' and the second (\code{$lambda}) is the shrinkage intensity coefficient.
#'
#' @references
#' \enc{Schäfer}{Schafer}, J.L., Strimmer, K. (2005), A Shrinkage Approach to Large-Scale Covariance Matrix
#' Estimation and Implications for Functional Genomics, \emph{Statistical Applications in Genetics
#' and Molecular Biology}, 4, 1
#'
#' Hyndman, R. J., Lee, A., Wang, E., and Wickramasuriya, S. (2020).
#' hts: Hierarchical and Grouped Time Series, \emph{R package version 6.0.1},
#' \href{https://CRAN.R-project.org/package=hts}{https://CRAN.R-project.org/package=hts}.
#'
#' @keywords utilities
#' @family utilities
#'
#' @export
shrink_estim <- function(x, minT = T) {
  if (is.matrix(x) == TRUE && is.numeric(x) == FALSE) {
    stop("The data matrix must be numeric!", call. = FALSE)
  }

  x <- stats::na.omit(x)
  p1 <- ncol(x)
  n2 <- nrow(x)

  if (minT == T) {
    covm <- crossprod(x) / n2
  } else {
    covm <- stats::cov(x)
  }

  tar <- diag(diag(covm))
  corm <- stats::cov2cor(covm)
  xs <- scale(x, center = FALSE, scale = sqrt(diag(covm)))
  xs <- xs[stats::complete.cases(xs), ]
  v <- (1 / (n2 * (n2 - 1))) * (crossprod(xs^2) - 1 / n2 * (crossprod(xs))^2)
  diag(v) <- 0
  corapn <- stats::cov2cor(tar)
  d <- (corm - corapn)^2
  lambda <- sum(v) / sum(d)
  lambda <- max(min(lambda, 1), 0)
  shrink.cov <- lambda * tar + (1 - lambda) * covm
  return(list(scov = shrink.cov, lambda = lambda))
}
