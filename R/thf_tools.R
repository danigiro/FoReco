#' Temporal reconciliation tools
#'
#' Some useful tools for forecast reconciliation through temporal hierarchies.
#'
#' @param m Highest available sampling frequency per seasonal cycle (max. order of temporal aggregation).
#' @param h Forecast horizon for the lowest frequency (most temporally aggregated) time series (\emph{default} is \code{1}).
#' @param sparse Option to return sparse object (\emph{default} is \code{TRUE}).
#'
#' @return A list of seven elements:
#' \item{\code{K}}{Temporal aggregation matrix.}
#' \item{\code{R}}{Temporal summing matrix.}
#' \item{\code{Zt}}{Zero constraints temporal kernel matrix\eqn{,\; \textbf{Z}_h'\textbf{Y}' =
#' \mathbf{0}_{\left[hk^* \times n \right]}}{}.}
#' \item{\code{kset}}{Set of factors (\code{p}) of \code{m} in descending order (from \code{m}
#' to 1)\eqn{,\;{\cal K} = \left\{k_p, k_{p-1}, \ldots,\right.}{} \eqn{\left. k_2, k_1\right\},}{} \eqn{k_p=m, \; k_1=1}{}.}
#' \item{\code{p}}{Number of elements of kset\eqn{,\;({\cal K})}{}.}
#' \item{\code{ks}}{Sum of \code{p-1} factors of \code{m} (out of \code{m} itself), \eqn{k^*}{k*}.}
#' \item{\code{kt}}{Sum of all factors of m (\eqn{k^{tot} = k^*+m}{kt = ks + m}).}
#'
#' @examples
#' # quarterly data
#' obj <- thf_tools(m = 4, sparse = FALSE)
#' @keywords utilities
#' @import Matrix
#' @export
#'
thf_tools <- function(m, h = 1, sparse = TRUE) {
  if (m < 2) {
    stop("m must be > 1", call. = FALSE)
  }

  out <- list()
  kset <- rev(divisors(m))
  p <- length(kset)
  ks <- sum(kset) - m
  kt <- sum(kset)

  if (sparse) {
    out$K <- do.call("rbind", Map(
      function(x, y) kronecker(x, y), lapply(rev(kset[-1]), function(x) Diagonal(x * h)),
      lapply(kset[-p], function(x) t(rep(1, x)))
    ))
    out$Zt <- cbind(Diagonal(h * ks), -out$K)
    out$R <- rbind(out$K, Diagonal(m * h))
  } else {
    out$K <- do.call("rbind", Map(
      function(x, y) kronecker(x, y), lapply(rev(kset[-1]), function(x) diag(x * h)),
      lapply(kset[-p], function(x) t(rep(1, x)))
    ))
    out$Zt <- cbind(diag(h * ks), -out$K)
    out$R <- rbind(out$K, diag(m * h))
  }

  out$kset <- kset
  out$p <- p
  out$ks <- ks
  out$kt <- kt
  return(out)
}

# x is a int number
# return: all factors of x
divisors <- function(x) {
  x <- as.integer(x)
  div <- seq_len(abs(x))
  factors <- div[x %% div == 0L]
  return(factors)
}
