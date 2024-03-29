#' Temporal reconciliation tools
#'
#' @description
#' \loadmathjax
#' Some useful tools for forecast reconciliation through temporal hierarchies.
#'
#' @param m Highest available sampling frequency per seasonal cycle (max. order of temporal aggregation, \mjseqn{m}),
#' or a subset of the \mjseqn{p} factors of \mjseqn{m}.
#' @param h Forecast horizon for the lowest frequency (most temporally aggregated) time series (\emph{default} is \mjseqn{1}).
#' @param sparse Option to return sparse object (\emph{default} is \code{TRUE}).
#'
#' @return A list of seven elements:
#' \item{\code{K}}{Temporal aggregation matrix.}
#' \item{\code{R}}{Temporal summing matrix.}
#' \item{\code{Zt}}{Zero constraints temporal kernel matrix, \mjseqn{\mathbf{Z}_h'\mathbf{Y}' =
#' \mathbf{0}_{\left[hk^* \times n \right]}}.}
#' \item{\code{kset}}{Set of factors (\mjseqn{p}) of \mjseqn{m} in descending order (from \mjseqn{m}
#' to 1), \mjseqn{{\cal K} = \left\lbrace k_p, k_{p-1}, \ldots, k_2, k_1\right\rbrace}, \mjseqn{k_p=m}, \mjseqn{k_1=1}.}
#' \item{\code{m}}{Highest available sampling frequency per seasonal cycle (max. order of temporal aggregation).}
#' \item{\code{p}}{Number of elements of kset, \mjseqn{{\cal K}}.}
#' \item{\code{ks}}{Sum of \mjseqn{p-1} factors of \mjseqn{m} (out of \mjseqn{m} itself), \mjseqn{k^*}.}
#' \item{\code{kt}}{Sum of all factors of m, \mjseqn{k^{tot} = k^*+m}.}
#'
#' @examples
#' # quarterly data
#' obj <- thf_tools(m = 4, sparse = FALSE)
#'
#' @keywords utilities
#' @family utilities
#'
#' @import Matrix
#' @export
#'
thf_tools <- function(m, h = 1, sparse = TRUE) {

  out <- list()
  if(length(m)>1){
    kset <- sort(m, decreasing = TRUE)
    m <- max(kset)

    if(m < 2){
      stop("m must be > 1", call. = FALSE)
    }

    if(min(kset)!=1){
      kset <- sort(c(kset, 1), decreasing = TRUE)
    }

    if(!all(kset %in% rev(divisors(m)))){
      stop("(", paste(kset, collapse = ", "), ") is not a subset of (", paste(rev(divisors(m)), collapse = ", "), ")", call. = FALSE)
    }
  }else{
    if(m < 2){
      stop("m must be > 1", call. = FALSE)
    }
    kset <- rev(divisors(m))
  }

  p <- length(kset)
  ks <- sum(m/kset[-p])
  kt <- sum(m/kset)

  if (sparse) {
    out$K <- do.call("rbind", Map(
      function(x, y) kronecker(x, y), lapply(m/kset[-p], function(x) .sparseDiagonal(x * h)),
      lapply(kset[-p], function(x) t(rep(1, x)))
    ))
    out$Zt <- cbind(.sparseDiagonal(h * ks), -out$K)
    out$R <- rbind(out$K, .sparseDiagonal(m * h))
  } else {
    out$K <- do.call("rbind", Map(
      function(x, y) kronecker(x, y), lapply(rev(kset[-1]), function(x) diag(x * h)),
      lapply(kset[-p], function(x) t(rep(1, x)))
    ))
    out$Zt <- cbind(diag(h * ks), -out$K)
    out$R <- rbind(out$K, diag(m * h))
  }

  out$kset <- kset
  out$m <- m
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
