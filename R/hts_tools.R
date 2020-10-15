#' Cross-sectional reconciliation tools
#'
#' Some useful tools for the cross-sectional reconciliation of linearly and hierarchically constrained time series
#'
#' @param C (\code{na x nb}) cross-sectional (contemporaneous) matrix mapping the bottom
#' level series into the higher level ones.
#' @param h Forecast horizon (\emph{default} is \code{1}).
#' @param Ut (\code{na x n}) zero constraints cross-sectional (contemporaneous) kernel
#' matrix \eqn{\textbf{U}'\textbf{Y} = \mathbf{0}_{\left[n_a \times (k^*+m)\right]}}{}
#' spanning the null space valid for the reconciled forecasts. It can be used instead of
#' parameter \code{C}, but needs \code{nb} (n = na + nb).
#' @param nb Number of bottom time series; if \code{C} is present, \code{nb} is not used.
#' @param sparse Option to return sparse object (\emph{default} is \code{TRUE}).
#'
#' @return A list of five elements:
#' \item{\code{S}}{(\code{n x nb}) cross-sectional (contemporaneous) summing matrix.}
#' \item{\code{Ut}}{(\code{na x n}) zero constraints cross-sectional (contemporaneous)
#' kernel matrix.}
#' \item{\code{n}}{Number of variables \code{na + nb}.}
#' \item{\code{na}}{Number of upper level variables.}
#' \item{\code{nb}}{Number of bottom level variables.}
#'
#' @examples
#' # One level hierarchy (na = 1, nb = 2)
#' obj <- hts_tools(C = matrix(c(1, 1), 1), sparse = FALSE)
#' @keywords utilities
#' @import Matrix
#'
#' @export
#'
hts_tools <- function(C, h = 1, Ut, nb, sparse = TRUE) {
  if (missing(C)) {
    n <- ncol(Ut)
    na <- n - nb
  } else {
    nb <- ncol(C)
    na <- nrow(C)
    n <- nb + na
    S <- rbind(C, Diagonal(nb))
    Ut <- cbind(Diagonal(na), -C)
  }

  out2 <- list()
  if (!missing(C)) {
    out2$S <- kronecker(S, Diagonal(h))
    if (!sparse) {
      out2$S <- as.matrix(out2$S)
    }
  }
  out2$Ut <- kronecker(Ut, Diagonal(h))

  if (!sparse) {
    out2$Ut <- as.matrix(out2$Ut)
  }

  out2$n <- n
  out2$na <- na
  out2$nb <- nb

  return(out2)
}
