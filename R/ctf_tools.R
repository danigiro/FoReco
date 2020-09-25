#' Cross-temporal reconciliation tools
#'
#' Some useful tools for the cross-temporal forecast reconciliation of cross-sectionally and
#' temporally constrained time series
#'
#' @param C (\code{na x nb}) cross-sectional (contemporaneous) matrix mapping the bottom
#' level series into the higher level ones.
#' @param m Highest available sampling frequency per seasonal cycle (max. order of
#' temporal aggregation).
#' @param h Forecast horizon for the lowest frequency (most temporally aggregated) time
#' series (\emph{default} is \code{1}).
#' @param Ut Zero constraints cross-sectional (contemporaneous) kernel matrix
#' \eqn{(\textbf{U}'\textbf{Y} = \mathbf{0})}{} spanning the null space valid for the reconciled
#' forecasts. It can be used instead of parameter \code{C}, but in this case \code{nb} (n = na + nb) is needed. If
#' the hierarchy admits a structural representation, \code{Ut} has dimension (\code{na x n}).
#' @param nb Number of bottom time series; if \code{C} is present, \code{nb} is not used.
#' @param Sstruc If \code{Sstruc = TRUE} the function returns also the structural representation matrix of
#' a cross-temporal system (\emph{default} is \code{FALSE}).
#'
#'
#' @return
#' \strong{ctf} list with:
#' \item{\code{Ht}}{Full row-rank cross-temporal zero-valued constraints (kernel)
#' matrix\eqn{,\; \textbf{H}'\textbf{y} = \mathbf{0}}{}.}
#' \item{\code{Htbreve}}{Complete, not full row-rank cross-temporal zero-valued
#' constraints (kernel) matrix.}
#' \item{\code{Htstruc}}{Zero constraints full row-rank cross-temporal kernel matrix
#' (structural representation) \eqn{,\; \check{\textbf{H}}'}{}.}
#' \item{\code{Sstruc}}{Cross-temporal summing matrix (structural
#' representation)\eqn{,\; \check{\textbf{S}}}{}.}
#'
#' \strong{hts} list from \code{\link{hts_tools}}
#'
#' \strong{thf} list from \code{\link{thf_tools}}
#'
#' @examples
#' # One level hierarchy (na = 1, nb = 2) with quarterly data
#' obj <- ctf_tools(C = matrix(c(1, 1), 1), m = 4, Sstruc = TRUE)
#' @keywords utilities
#' @import Matrix
#'
#' @export
#'
ctf_tools <- function(C, m, h = 1, Ut, nb, Sstruc = FALSE) {
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

  tmp <- thf_tools(m, h = h, sparse = TRUE)
  Zt <- tmp$Zt
  ks <- tmp$ks
  kt <- tmp$kt

  Htbreve <- rbind(kronecker(Ut, Diagonal(h * kt)), kronecker(Diagonal(n), Zt))
  D <- commat(h * kt, n)
  Us <- cbind(Matrix(0, h * NROW(Ut) * m, n * h * ks), kronecker(Diagonal(h * m), Ut)) %*% D
  Ht <- rbind(Us, kronecker(Diagonal(n), Zt))

  out1 <- list()
  out1$Ht <- Ht
  out1$Htbreve <- Htbreve

  if (Sstruc) {
    if (missing(C)) {
      stop("The argument C is not specified, but it's necessary for the structural representation of a cross-temporal system")
    }
    Qtilde <- commat(nb, kt) %*% bdiag(commat(ks, nb), commat(m, nb))
    Q <- bdiag(Diagonal(na * kt), Qtilde)

    Htstruc <- Matrix(pracma::rref(as.matrix(Ht %*% Q)), sparse = TRUE)
    r <- nrow(Htstruc)
    c <- ncol(Htstruc)
    Cstruc <- -Htstruc[, (r + 1):c]
    Sstruc <- rbind(Cstruc, Diagonal(m * nb))

    # Trivial checks on matrix Hstruc
    # 1. The left square block must be an identity matrix
    # 2. The last m*nb elements of the first row mast be all equal to -1
    flag <- 0
    Htcheck1 <- all(abs(Htstruc[, 1:r] - Diagonal(r)) < 1e-6)
    Htcheck2 <- all(abs(Htstruc[1, (r + 1):c] + rep(1, m * nb)) < 1e-6)
    if (!Htcheck1 | !Htcheck2) {
      flag <- 1
      warning("Failed checks on matrix Htstruc. Check the results!", call. = FALSE)
    }

    out1$Htstruc <- Htstruc
    out1$Sstruc <- Sstruc
  }


  out2 <- list()
  if (!missing(C)) {
    out2$S <- kronecker(S, Diagonal(h))
  }
  out2$Ut <- kronecker(Ut, Diagonal(h))
  out2$n <- n
  out2$na <- na
  out2$nb <- nb

  return(list(ctf = out1, hts = out2, thf = tmp))
}
