#' Cross-temporal reconciliation tools
#'
#' @description
#' \loadmathjax
#' Some useful tools for the cross-temporal forecast reconciliation of a linearly constrained
#' (hierarchical/grouped) multiple time series.
#'
#' @param m Highest available sampling frequency per seasonal cycle (max. order
#' of temporal aggregation, \mjseqn{m}), or a subset of the \mjseqn{p} factors
#' of \mjseqn{m}.
#' @param h Forecast horizon for the lowest frequency (most temporally aggregated) time
#' series (\emph{default} is \code{1}).
#' @param C (\mjseqn{n_a \times n_b}) cross-sectional (contemporaneous) matrix
#' mapping the bottom level series into the higher level ones.
#' @param Ut Zero constraints cross-sectional (contemporaneous) kernel matrix
#' \mjseqn{(\mathbf{U}'\mathbf{y} = \mathbf{0})} spanning the null space valid
#' for the reconciled forecasts. It can be used instead of parameter
#' \code{C}, but \code{nb} (\mjseqn{n = n_a + n_b}) is needed if
#' \mjseqn{\mathbf{U}' \neq [\mathbf{I} \ -\mathbf{C}]}{}. If the hierarchy
#' admits a structural representation, \mjseqn{\mathbf{U}'} has dimension
#' (\mjseqn{n_a \times n}).
#' @param nb Number of bottom time series; if \code{C} is present, \code{nb}
#' and \code{Ut} are not used.
#' @param sparse Option to return sparse object (\emph{default} is \code{TRUE}).
#'
#' @return
#' \strong{ctf} list with:
#' \item{\code{Ht}}{Full row-rank cross-temporal zero constraints (kernel)
#' matrix coherent with \mjseqn{\mathbf{y} = \mbox{vec}(\mathbf{Y}')}: \mjseqn{\mathbf{H}'\mathbf{y} = \mathbf{0}}.}
#' \item{\code{Hbrevet}}{Complete, not full row-rank cross-temporal zero
#' constraints (kernel) matrix coherent with \mjseqn{\mathbf{y} = \mbox{vec}(\mathbf{Y}')}:
#' \mjseqn{\breve{\mathbf{H}}'\mathbf{y} = \mathbf{0}}.}
#' \item{\code{Hcheckt}}{Full row-rank cross-temporal zero constraints (kernel) matrix coherent with
#' \mjseqn{\check{\mathbf{y}}} (structural representation):
#' \mjseqn{\check{\mathbf{H}}' \check{\mathbf{y}} = \mathbf{0}}.}
#' \item{\code{Ccheck}}{Cross-temporal aggregation matrix \mjseqn{\check{\mathbf{C}}}
#' coherent with \mjseqn{\check{\mathbf{y}}} (structural representation).}
#' \item{\code{Scheck}}{Cross-temporal summing matrix \mjseqn{\check{\mathbf{S}}}
#' coherent with \mjseqn{\check{\mathbf{y}}} (structural representation).}
#' \item{\code{Stilde}}{Cross-temporal summing matrix \mjseqn{\widetilde{\mathbf{S}}}
#' coherent with \mjseqn{\mathbf{y} = \mbox{vec}(\mathbf{Y}')}.}
#'
#' \strong{hts} list from \code{\link{hts_tools}} .
#'
#' \strong{thf} list from \code{\link{thf_tools}} .
#'
#' @examples
#' # One level hierarchy (na = 1, nb = 2) with quarterly data
#' obj <- ctf_tools(C = matrix(c(1, 1), 1), m = 4)
#'
#' @keywords utilities
#' @family utilities
#'
#' @import Matrix
#'
#' @export
ctf_tools <- function(C, m, h = 1, Ut, nb, sparse = TRUE) {
  # Using Ut or C
  if (missing(C)) {
    if (missing(Ut)) {
      stop("Please, give C or Ut.", call. = FALSE)
    } else if(missing(nb)){
      hts <- hts_tools(Ut = Ut, h = 1, sparse = TRUE)
    } else {
      hts <- hts_tools(Ut = Ut, nb = nb, h = 1, sparse = TRUE)
    }
  } else {
    hts <- hts_tools(C = C, h = 1, sparse = TRUE)
  }

  n <- hts$n
  na <- hts$na
  nb <- hts$nb
  C <- hts$C
  S <- hts$S
  Ut <- hts$Ut

  tmp <- thf_tools(m, h = h, sparse = TRUE)
  m <- tmp$m
  Zt <- tmp$Zt
  K <- tmp$K
  R <- tmp$R
  ks <- tmp$ks
  kt <- tmp$kt

  Hbrevet <- rbind(kronecker(Ut, .sparseDiagonal(h * kt)), kronecker(.sparseDiagonal(n), Zt))
  P <- commat(h * kt, n)
  Us <- cbind(Matrix(0, h * NROW(Ut) * m, n * h * ks), kronecker(.sparseDiagonal(h * m), Ut)) %*% P
  Ht <- rbind(Us, kronecker(.sparseDiagonal(n), Zt))

  out1 <- list()
  if(!sparse){
    out1$Ht <- as.matrix(Ht)
    out1$Hbrevet <- as.matrix(Hbrevet)
  }else{
    out1$Ht <- Ht
    out1$Hbrevet <- Hbrevet
  }

  if(!is.null(C)){
    Ccheck <- rbind(kronecker(C, R),kronecker(.sparseDiagonal(nb), K))
    Hcheckt <- cbind(.sparseDiagonal(h*(na*m + n*ks)), -Ccheck)
    Scheck <- rbind(Ccheck, .sparseDiagonal(nb*m*h))
    Stilde <- kronecker(S, R)

    if(!sparse){
      out1$Hcheckt <- as.matrix(Hcheckt)
      out1$Ccheck <- as.matrix(Ccheck)
      out1$Scheck <- as.matrix(Scheck)
      out1$Stilde <- as.matrix(Stilde)
    }else{
      out1$Hcheckt <- Hcheckt
      out1$Ccheck <- Ccheck
      out1$Scheck <- Scheck
      out1$Stilde <- Stilde
    }
  }

  return(list(ctf = out1, hts = hts, thf = tmp))
}
