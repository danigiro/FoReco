#' Cross-sectional reconciliation tools
#'
#' @description
#' \loadmathjax
#' Some useful tools for the cross-sectional forecast reconciliation of a
#' linearly constrained (e.g., hierarchical/grouped) multiple time series.
#'
#' @param C (\mjseqn{n_a \times n_b}) cross-sectional (contemporaneous) matrix
#' mapping the bottom level series into the higher level ones.
#' @param Ut Zero constraints cross-sectional (contemporaneous) kernel matrix
#' \mjseqn{(\mathbf{U}'\mathbf{y} = \mathbf{0})} spanning the null space valid
#' for the reconciled forecasts. It can be used instead of parameter
#' \code{C}, but \code{nb} is needed if
#' \mjseqn{\mathbf{U}' \neq [\mathbf{I} \ -\mathbf{C}]}. If the hierarchy
#' admits a structural representation, \mjseqn{\mathbf{U}'} has dimension
#' (\mjseqn{n_a \times n}).
#' @param nb Number of bottom time series; if \code{C} is present, \code{nb}
#' and \code{Ut} are not used.
#' @param h Forecast horizon (\emph{default} is \code{1}).
#' @param sparse Option to return sparse matrices (\emph{default} is \code{TRUE}).
#'
#' @return A list of five elements:
#' \item{\code{C}}{(\mjseqn{n \times n_b}) cross-sectional (contemporaneous) aggregation matrix.}
#' \item{\code{S}}{(\mjseqn{n \times n_b}) cross-sectional (contemporaneous) summing matrix,
#' \mjseqn{\mathbf{S} = \left[\begin{array}{c} \mathbf{C} \cr \mathbf{I}_{n_b}\end{array}\right].}}
#' \item{\code{Ut}}{(\mjseqn{n_a \times n}) zero constraints cross-sectional (contemporaneous)
#' kernel matrix. If the hierarchy admits a structural representation \mjseqn{\mathbf{U}' = [\mathbf{I} \ -\mathbf{C}]}}
#' \item{\code{n}}{Number of variables \mjseqn{n_a + n_b}.}
#' \item{\code{na}}{Number of upper level variables.}
#' \item{\code{nb}}{Number of bottom level variables.}
#'
#' @examples
#' # One level hierarchy (na = 1, nb = 2)
#' obj <- hts_tools(C = matrix(c(1, 1), 1), sparse = FALSE)
#'
#' @keywords utilities
#' @family utilities
#'
#' @import Matrix
#'
#' @export
#'
hts_tools <- function(C, h = 1, Ut, nb, sparse = TRUE) {
  if (missing(C)) {
    if (missing(Ut)) {
      stop("Please, give C or Ut.", call. = FALSE)
    } else {
      if (!(is.matrix(Ut) | is(Ut, "Matrix"))){
        stop("Ut must be a matrix.", call. = FALSE)
      }

      if(!is(Ut, "Matrix")){
        Ut <- Matrix(Ut, sparse = TRUE)
      }

      if(all(.sparseDiagonal(NROW(Ut))==Ut[,1:NROW(Ut)])){
        C <- -Ut[, -c(1:NROW(Ut))]
        if(any(rowSums(abs(C)) == 0)){
          message("Removed a zeros row in C matrix")
          C <- C[rowSums(abs(C)) != 0,]
          Ut <- cbind(.sparseDiagonal(NROW(C)), -C)
        }
        n <- NCOL(Ut)
        na <- NROW(Ut)
        nb <- n - na
        S <- rbind(C, .sparseDiagonal(nb))
      } else if(missing(nb)){
        stop("Ut is not in form [I -C], give also nb", call. = FALSE)
      }else{
        n <- NCOL(Ut)
        na <- n-nb
        C <- NULL
        S <- NULL
      }
    }
  } else {
    if (!(is.matrix(C) | is(C, "Matrix"))){
      stop("C must be a matrix.", call. = FALSE)
    }
    if(!is(C, "Matrix")){
      C <- Matrix(C, sparse = TRUE)
    }

    if(any(rowSums(abs(C)) == 0)){
      message("Removed a zeros row in C matrix")
      C <- C[rowSums(abs(C)) != 0,]
    }

    nb <- ncol(C)
    na <- nrow(C)
    n <- nb + na
    S <- rbind(C, .sparseDiagonal(nb))
    Ut <- cbind(.sparseDiagonal(na), -C)
  }

  if (n <= nb) {
    stop("n <= nb, total number of TS is less (or equal) than the number of bottom TS.", call. = FALSE)
  }

  out2 <- list()
  if (!is.null(C)) {
    if(any(rowSums(abs(C)) == 1)){
      message(paste0("There is only one non-zero value in ", sum(rowSums(abs(C)) == 1),
                     " row(s) of C.",
                     "\nRemember that Foreco can also work with unbalanced hierarchies (recommended)."))
    }

    if(h>1){
      out2$C <- kronecker(C, .sparseDiagonal(h))
      out2$S <- kronecker(S, .sparseDiagonal(h))
    }else{
      out2$S <- S
      out2$C <- C
    }

    if (!sparse) {
      out2$C <- as.matrix(out2$C)
      out2$S <- as.matrix(out2$S)
    }
  }
  out2$Ut <- kronecker(Ut, .sparseDiagonal(h))

  if (!sparse) {
    out2$Ut <- as.matrix(out2$Ut)
  }

  out2$n <- n
  out2$na <- na
  out2$nb <- nb
  return(out2)
}
