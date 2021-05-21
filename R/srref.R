#' (Sparse) Reduced Row Echelon Form of a matrix
#'
#' Returns the reduced row echelon form of a (possibly sparse) matrix through Gauss-Jordan
#' elimination. This function is used by the \code{\link{ut2c}} function
#' to obtain a 'structural representation' of a general linearly constrained
#' multiple time series, from its zero-constraints kernel representation
#' (Di Fonzo and Girolimetto, 2020).
#'
#' @param A A numeric (possibly sparse) matrix.
#' @param tol Tolerance for checking for 0 pivot.
#' @param verbose If \code{TRUE}, print intermediate steps.
#' @param sparse Option to return sparse matrices (\emph{default} is \code{TRUE}).
#'
#' @return A (sparse) matrix with the same dimension as \mjseqn{\mathbf{A}}
#'
#' @author Originally written by John Fox and modified by Geoffrey Brent to fix a bug,
#' extended to deal with sparse matrices.
#'
#' @family utilities
#'
#' @references
#' Di Fonzo, T., Girolimetto, D. (2020), Cross-Temporal Forecast Reconciliation:
#' Optimal Combination Method and Heuristic Alternatives, Department of Statistical
#' Sciences, University of Padua, \href{https://arxiv.org/abs/2006.08570}{arXiv:2006.08570}.
#'
#' @export
#'
#' @import Matrix
srref <- function(A, tol=sqrt(.Machine$double.eps), verbose=FALSE, sparse = TRUE){
  if(is.matrix(A) | is(A, "Matrix")){
    A <- as(A, "dgCMatrix")
  }else{
    stop("A must be a matrix", call. = FALSE)
  }

  n <- nrow(A)
  m <- ncol(A)

  # Verbose condition
  if(verbose){
    quant <- round(stats::quantile(1:min(n,m), probs = seq(0, 1, 0.25)[-1]))
    if(min(n,m)<100){
      quant <- n + m
    }else{
      message("Reduced Row Echelon Form (%): 0%...", appendLF = FALSE)
    }
  }

  x.position <- 1 # col
  y.position <- 1 # row
  # change loop:
  while((x.position<=m) & (y.position<=n)){
    col <- as(A[,x.position], "sparseVector")
    col[1:n < y.position] <- 0
    if(sum(abs(col)) < tol){
      pivot <- 0
    }else{ # find maximum pivot in current column at or below current row
      which <- col@i[which.max(abs(col@x))]
      pivot <- col[which]
    }
    if(abs(pivot) <= tol){
      x.position <- x.position+1     # check for 0 pivot
    }else{
      if(which > y.position){  # exchange rows
        id <- A@i
        lp <- A@i
        lp[id==which-1] <- y.position-1

        lp[id==y.position-1] <- which-1
        A <- sparseMatrix(i = lp+1, p=A@p, x=A@x)
      }

      A[y.position,]<-A[y.position,]/pivot # pivot

      row <-as(A[y.position,], "sparseVector")
      A <- A - A[,x.position, drop = FALSE] %*% t(row) # sweep
      A[y.position,]<-row # restore current row
      x.position<-x.position+1
      y.position<-y.position+1
    }


    if(verbose){
      if(min(n,m)==n){
        step <- y.position
      }else{
        step <- x.position
      }
      printstep <- (step%%quant == 0)
      if(any(printstep)){
        message(names(quant[printstep]), appendLF = FALSE)
        quant <- quant[-1]
        if(length(quant) > 0){
          message("...", appendLF = FALSE)
        }
      }

    }
  }

  A <- drop0(A, tol = tol)

  if(!sparse){
    A <- as.matrix(A)
  }

  return(A)
}
