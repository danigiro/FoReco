#' Linear combination (aggregation) matrix for a general linearly constrained multiple time series
#'
#' @description
#' This function transforms a general (possibly redundant) zero constraints matrix into a
#' linear combination (aggregation) matrix \eqn{\mathbf{A}_{cs}}.
#' When working with a general linearly constrained multiple (\eqn{n}-variate)
#' time series, getting a linear combination matrix \eqn{\mathbf{A}_{cs}} is a critical
#' step to obtain a structural-like representation such that
#' \deqn{\mathbf{C}_{cs} = [\mathbf{I} \quad -\mathbf{A}],}
#' where \eqn{\mathbf{C}_{cs}} is the full rank zero constraints matrix (Girolimetto and
#' Di Fonzo, 2023).
#'
#' @param cons_mat A (\eqn{r \times n}) numeric matrix representing the cross-sectional
#' zero constraints.
#' @param method  Method to use: "\code{rref}" for the Reduced Row Echelon
#' Form through Gauss-Jordan elimination (\emph{default}), or "\code{qr}"
#' for the (pivoting) QR decomposition (Strang, 2019).
#' @param tol Tolerance for the "\code{rref}" or "\code{qr}" method.
#' @param verbose If \code{TRUE}, intermediate steps are printed (\emph{default} is \code{FALSE}).
#' @param sparse Option to return a sparse matrix (\emph{default} is \code{TRUE}).
#'
#' @returns A list with two elements: (i) the linear combination (aggregation) matrix
#' (\code{agg_mat}) and (ii) the vector of the column permutations (\code{pivot}).
#'
#' @examples
#' ## Two hierarchy sharing the same top-level variable, but not sharing the bottom variables
#' #        X            X
#' #    |-------|    |-------|
#' #    A       B    C       D
#' #  |---|
#' # A1   A2
#' # 1) X = C + D,
#' # 2) X = A + B,
#' # 3) A = A1 + A2.
#' cons_mat <- matrix(c(1,-1,-1,0,0,0,0,
#'                1,0,0,-1,-1,0,0,
#'                0,0,0,1,0,-1,-1), nrow = 3, byrow = TRUE)
#' obj <- lcmat(cons_mat = cons_mat, verbose = TRUE)
#' agg_mat <- obj$agg_mat # linear combination matrix
#' pivot <- obj$pivot # Pivot vector
#'
#' @usage lcmat(cons_mat, method = "rref", tol = sqrt(.Machine$double.eps),
#'        verbose = FALSE, sparse = TRUE)
#'
#' @references
#' Girolimetto, D. and Di Fonzo, T. (2023), Point and probabilistic forecast reconciliation
#' for general linearly constrained multiple time series,
#' \emph{Statistical Methods & Applications}, in press. \doi{10.1007/s10260-023-00738-6}.
#'
#' Strang, G. (2019), \emph{Linear algebra and learning from data}, Wellesley, Cambridge Press.
#'
#' @family Utilities
#'
#' @export
#'
lcmat <- function(cons_mat, method = "rref", tol = sqrt(.Machine$double.eps),
                  verbose = FALSE, sparse = TRUE){

  method <- match.arg(method, c("rref", "qr"))

  # Check inputs
  if(is.matrix(cons_mat) | is(cons_mat, "Matrix")){
    cons_mat <- as(cons_mat, "CsparseMatrix")
  }else{
    cli_abort("{.arg cons_mat} is not a numeric matrix.", call = NULL)
  }
  cons_mat <- cons_mat
  if(method == "rref"){ # Reduced Row Echelon Form
    # Info matrix
    cnames <- colnames(cons_mat)
    n <- nrow(cons_mat)
    m <- ncol(cons_mat)

    # Verbose condition
    prog <- FALSE
    if(verbose){
      if(min(n, m) < 100){
        prog <- FALSE
      }else{
        prog <- TRUE
        message("Gauss-Jordan elimination (%): |", appendLF = FALSE)

        if(min(n, m)==n){
          stepcheck <- round(seq(2, n, length.out = 10))
        }else{
          stepcheck <- round(seq(2, m, length.out = 10))
        }
      }
    }

    xpos <- 1 # col
    ypos <- 1 # row

    # Gauss-Jordan elimination:
    while((xpos <= m) & (ypos <= n)){
      col <- as(cons_mat[,xpos], "sparseVector")
      col[1:n < ypos] <- 0

      if(sum(abs(col)) < tol){
        pivot <- 0
      }else{ # find maximum pivot in current column at or below current row
        whc <- col@i[which.max(abs(col@x))]
        pivot <- col[whc]
      }

      if(abs(pivot) <= tol){
        xpos <- xpos+1     # check for 0 pivot
      }else{
        if(whc > ypos){  # exchange rows
          id <- cons_mat@i
          lp <- cons_mat@i
          lp[id == (whc - 1)] <- ypos - 1

          lp[id == (ypos - 1)] <- whc - 1
          cons_mat <- sparseMatrix(i = lp + 1, p = cons_mat@p, x = cons_mat@x)
        }

        cons_mat[ypos,] <- cons_mat[ypos,]/pivot # pivot

        row <-as(cons_mat[ypos,], "sparseVector")
        cons_mat <- cons_mat - cons_mat[,xpos, drop = FALSE] %*% t(row)
        cons_mat[ypos,] <- row # restore current row
        xpos <- xpos+1
        ypos <- ypos+1
      }

      if(verbose & prog){
        if(min(n, m) == n){
          step <- ypos
        }else{
          step <- xpos
        }
        if(step%in% stepcheck){
          message("=", appendLF = FALSE)
        }
      }
    }
    if(verbose & prog){
      message("| done\n", appendLF = FALSE)
    }
    # Generalized Reduced Zero constraints cross-sectional kernel matrix
    cons_mat <- drop0(cons_mat, tol = tol)
    cons_mat <- cons_mat[rowSums(abs(cons_mat))!=0,, drop = FALSE]
    colnames(cons_mat) <- cnames

    # Generalized cross-sectional combination matrix
    pivot <- apply(cons_mat, 1, function(x) which(x == 1)[1])
    agg_mat <- -cons_mat[, setdiff(1:NCOL(cons_mat), pivot), drop = FALSE]
    pivot <- c(pivot, setdiff(1:NCOL(cons_mat), pivot))
    rownames(agg_mat) <- cnames[pivot[1:NROW(agg_mat)]]
    agg_mat <- sparse2dense(agg_mat, sparse = sparse)

  }else{ # QR decomposition
    QR <- base::qr(cons_mat, tol = tol)
    id_qr <- order(QR$pivot)
    R <- qr.R(QR) # Row Echelon Form

    # Zeros rows
    indep_rows <- as.logical(apply(drop0(R, tol = tol)!=0, 1, max))
    # First not zero number by row
    i_dep <- apply(drop0(R, tol = tol)!=0, 1, which.max)

    # basic variable
    i_dep <- sort(i_dep[indep_rows])
    # free variable
    i_indep <- setdiff(1:NCOL(cons_mat), i_dep)

    # R = [R1 | R2]
    R1 <- Matrix(R[indep_rows, i_dep, drop = FALSE], sparse = TRUE)
    R2 <- Matrix(R[indep_rows, i_indep, drop = FALSE], sparse = TRUE)

    # Generalized cross-sectional combination matrix
    agg_mat <- -solve(R1, R2)
    agg_mat <- agg_mat[, id_qr[id_qr>NROW(agg_mat)]-NROW(agg_mat), drop = FALSE]

    # Adjust output
    agg_mat <- drop0(agg_mat, tol = tol)
    agg_mat <- sparse2dense(agg_mat, sparse = sparse)

    # QR pivot
    pivot <- c(QR$pivot[1:NROW(agg_mat)],
               sort(QR$pivot[-c(1:NROW(agg_mat))]))
  }
  if(any(pivot != 1:length(pivot)) & verbose){
    cli_bullets(c("!"="A pivot is performed. Remember to apply the pivot also to the base forecast.",
                "i"="E.g. {.code base[, pivot]} in cross-sectional or {.code base[pivot, ]} in cross-temporal."))
  }
  return(list(agg_mat = agg_mat, pivot = pivot))
}
