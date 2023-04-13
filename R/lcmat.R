#' Linear Combination Matrix for a general linearly constrained multiple time series
#'
#' @description
#' \loadmathjax
#' When working with a general linearly constrained multiple (\mjseqn{n}-variate) time series
#' (\mjseqn{\mathbf{x}_t}), getting a linear combination matrix
#' \mjseqn{\bar{\mathbf{C}}} is a critical step to obtain a \emph{structural-like}
#' representation such that, for \mjseqn{t = 1, ..., T},
#' \mjsdeqn{\bar{\mathbf{U}}'= [\mathbf{I} \quad -\bar{\mathbf{C}}] \quad \Rightarrow
#' \quad \mathbf{y}_t = \mathbf{P}\mathbf{x}_t = \left[\begin{array}{c}
#' \mathbf{v}_t\cr
#' \mathbf{f}_t
#' \end{array}\right] = \left[\begin{array}{c}
#' \mathbf{\bar{\mathbf{C}}}\cr
#' \mathbf{I}
#' \end{array}\right]\mathbf{f}_t = \mathbf{\bar{\mathbf{S}}}\mathbf{f}_t,}
#' where \mjseqn{\bar{\mathbf{U}}'} is the (\mjseqn{n_v \times n}) full rank zero constraints matrix,
#' \mjseqn{\bar{\mathbf{S}}} is the (\mjseqn{n \times n_f}) matrix analogous of the summing matrix
#' \mjseqn{\mathbf{S}} for a genuine hierarchical/groupped times series,
#' \mjseqn{\bar{\mathbf{C}}} is the (\mjseqn{n_v \times n_f}) linear combination matrix
#' such that \mjseqn{\mathbf{v}_t = \bar{\mathbf{C}}\mathbf{f}_t},
#' \mjseqn{\mathbf{v}_t} is the (\mjseqn{n_v \times 1}) vector of ‘basic’ variables, and
#' \mjseqn{\mathbf{f}_t} is the (\mjseqn{n_f \times 1}) vector of ‘free’ variables
#' (Di Fonzo and Girolimetto, 2022).
#'
#' @param Gt (\mjseqn{r \times n}) coefficient matrix (\mjseqn{\mathbf{\Gamma}'})
#' for a general linearly constrained multiple time series (\mjseqn{\mathbf{x}_t})
#' such that \mjseqn{\mathbf{\Gamma}'\mathbf{x}_t = \mathbf{0}_{(r \times 1)}}.
#' @param alg Technique used to trasform \mjseqn{\mathbf{\Gamma}'} in
#' \mjseqn{\bar{\mathbf{U}}' = [\mathbf{I} \quad -\bar{\mathbf{C}}]}, such that
#' \mjseqn{\bar{\mathbf{U}'}\mathbf{y}_t = \mathbf{0}_{(n_v \times 1)}}. Use
#' \code{"rref"} for the Row Reduced Echelon Form through Gauss-Jordan elimination
#' (\emph{default}), or \code{"qr"} for the (pivoting) QR decomposition (Strang, 2019).
#' @param tol Tolerance for the \code{"rref"} or  \code{"qr"} algorithm.
#' @param verbose If \code{TRUE}, intermediate steps are printed (\emph{default} is \code{FALSE}).
#' @param sparse Option to return a sparse \mjseqn{\bar{\mathbf{C}}}
#' matrix (\emph{default} is \code{TRUE}).
#'
#' @details
#' Looking for an analogous of the summing matrix \mjseqn{\mathbf{S}}, say
#' \mjseqn{\bar{\mathbf{S}} = \left[\begin{array}{c}
#' \mathbf{\bar{\mathbf{C}}}\cr
#' \mathbf{I}
#' \end{array}\right]}, the \code{lcmat} function transforms  \mjseqn{\mathbf{\Gamma}'} into
#' \mjseqn{\bar{\mathbf{U}}' = [\mathbf{I} \quad -\bar{\mathbf{C}}]}, such that
#' \mjseqn{\bar{\mathbf{U}}'\mathbf{y}_t = \mathbf{0}_{(n_v \times 1)}}.
#' Consider the simple example of a linearly constrained multiple time series consisting
#' of two hierarchies, each with distinct bottom time series,
#' with a common top-level series (\mjseqn{X}):
#' \mjsdeqn{\begin{array}{l}
#' 1)\; X = C + D,\cr
#' 2)\; X = A + B, \cr
#' 3)\; A = A1 + A2.
#' \end{array}}
#' The coefficient matrix \mjseqn{\mathbf{\Gamma}'} of the linear system
#' \mjseqn{\mathbf{\Gamma}'\mathbf{x}_t=\mathbf{0}}
#' (\mjseqn{\mathbf{x}_t = [X\; C\; D\; A\; B\; A1\; A2]}) is
#' \mjsdeqn{\mathbf{\Gamma}' = \left[\begin{array}{ccccccc}
#' 1 & -1 & -1 & 0 & 0 & 0 & 0 \cr
#' 1 & 0 & 0 & -1 & -1 & 0 & 0 \cr
#' 0 & 0 & 0 & 1 & 0 & -1 & -1
#' \end{array}\right].}
#' The \link[FoReco]{lcmat} function returns
#' \mjsdeqn{\bar{\mathbf{C}} = \left[\begin{array}{cccc}
#' 0 & 1 & 1 & 1 \cr
#' -1 & 1 & 1 & 1 \cr
#' 0 & 0 & -1 & -1
#' \end{array}\right].}
#' Then
#' \mjsdeqn{\bar{\mathbf{U}}' = \left[\begin{array}{ccc|cccc}
#' 1 & 0 & 0 & 0 & -1 & -1 & -1 \cr
#' 0 & 1 & 0 & 1 & -1 & -1 & -1 \cr
#' 0 & 0 & 1 & 0 & 0 & 1 & 1
#' \end{array}\right], \quad \mbox{with} \quad
#' \bar{\mathbf{U}}'\mathbf{y}_t = \bar{\mathbf{U}}' \left[\begin{array}{c}
#' \mathbf{v}_t \cr
#' \mathbf{f}_t
#' \end{array}\right] = \mathbf{0},}
#' where \mjseqn{\mathbf{v}_t = [X\; C\; A]}, and
#' \mjseqn{\mathbf{f}_t = [D\; B\; A1\; A2]}.
#'
#' @return A list with
#' \item{\code{Cbar}}{(\mjseqn{n_v \times n_f}) linear combination matrix \mjseqn{\bar{\mathbf{C}}}}
#' \item{\code{pivot}}{(\mjseqn{n \times 1}) vector of the column permutations
#' s.t. \mjseqn{\mathbf{P} = \mathbf{I}[,\mbox{pivot}]}}
#'
#' @examples
#' Gt <- matrix(c(1,-1,-1,0,0,0,0,
#'                1,0,0,-1,-1,0,0,
#'                0,0,0,1,0,-1,-1), nrow = 3, byrow = TRUE)
#' Cbar <- lcmat(Gt = Gt)$Cbar
#' P <- diag(1, NCOL(Gt))[,lcmat(Gt = Gt)$pivot]
#'
#' @usage lcmat(Gt, alg = "rref", tol = sqrt(.Machine$double.eps),
#'        verbose = FALSE, sparse = TRUE)
#'
#' @references
#' Di Fonzo, T., Girolimetto, D. (2022), \emph{Point and probabilistic forecast reconciliation
#' for general linearly constrained multiple time series} (mimeo).
#'
#' Strang, G., (2019), \emph{Linear algebra and learning from data}, Wellesley, Cambridge Press.
#'
#' @keywords utilities
#' @family utilities
#'
#' @import Matrix
#'
#' @export
#'
lcmat <- function(Gt, alg = "rref", tol = sqrt(.Machine$double.eps),
                  verbose = FALSE, sparse = TRUE){

  alg <- match.arg(alg, c("rref", "qr"))

  # Check inputs
  if(is.matrix(Gt) | is(Gt, "Matrix")){
    Gt <- as(Gt, "CsparseMatrix")
  }else{
    stop("Gt must be a matrix", call. = FALSE)
  }
  Ut <- Gt
  if(alg == "rref"){
    # Reduced Row Echelon Form
    if(verbose){
      message("RREF step")
    }

    # Info matrix
    cnames <- colnames(Ut)
    n <- nrow(Ut)
    m <- ncol(Ut)

    # Verbose condition
    if(verbose){
      if(min(n, m) < 100){
        quant <- n + m
      }else{
        quant <- round(stats::quantile(1:min(n, m),
                                       probs = seq(0, 1, 0.25)[-1]))
        message("Gauss-Jordan elimination (%): 0% ", appendLF = FALSE)
      }
    }

    xpos <- 1 # col
    ypos <- 1 # row

    # Gauss-Jordan elimination:
    while((xpos <= m) & (ypos <= n)){
      col <- as(Ut[,xpos], "sparseVector")
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
          id <- Ut@i
          lp <- Ut@i
          lp[id == (whc - 1)] <- ypos - 1

          lp[id == (ypos - 1)] <- whc - 1
          Ut <- sparseMatrix(i = lp + 1, p = Ut@p, x = Ut@x)
        }

        Ut[ypos,] <- Ut[ypos,]/pivot # pivot

        row <-as(Ut[ypos,], "sparseVector")
        Ut <- Ut - Ut[,xpos, drop = FALSE] %*% t(row)
        Ut[ypos,] <- row # restore current row
        xpos <- xpos+1
        ypos <- ypos+1
      }

      if(verbose){
        if(min(n, m) == n){
          step <- ypos
        }else{
          step <- xpos
        }

        printstep <- (step%%quant == 0)
        if(any(printstep)){
          message(names(quant[printstep]), appendLF = FALSE)
          quant <- quant[-1]

          if(length(quant) > 0){
            message(" ", appendLF = FALSE)
          }
        }

      }
    }

    # Generalized Reduced Zero constraints cross-sectional kernel matrix
    Ut <- drop0(Ut, tol = tol)
    Ut <- Ut[rowSums(abs(Ut))!=0,, drop = FALSE]
    colnames(Ut) <- cnames

    # Generalized cross-sectional combination matrix
    pivot <- apply(Ut, 1, function(x) which(x == 1)[1])
    Cbar <- -Ut[, setdiff(1:NCOL(Ut), pivot), drop = FALSE]
    pivot <- c(pivot, setdiff(1:NCOL(Ut), pivot))
    rownames(Cbar) <- cnames[pivot[1:NROW(Cbar)]]
    if(!sparse){
      Cbar <- as.matrix(Cbar)
    }

    # Ite return
    return(list(Cbar = Cbar, pivot = pivot))
  }else{
    # QR decomposition
    if(verbose){
      message("QR decomposition step")
    }

    QR <- base::qr(Ut, tol = tol)
    id_qr <- order(QR$pivot)
    R <- qr.R(QR) # Row Echelon Form

    # Zeros rows
    indep_rows <- as.logical(apply(drop0(R, tol = tol)!=0, 1, max))
    # First not zero number by row
    i_dep <- apply(drop0(R, tol = tol)!=0, 1, which.max)

    # basic variable
    i_dep <- sort(i_dep[indep_rows])
    # free variable
    i_indep <- setdiff(1:NCOL(Ut), i_dep)

    # R = [R1 | R2]
    R1 <- Matrix(R[indep_rows, i_dep, drop = FALSE], sparse = TRUE)
    R2 <- Matrix(R[indep_rows, i_indep, drop = FALSE], sparse = TRUE)

    # Generalized cross-sectional combination matrix
    Cbar <- -solve(R1, R2)
    Cbar <- Cbar[, id_qr[id_qr>NROW(Cbar)]-NROW(Cbar), drop = FALSE]

    # Adjust output
    Cbar <- drop0(Cbar, tol = tol)
    if(!sparse){
      Cbar <- as.matrix(Cbar)
    }

    # QR return
    return(list(Cbar = Cbar,
                pivot = c(QR$pivot[1:NROW(Cbar)],
                          sort(QR$pivot[-c(1:NROW(Cbar))]))))
  }
}
