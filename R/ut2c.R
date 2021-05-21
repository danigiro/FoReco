#' Cross-sectional 'structural representation' of a general linearly constrained
#' multiple time series
#'
#' @description
#' \loadmathjax
#' Switching from a zero-constraints kernel representation of a linearly constrained
#' multiple time series, \mjseqn{\mathbf{U}'\mathbf{y}=\mathbf{0}},
#' to a 'structural representation', \mjseqn{\mathbf{y}=\mathbf{S}\mathbf{b}}.
#' (\emph{Experimental version})
#'
#' @param Ut (\mjseqn{r \times c}) zero constraints cross-sectional
#' (contemporaneous) kernel matrix \mjseqn{(\textbf{U}'\textbf{y} = \mathbf{0})}
#' spanning the null space valid for the target forecasts.
#' @param sparse Option to return sparse object (\emph{default} is \code{TRUE}).
#' @param verbose If \code{TRUE}, print intermediate steps of \code{\link{srref}}.
#'
#' @details
#' Consider the simple example of linearly constrained multiple time series consisting
#' of two hierarchies, each with distinct bottom time series,
#' with a common top-level series (\mjseqn{T}):
#' \mjsdeqn{\begin{array}{ll}
#' 1)\; T = A + B        & 4)\; T = X + Y \cr
#' 2)\; A = AA + AB      & 5)\; X = XX + XY \cr
#' 3)\; B = BA + BB + BC & 6)\; Y = YX + YY
#' \end{array}.}
#' Given the cross-sectional aggregation matrices of each hierarchy,
#' \mjsdeqn{\mathbf{C}_1 = \left[\begin{array}{ccccc}
#' 1 & 1 & 0 & 0 & 0\cr
#' 0 & 0 & 1 & 1 & 1
#' \end{array}\right]\quad \mathrm{and} \quad \mathbf{C}_2 = \left[\begin{array}{ccccc}
#' 1 & 1 & 0 & 0\cr
#' 0 & 0 & 1 & 1
#' \end{array}\right],}
#' the zero constraints cross-sectional kernel matrix \mjseqn{\mathbf{U}'},
#' which accounts for the constraints of both hierarchies,
#' can be built as follows:
#' \mjsdeqn{\footnotesize\mathbf{U}' = \left[\begin{array}{cccccccccccccc|c}
#' 1 & 0 & 0 &-1 &-1 &-1 &-1 &-1 & 0 & 0 & 0 & 0 & 0 & 0 & T\cr
#' 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 &-1 &-1 &-1 &-1 & T\cr
#' 0 & 1 & 0 &-1 &-1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & A\cr
#' 0 & 0 & 1 & 0 & 0 &-1 &-1 &-1 & 0 & 0 & 0 & 0 & 0 & 0 & B\cr
#' 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 &-1 &-1 & 0 & 0 & X\cr
#' 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 &-1 &-1 & Y\cr
#' \hline
#' T & A & B & AA& AB& BA& BB& BC& X & Y & XX& XY& YX& YY&
#' \end{array}\right].}
#' Function \code{\link{ut2c}} returns a matrix
#' \mjseqn{\footnotesize\mathbf{U}'_{rref} = \left[\mathbf{I}_6 \; -\mathbf{C} \right]}:
#' \mjsdeqn{\footnotesize\mathbf{U}'_{rref} = \left[\begin{array}{cccccccccccccc|c}
#' 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 &-1 &-1 &-1 &-1 & T\cr
#' 0 & 1 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 1 &-1 &-1 &-1 &-1 & A\cr
#' 0 & 0 & 1 & 0 & 0 & 0 & 0 &-1 &-1 &-1 & 0 & 0 & 0 & 0 & B\cr
#' 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 &-1 &-1 & 0 & 0 & X\cr
#' 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 &-1 &-1 & Y\cr
#' 0 & 0 & 0 & 0 & 0 & 1 & 1 & 1 & 1 & 1 &-1 &-1 &-1 &-1 & AA\cr
#' \hline
#' T & A & B & X & Y & AA& AB& BA& BB& BC& XX& XY& YX& YY&
#' \end{array}\right],}
#' from which the 'structural sum matrix'
#' \mjsdeqn{\mathbf{S} = \left[\mathbf{C}' \; \mathbf{I}_8 \right]'}
#' may be easily obtained.
#' @family utilities
#'
#' @return A list with
#' \item{\code{Ut_reshape}}{matrix with rows and columns rearranged so that
#' each row has 1 as the first non-null element starting from the left,
#' and that the first r columns have 1 as the first element starting from
#' the bottom, \mjseqn{\mathbf{U}_{reshape}' = \mathbf{U}'[,\mbox{cid}]}.}
#' \item{\code{Ut_rref}}{reduced row echelon form of \mjseqn{\mathbf{U}'_{reshape}} with
#' \link{srref}.}
#' \item{\code{C}}{(\mjseqn{n_a \times n_b}) cross-sectional (contemporaneous) matrix
#' mapping the bottom level series into the higher level ones,
#' \mjseqn{\mathbf{U}_{srref}' = [\mathbf{I} \ -\mathbf{C}]}.}
#' \item{\code{cid}}{(\mjseqn{c \times 1}) vector of the column permutations of the matrix \mjseqn{\mathbf{U}'}.}
#' \item{\code{nb}}{number of bottom time series, \mjseqn{n_b}.}
#' \item{\code{na}}{number of upper time series, \mjseqn{n_a = r}.}
#' \item{\code{n}}{number of time series, \mjseqn{n_a + n_b = n = c}.}
#'
#' @references
#' Di Fonzo, T., Girolimetto, D. (2020), Cross-Temporal Forecast Reconciliation:
#' Optimal Combination Method and Heuristic Alternatives, Department of Statistical
#' Sciences, University of Padua, \href{https://arxiv.org/abs/2006.08570}{arXiv:2006.08570}.
#'
#' @export
#'
#' @examples
#' C1 <- Matrix(c(1,1,0,0,0,
#'                0,0,1,1,1), 2, byrow = TRUE, sparse = TRUE)
#' C2 <- Matrix(c(1,1,0,0,
#'                0,0,1,1), 2, byrow = TRUE, sparse = TRUE)
#' Ut <- rbind(c(1, rep(0, NROW(C1)), rep(-1, NCOL(C1)), rep(0, NROW(C2)), rep(0, NCOL(C2))),
#'             c(1, rep(0, NROW(C1)), rep(0, NCOL(C1)), rep(0, NROW(C2)), rep(-1, NCOL(C2))),
#'             cbind(0, Diagonal(NROW(C1)), -C1, Matrix(0, NROW(C1), NROW(C2)+NCOL(C2))),
#'             cbind(0, Matrix(0, NROW(C2), NROW(C1)+NCOL(C1)), Diagonal(NROW(C2)), -C2))
#' colnames(Ut) <- c("T", "A", "B", "AA", "AB", "BA", "BB", "BC",
#'                   "X", "Y", "XX", "XY", "YX", "YY")
#' rownames(Ut) <- c("T", "T", "A", "B", "X", "Y")
#' obj_ut2c <- ut2c(Ut, sparse = FALSE)
#' Ut_struc <- obj_ut2c$Ut_rref
#' S <- rbind(obj_ut2c$C, diag(1, obj_ut2c$nb))
#'
ut2c <- function(Ut, sparse = TRUE, verbose = FALSE){
  if(!is(Ut, "dgCMatrix")){
    Ut <- as(Ut, "dgCMatrix")
  }

  obj <- reshape_mat(mat = Ut)
  Ut_reshape <- obj$mat
  cid <- obj$cid

  out <- list()
  out$Ut_reshape <- Ut_reshape

  if(!isDiagonal(Ut_reshape[,c(1:NROW(Ut_reshape))]) | !all(Ut_reshape@x[1:NROW(Ut_reshape)]==1)){
    Ut_rref <- srref(A = Ut_reshape, sparse = TRUE, verbose = verbose)

    rs <- rowSums(abs(Ut_rref))==0
    if(any(rs)){
      message("rref(Ut) has ", sum(rs), " zeros rows. These rows have been removed")
      Ut_rref <- Ut_rref[!rs,,drop=FALSE]
    }

    colnames(Ut_rref) <- colnames(Ut_reshape)
    rownames(Ut_rref) <- colnames(Ut_reshape)[1:NROW(Ut_rref)]
    out$Ut_rref <- Ut_rref

    out$C <- -Ut_rref[,-c(1:NROW(Ut_rref)),drop=FALSE]
    colnames(out$C) <- colnames(Ut_reshape)[-c(1:NROW(Ut_rref))]
    rownames(out$C) <- colnames(Ut_reshape)[1:NROW(Ut_rref)]
    if(!sparse){
      out$Ut_reshape <- as.matrix(out$Ut_reshape)
      out$Ut_rref <- as.matrix(out$Ut_rref)
      out$C <- as.matrix(out$C)
    }
  }else{
    out$C <- -Ut_reshape[,-c(1:NROW(Ut_reshape)),drop=FALSE]
    colnames(out$C) <- colnames(Ut_reshape)[-c(1:NROW(Ut_reshape))]
    rownames(out$C) <- colnames(Ut_reshape)[1:NROW(Ut_reshape)]

    if(!sparse){
      out$Ut_reshape <- as.matrix(out$Ut_reshape)
      out$C <- as.matrix(out$C)
    }
  }

  names(cid) <- colnames(Ut)[cid]
  out$cid <- cid
  out$nb <- NCOL(out$C)
  out$na <- NROW(out$C)
  out$n <- NCOL(out$C)+NROW(out$C)
  return(out)
}

reshape_mat <- function(mat){
  if(!is(mat, "dgCMatrix")){
    mat <- as(mat, "dgCMatrix")
  }
  r <- mat@i+1
  c <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  c_new <- c[which(mat@x==1)][!duplicated(r[which(mat@x==1)])]
  id <- c(unique(c_new), c(1:NCOL(mat))[!c(1:NCOL(mat)) %in% unique(c_new)])
  #r_new <- r[which(mat@x==1)][!duplicated(r[which(mat@x==1)])] # not now
  out <- list()
  out$mat <- mat[, id]
  out$cid <- id
  return(out)
}


