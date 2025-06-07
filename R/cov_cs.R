#' Cross-sectional covariance matrix approximation
#'
#' @description
#' This function provides an approximation of the cross-sectional base forecasts errors
#' covariance matrix using different reconciliation methods (see Wickramasuriya et al.,
#' 2019 and Di Fonzo and Girolimetto, 2023).
#'
#' @usage
#' cscov(comb = "ols", n = NULL, agg_mat = NULL, res, mse = TRUE,
#'       shrink_fun = shrink_estim, ...)
#'
#' @param n Number of variables (\eqn{n = n_a + n_b}).
#' @inheritParams csrec
#' @param comb A string specifying the reconciliation method.
#'   \itemize{
#'      \item Ordinary least squares:
#'      \itemize{
#'      \item "\code{ols}" (\emph{default}) - identity error covariance matrix.
#'      }
#'     \item Weighted least squares:
#'      \itemize{
#'      \item "\code{str}" - structural variances.
#'      \item "\code{wls}" - series variances (uses \code{res}).
#'      }
#'     \item Generalized least squares (uses \code{res}):
#'      \itemize{
#'      \item "\code{shr}" - shrunk covariance (Wickramasuriya et al., 2019).
#'      \item "\code{oasd}" - oracle shrunk covariance (Ando and Xiao, 2023).
#'      \item "\code{sam}" - sample covariance.
#'      }
#'   }
#' @param mse If \code{TRUE} (\emph{default}) the residuals used to compute the covariance
#' matrix are not mean-corrected.
#' @param shrink_fun Shrinkage function of the covariance matrix, [shrink_estim] (\emph{default}).
#' @param ... Not used.
#'
#' @returns A (\eqn{n \times n}) symmetric positive (semi-)definite matrix.
#' @aliases cscov.default cscov.ols cscov.str cscov.wls cscov.shr cscov.sam
#' @family Framework: cross-sectional
#'
#' @examples
#' # Aggregation matrix for Z = X + Y
#' A <- t(c(1,1))
#' # (10 x 3) in-sample residuals matrix (simulated)
#' res <- t(matrix(rnorm(n = 30), nrow = 3))
#'
#' cov1 <- cscov("ols", n = 3)          # OLS methods
#' cov2 <- cscov("str", agg_mat = A)    # STR methods
#' cov3 <- cscov("wls", res = res)      # WLS methods
#' cov4 <- cscov("shr", res = res)      # SHR methods
#' cov5 <- cscov("sam", res = res)      # SAM methods
#'
#' # Custom covariance matrix
#' cscov.ols2 <- function(comb, x) diag(x)
#' cscov(comb = "ols2", x = 3) # == cscov("ols", n = 3)
#'
#' @references
#' Ando, S., and Xiao, M. (2023), High-dimensional covariance matrix estimation:
#' shrinkage toward a diagonal target. \emph{IMF Working Papers}, 2023(257), A001.
#'
#' Di Fonzo, T. and Girolimetto, D. (2023), Cross-temporal forecast reconciliation:
#' Optimal combination method and heuristic alternatives, \emph{International Journal
#' of Forecasting}, 39, 1, 39-57. \doi{10.1016/j.ijforecast.2021.08.004}
#'
#' Wickramasuriya, S.L., Athanasopoulos, G. and Hyndman, R.J. (2019), Optimal forecast
#' reconciliation for hierarchical and grouped time series through trace minimization,
#' \emph{Journal of the American Statistical Association}, 114, 526, 804-819.
#' \doi{10.1080/01621459.2018.1448825}
#'
#' @export
cscov <- function(comb = "ols", n = NULL, agg_mat = NULL, res = NULL, mse = TRUE,
                  shrink_fun = shrink_estim, ...){
  if(is(comb, "Matrix") | is(comb, "matrix")){
    comb
  }else if(is.character(comb)){
    class(comb) <- comb
    UseMethod("cscov", comb)
  }else{
    cli_abort("{.arg comb} is not a character.", call = NULL)
  }
}

#.cscov <- function(comb, n, agg_mat, res, mse = TRUE, shrink_fun = shrink_estim, ...){
#  class(comb) <- comb
#  UseMethod(".cscov", comb)
#}

#' @export
cscov.default <- function(comb, n = NULL, agg_mat = NULL, res = NULL, mse = TRUE,
                          shrink_fun = shrink_estim, ...){
  cli_abort("Please, provide a valid covariance matrix approach ({.arg comb}).", call = NULL)
}

#' @export
cscov.ols <- function(comb = "ols", ..., n = NULL){
  if(is.null(n)){
    cli_abort("Argument {.arg n} is NULL.", call = NULL)
  }
  .sparseDiagonal(n)
}

#' @export
cscov.str <- function(comb = "str", ..., agg_mat = NULL, strc_mat = NULL){
  if(is.null(strc_mat)){
    if(is.null(agg_mat)){
      cli_abort("Argument {.arg agg_mat} is NULL.", call = NULL)
    }

    strc_mat <- cstools(agg_mat = agg_mat)$strc_mat
  }

  .sparseDiagonal(x = rowSums(abs(strc_mat)))
}

#' @export
cscov.wls <- function(comb = "wls", ..., res = NULL, mse = TRUE){
  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }
  .sparseDiagonal(x = apply(res, 2, function(x) ifelse(mse, sum(x^2, na.rm = TRUE)/sum(!is.na(x)),
                                                       var(x, na.rm = TRUE))))
}

#' @export
cscov.shr <- function(comb = "shr", ..., res = NULL, mse = TRUE,
                      shrink_fun = shrink_estim){
  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }
  shrink_fun(res, mse = mse)
}

#' @export
cscov.oasd <- function(comb = "shr", ..., res = NULL, mse = TRUE){
  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }
  shrink_oasd(res, mse = mse)
}

#' @export
cscov.sam <- function(comb = "sam", ..., res = NULL, mse = TRUE){
  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }
  sample_estim(res, mse = mse)
}
