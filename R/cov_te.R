#' Temporal covariance matrix approximation
#'
#' @description
#' This function provides an approximation of the temporal base forecasts errors
#' covariance matrix using different reconciliation methods (see Di Fonzo and Girolimetto, 2023).
#'
#' @usage
#' tecov(comb, agg_order = NULL, res = NULL, tew = "sum",
#'       mse = TRUE, shrink_fun = shrink_estim, ...)
#'
#' @inheritParams terec
#' @param comb A string specifying the reconciliation method.
#'   \itemize{
#'      \item Ordinary least squares:
#'      \itemize{
#'      \item "\code{ols}" (\emph{default}) - identity error covariance.
#'      }
#'     \item Weighted least squares:
#'      \itemize{
#'      \item "\code{str}" - structural variances.
#'      \item "\code{wlsh}" - hierarchy variances (uses \code{res}).
#'      \item "\code{wlsv}" - series variances (uses \code{res}).
#'      }
#'     \item Generalized least squares (uses \code{res}):
#'      \itemize{
#'      \item "\code{acov}" - series auto-covariance.
#'      \item "\code{strar1}" - structural Markov covariance.
#'      \item "\code{sar1}" - series Markov covariance.
#'      \item "\code{har1}" - hierarchy Markov covariance.
#'      \item "\code{shr}"/"\code{sam}" - shrunk/sample covariance.
#'      }
#'   }
#' @param mse If \code{TRUE} (\emph{default}) the residuals used to compute the covariance
#' matrix are not mean-corrected.
#' @param shrink_fun Shrinkage function of the covariance matrix, [shrink_estim] (\emph{default})
#' @param ... Not used.
#'
#' @returns A (\eqn{(k^\ast+m) \times (k^\ast+m)}) symmetric matrix.
#'
#' @references
#' Di Fonzo, T. and Girolimetto, D. (2023), Cross-temporal forecast reconciliation:
#' Optimal combination method and heuristic alternatives, \emph{International Journal
#' of Forecasting}, 39, 1, 39-57. \doi{10.1016/j.ijforecast.2021.08.004}
#'
#' @examples
#' # (7 x 70) in-sample residuals matrix (simulated), agg_order = 4
#' res <- rnorm(70)
#'
#' cov1 <- tecov("ols", agg_order = 4)                 # OLS methods
#' cov2 <- tecov("str", agg_order = 4)                 # STRC methods
#' cov3 <- tecov("wlsv", agg_order = 4, res = res)     # WLSv methods
#' cov4 <- tecov("wlsh", agg_order = 4, res = res)     # WLSh methods
#' cov5 <- tecov("acov", agg_order = 4, res = res)     # ACOV methods
#' cov6 <- tecov("strar1", agg_order = 4, res = res)   # STRAR1 methods
#' cov7 <- tecov("har1", agg_order = 4, res = res)     # HAR1 methods
#' cov8 <- tecov("sar1", agg_order = 4, res = res)     # SAR1 methods
#' cov9 <- tecov("shr", agg_order = 4, res = res)      # SHR methods
#' cov10 <- tecov("sam", agg_order = 4, res = res)     # SAM methods
#'
#' # Custom covariance matrix
#' tecov.ols2 <- function(comb, x) diag(x)
#' tecov(comb = "ols2", x = 7) # == tecov("ols", agg_order = 4)
#'
#' @aliases tecov.default tecov.ols tecov.str tecov.wlsv tecov.wlsh tecov.acov tecov.strar1
#' @aliases tecov.sar1 tecov.har1 tecov.shr tecov.sam
#'
#' @family Framework: temporal
#' @export
tecov <- function(comb = "ols", agg_order = NULL, res = NULL, tew = "sum",
                  mse = TRUE, shrink_fun = shrink_estim, ...){
  if(is(comb, "Matrix") | is(comb, "matrix")){
    comb
  }else if(is.character(comb)){
    class(comb) <- comb
    UseMethod("tecov", comb)
  }else{
    cli_abort("{.arg comb} is not a character.", call = NULL)
  }
}

#' @export
tecov.default <- function(comb, agg_order = NULL, res = NULL, tew = "sum",
                          mse = TRUE, shrink_fun = shrink_estim, ...){
  cli_abort("Please, provide a valid covariance matrix approach ({.arg comb}).", call = NULL)
}

#' @export
tecov.ols <- function(comb = "ols", ..., agg_order = NULL){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg agg_order} is NULL.", call = NULL)
  }
  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }
  .sparseDiagonal(sum(max(kset)/kset))
}

#' @export
tecov.str <- function(comb = "str", ..., agg_order = NULL, tew = "sum", strc_mat = NULL){
  if(is.null(strc_mat)){
    if(is.null(agg_order)){
      cli_abort("Argument {.arg agg_order} is NULL.", call = NULL)
    }
    strc_mat <- tetools(agg_order = agg_order, tew = tew)$strc_mat
  }

  .sparseDiagonal(x = rowSums(strc_mat))
}

#' @export
tecov.wlsv <- function(comb = "wlsv", ..., agg_order = NULL, res = NULL, mse = TRUE){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg agg_order} is NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) != 1){
    cli_abort("{.arg res} is not a vector.", call = NULL)
  }else if(length(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} length.", call = NULL)
  }

  m <- max(kset)
  N <- length(res) / sum(max(kset)/kset)
  var_freq <- sapply(kset, function(x)
    sample_estim(res[rep(kset, (m/kset) * N) == x], mse = mse))

  .sparseDiagonal(x = rep(var_freq, (m/kset)))
}

#' @export
tecov.wlsh <- function(comb = "wlsh", ..., agg_order = NULL, res = NULL, mse = TRUE){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg agg_order} is NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) != 1){
    cli_abort("{.arg res} is not a vector.", call = NULL)
  }else if(length(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} length.", call = NULL)
  }

  N <- length(res) / sum(max(kset)/kset)
  res_mat <- vec2hmat(vec = res, h = N, kset = kset)
  .sparseDiagonal(x = apply(res_mat, 2, function(x) ifelse(mse, sum(x^2)/length(x), var(x))))
}

#' @export
tecov.acov <- function(comb = "acov", ..., agg_order = NULL, res = NULL, mse = TRUE){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg agg_order} is NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) != 1){
    cli_abort("{.arg res} is not a vector.", call = NULL)
  }else if(length(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} length.", call = NULL)
  }

  N <- length(res) / sum(max(kset)/kset)
  m <- max(kset)
  res_mat <- vec2hmat(vec = res, h = N, kset = kset)

  bdiag(lapply(kset, function(x)
    sample_estim(res_mat[, rep(kset, (m/kset)) == x], mse = mse)))
}

#' @export
tecov.strar1 <- function(comb = "strar1", ..., agg_order = NULL, res = NULL, tew = "sum"){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg agg_order} is NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  tmp <- tetools(agg_order = agg_order, tew = tew)
  kset <- tmp$set
  strc_mat <- tmp$strc_mat

  if(NCOL(res) != 1){
    cli_abort("{.arg res} is not a vector.", call = NULL)
  }else if(length(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} length.", call = NULL)
  }

  N <- length(res) / sum(max(kset)/kset)
  m <- max(kset)
  rho <- lapply(kset, function(x)
    stats::acf(stats::na.omit(res[rep(kset, (m/kset) * N) == x]), 1, plot = F)$acf[2, 1, 1])
  expo <- lapply((m/kset), function(x) toeplitz(1:x) - 1)

  Gam <- Matrix::bdiag(Map(function(x, y) x^y, x = rho, y = expo))
  Ostr2 <- .sparseDiagonal(x = apply(strc_mat, 1, sum))^0.5

  return(Ostr2 %*% Gam %*% Ostr2)
}

#' @export
tecov.sar1 <- function(comb = "sar1", ..., agg_order = NULL, res = NULL,
                       tew = "sum", mse = TRUE){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg agg_order} is NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) != 1){
    cli_abort("{.arg res} is not a vector.", call = NULL)
  }else if(length(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} length.", call = NULL)
  }

  N <- length(res) / sum(max(kset)/kset)
  m <- max(kset)
  rho <- lapply(kset, function(x)
    stats::acf(stats::na.omit(res[rep(kset, (m/kset) * N) == x]), 1, plot = F)$acf[2, 1, 1])
  expo <- lapply((m/kset), function(x) toeplitz(1:x) - 1)

  Gam <- Matrix::bdiag(Map(function(x, y) x^y, x = rho, y = expo))
  var_freq <- sapply(kset, function(x) sample_estim(res[rep(kset, (m/kset) * N) == x],
                                                    mse = mse))
  Os2 <- .sparseDiagonal(x = rep(var_freq, (m/kset)))^0.5
  return(Os2 %*% Gam %*% Os2)
}

#' @export
tecov.har1 <- function(comb = "har1", ..., agg_order = NULL, res = NULL, mse = TRUE){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg agg_order} is NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) != 1){
    cli_abort("{.arg res} is not a vector.", call = NULL)
  }else if(length(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} length.", call = NULL)
  }

  N <- length(res) / sum(max(kset)/kset)
  m <- max(kset)
  res_mat <- vec2hmat(vec = res, h = N, kset = kset)

  rho <- lapply(kset, function(x)
    stats::acf(stats::na.omit(res[rep(kset, (m/kset) * N) == x]), 1, plot = F)$acf[2, 1, 1])
  expo <- lapply((m/kset), function(x) toeplitz(1:x) - 1)

  Gam <- Matrix::bdiag(Map(function(x, y) x^y, x = rho, y = expo))
  diagO <- diag(sample_estim(res_mat, mse = mse))
  Oh2 <- .sparseDiagonal(x = diagO)^0.5
  return(Oh2 %*% Gam %*% Oh2)
}

#' @export
tecov.shr <- function(comb = "shr", ..., agg_order = NULL, res = NULL,
                      mse = TRUE, shrink_fun = shrink_estim){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg agg_order} is NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) != 1){
    cli_abort("{.arg res} is not a vector.", call = NULL)
  }else if(length(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} length.", call = NULL)
  }

  N <- length(res) / sum(max(kset)/kset)
  res_mat <- vec2hmat(vec = res, h = N, kset = kset)
  shrink_fun(res_mat, mse = mse)
}

#' @export
tecov.sam <- function(comb = "sam", ..., agg_order = NULL, res = NULL, mse = TRUE){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg agg_order} is NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) != 1){
    cli_abort("{.arg res} is not a vector.", call = NULL)
  }else if(length(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} length.", call = NULL)
  }

  N <- length(res) / sum(max(kset)/kset)
  res_mat <- vec2hmat(vec = res, h = N, kset = kset)
  sample_estim(res_mat, mse = mse)
}
