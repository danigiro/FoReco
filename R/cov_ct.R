#' Cross-temporal covariance matrix approximation
#'
#' @description
#' This function provides an approximation of the cross-temporal base forecasts errors
#' covariance matrix using different reconciliation methods (Di Fonzo and Girolimetto, 2023,
#' Girolimetto et al., 2023).
#'
#' @usage
#' ctcov(comb = "ols", n = NULL, agg_mat = NULL, agg_order = NULL, res = NULL,
#'       tew = "sum", mse = TRUE, shrink_fun = shrink_estim, ...)
#'
#' @inheritParams ctrec
#' @param comb A string specifying the reconciliation method.
#'   \itemize{
#'      \item Ordinary least squares:
#'      \itemize{
#'      \item "\code{ols}" (\emph{default}) - identity error covariance.
#'      }
#'     \item Weighted least squares:
#'      \itemize{
#'      \item "\code{str}" - structural variances.
#'      \item "\code{csstr}" - cross-sectional structural variances.
#'      \item "\code{testr}" - temporal structural variances.
#'      \item "\code{wlsh}" - hierarchy variances (uses \code{res}).
#'      \item "\code{wlsv}" - series variances (uses \code{res}).
#'      }
#'     \item Generalized least squares (uses \code{res}):
#'      \itemize{
#'      \item "\code{acov}" - series auto-covariance.
#'      \item "\code{bdshr}"/"\code{bdsam}" - shrunk/sample block diagonal cross-sectional covariance.
#'      \item "\code{Sshr}"/"\code{Ssam}" - series shrunk/sample covariance.
#'      \item "\code{shr}"/"\code{sam}" - shrunk/sample covariance.
#'      \item "\code{hbshr}"/"\code{hbsam}" - shrunk/sample high frequency bottom time series covariance.
#'      \item "\code{bshr}"/"\code{bsam}" - shrunk/sample bottom time series covariance.
#'      \item "\code{hshr}"/"\code{hsam}" - shrunk/sample high frequency covariance.
#'      }
#'   }
#' @param n Cross-sectional number of variables.
#' @param mse If \code{TRUE} (\emph{default}) the residuals used to compute the covariance
#' matrix are not mean-corrected.
#' @param shrink_fun Shrinkage function of the covariance matrix, [shrink_estim] (\emph{default}).
#' @param ... Not used.
#'
#' @returns A (\eqn{n(k^\ast+m) \times n(k^\ast+m)}) symmetric matrix.
#'
#' @references
#' Di Fonzo, T. and Girolimetto, D. (2023), Cross-temporal forecast reconciliation:
#' Optimal combination method and heuristic alternatives, \emph{International Journal
#' of Forecasting}, 39, 1, 39-57. \doi{10.1016/j.ijforecast.2021.08.004}
#'
#' Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J. (2024),
#' Cross-temporal probabilistic forecast reconciliation: Methodological and
#' practical issues. \emph{International Journal of Forecasting}, 40, 3, 1134-1151.
#' \doi{10.1016/j.ijforecast.2023.10.003}
#'
#' @examples
#' set.seed(123)
#' # Aggregation matrix for Z = X + Y
#' A <- t(c(1,1))
#' # (3 x 70) in-sample residuals matrix (simulated),
#' # agg_order = 4 (annual-quarterly)
#' res <- rbind(rnorm(70), rnorm(70), rnorm(70))
#'
#' cov1 <- ctcov("ols", n = 3, agg_order = 4)                     # OLS methods
#' cov2 <- ctcov("str", agg_mat = A, agg_order = 4)               # STR methods
#' cov3 <- ctcov("csstr", agg_mat = A, agg_order = 4)             # CSSTR methods
#' cov4 <- ctcov("testr", n = 3, agg_order = 4)                   # TESTR methods
#' cov5 <- ctcov("wlsv", agg_order = 4, res = res)                # WLSv methods
#' cov6 <- ctcov("wlsh", agg_order = 4, res = res)                # WLSh methods
#' cov7 <- ctcov("shr", agg_order = 4, res = res)                 # SHR methods
#' cov8 <- ctcov("sam", agg_order = 4, res = res)                 # SAM methods
#' cov9 <- ctcov("acov", agg_order = 4, res = res)                # ACOV methods
#' cov10 <- ctcov("Sshr", agg_order = 4, res = res)               # Sshr methods
#' cov11 <- ctcov("Ssam", agg_order = 4, res = res)               # Ssam methods
#' cov12 <- ctcov("hshr", agg_order = 4, res = res)               # Hshr methods
#' cov13 <- ctcov("hsam", agg_order = 4, res = res)               # Hsam methods
#' cov14 <- ctcov("hbshr", agg_mat = A, agg_order = 4, res = res) # HBshr methods
#' cov15 <- ctcov("hbsam", agg_mat = A, agg_order = 4, res = res) # HBsam methods
#' cov16 <- ctcov("bshr", agg_mat = A, agg_order = 4, res = res)  # Bshr methods
#' cov17 <- ctcov("bsam", agg_mat = A, agg_order = 4, res = res)  # Bsam methods
#' cov18 <- ctcov("bdshr", agg_order = 4, res = res)              # BDshr methods
#' cov19 <- ctcov("bdsam", agg_order = 4, res = res)              # BDsam methods
#'
#' # Custom covariance matrix
#' ctcov.ols2 <- function(comb, x) diag(x)
#' cov20 <- ctcov(comb = "ols2", x = 21) # == ctcov("ols", n = 3, agg_order = 4)
#'
#' @aliases ctcov.default ctcov.ols ctcov.str ctcov.csstr ctcov.testr ctcov.sam ctcov.wlsh
#' @aliases ctcov.wlsv ctcov.shr ctcov.acov ctcov.Ssam ctcov.Sshr ctcov.bdshr ctcov.bdsam
#' @aliases ctcov.hbshr ctcov.hbsam ctcov.hshr ctcov.hsam ctcov.bshr ctcov.bsam
#'
#' @family Framework: cross-temporal
#' @export
ctcov <- function(comb = "ols", n = NULL, agg_mat = NULL, agg_order = NULL, res = NULL,
                  tew = "sum", mse = TRUE, shrink_fun = shrink_estim, ...){
  if(is(comb, "Matrix") | is(comb, "matrix")){
    comb
  }else if(is.character(comb)){
    class(comb) <- comb
    UseMethod("ctcov", comb)
  }else{
    cli_abort("{.arg comb} is not a character.", call = NULL)
  }
}

#' @export
ctcov.default <- function(comb, n = NULL, agg_mat = NULL, agg_order = NULL, res = NULL,
                          tew = "sum", mse = TRUE, shrink_fun = shrink_estim, ...){
  cli_abort("Please, provide a valid covariance matrix approach ({.arg comb}).", call = NULL)
}

#' @export
ctcov.ols <- function(comb = "ols", ..., n = NULL, agg_order = NULL){
  if(is.null(n) || is.null(agg_order)){
    cli_abort("Argument {.arg n} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  .sparseDiagonal(n * sum(max(kset)/kset)) # Omega
}

#' @export
ctcov.str <- function(comb = "str", ..., agg_mat = NULL, agg_order = NULL, tew = "sum", strc_mat = NULL){
  if(is.null(strc_mat)){
    if(is.null(agg_mat) || is.null(agg_order)){
      cli_abort("Argument {.arg agg_mat} and/or {.arg agg_order} are NULL.", call = NULL)
    }

    strc_mat <- cttools(agg_mat = agg_mat, agg_order = agg_order, tew = tew)$strc_mat
  }

  .sparseDiagonal(x = rowSums(strc_mat)) # Omega
}

#' @export
ctcov.csstr <- function(comb = "csstr", ..., agg_mat = NULL, agg_order = NULL){
  if(is.null(agg_mat) || is.null(agg_order)){
    cli_abort("Argument {.arg agg_mat} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  Scs <- cstools(agg_mat = agg_mat)$strc_mat
  .sparseDiagonal(x = rep(rowSums(Scs), each = sum(max(kset)/kset))) # Omega
}

#' @export
ctcov.testr <- function(comb = "testr", ..., n = NULL, agg_order = NULL, tew = "sum"){
  if(is.null(n) || is.null(agg_order)){
    cli_abort("Argument {.arg n} and/or {.arg agg_order} are NULL.", call = NULL)
  }
  Ste <- tetools(agg_order = agg_order, tew = tew)$strc_mat
  .sparseDiagonal(x = rep(rowSums(Ste), n))
}

#' @export
ctcov.sam <- function(comb = "sam", ..., agg_order = NULL, res = NULL, mse = TRUE){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg n} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  n <- NROW(res)

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
  }

  res_mat <- mat2hmat(res, h = NCOL(res) / sum(max(kset)/kset), kset = kset, n = n)
  sample_estim(res_mat, mse = mse)
}

#' @export
ctcov.wlsh <- function(comb = "wlsh", ..., agg_order = NULL, res = NULL, mse = TRUE){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg n} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  n <- NROW(res)

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
  }

  res_mat <- mat2hmat(res, h = NCOL(res) / sum(max(kset)/kset), kset = kset, n = n)
  res_mat <- na.omit(res_mat)
  .sparseDiagonal(x = apply(res_mat, 2, function(x) ifelse(mse, sum(x^2)/length(x), var(x))))
}

#' @export
ctcov.wlsv <- function(comb = "wlsv", ..., agg_order = NULL, res = NULL, mse = TRUE){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg n} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  n <- NROW(res)

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
  }

  N <- NCOL(res) / sum(max(kset)/kset)
  m <- max(kset)
  var_freq <- apply(res, 1, function(z) sapply(kset, function(x)
    sample_estim(z[rep(kset, (m/kset) * N) == x], mse = mse)))
  .sparseDiagonal(x = rep(as.vector(var_freq), rep((m/kset), n)))
}

#' @export
ctcov.shr <- function(comb = "shr", ..., agg_order = NULL, res = NULL,
                      mse = TRUE, shrink_fun = shrink_estim){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg n} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  n <- NROW(res)

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
  }

  res_mat <- mat2hmat(res, h = NCOL(res) / sum(max(kset)/kset), kset = kset, n = n)
  shrink_fun(res_mat, mse = mse)
}

#' @export
ctcov.acov <- function(comb = "acov", ..., agg_order = NULL, res = NULL, mse = TRUE){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg n} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  n <- NROW(res)

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
  }

  m <- max(kset)
  res_mat <- mat2hmat(res, h = NCOL(res) / sum(max(kset)/kset), kset = kset, n = n)
  mat1 <- bdiag(rep(lapply((m/kset), function(x) matrix(1, nrow = x, ncol = x)), n))
  sample_estim(res_mat, mse = mse) * mat1
}

#' @export
ctcov.Ssam <- function(comb = "Ssam", ..., agg_order = NULL, res = NULL, mse = TRUE){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg n} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  n <- NROW(res)

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
  }

  m <- max(kset)
  res_mat <- mat2hmat(res, h = NCOL(res) / sum(max(kset)/kset), kset = kset, n = n)
  mat1 <- bdiag(replicate(n, matrix(1, nrow = sum(max(kset)/kset),
                                    ncol = sum(max(kset)/kset)), simplify = FALSE))
  sample_estim(res_mat, mse = mse) * mat1
}

#' @export
ctcov.Sshr <- function(comb = "Sshr", ..., agg_order = NULL, res = NULL,
                       mse = TRUE, shrink_fun = shrink_estim){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg n} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  n <- NROW(res)

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
  }

  res_mat <- mat2hmat(res, h = NCOL(res) / sum(max(kset)/kset), kset = kset, n = n)
  shrink <- lapply(1:n, function(x)
    shrink_fun(res_mat[, rep(1:n, each = sum(max(kset)/kset)) == x], mse = mse))
  bdiag(shrink)
}

#' @export
ctcov.bdshr <- function(comb = "bdshr", ..., agg_order = NULL, res = NULL,
                        mse = TRUE, shrink_fun = shrink_estim){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg n} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  n <- NROW(res)

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
  }

  m <- max(kset)
  N <- NCOL(res) / sum(max(kset)/kset)
  Wlist <- lapply(kset, function(x) shrink_fun(t(res[, rep(kset, N * (m/kset)) == x]), mse = mse))
  W <- rep(Wlist, (m/kset))
  P <- commat(n, sum(max(kset)/kset))
  P %*% bdiag(W) %*% t(P)
}

#' @export
ctcov.bdsam <- function(comb = "bdsam", ..., agg_order = NULL, res = NULL, mse = TRUE){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg n} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  n <- NROW(res)

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
  }

  m <- max(kset)
  N <- NCOL(res) / sum(max(kset)/kset)
  Wlist <- lapply(kset, function(x) sample_estim(t(res[, rep(kset, N * (m/kset)) == x]), mse = mse))
  W <- rep(Wlist, (m/kset))
  P <- commat(n, sum(max(kset)/kset))
  P %*% bdiag(W) %*% t(P)
}

#' @export
ctcov.hbshr <- function(comb = "hbshr", ..., agg_mat = NULL, agg_order = NULL, res = NULL,
                        tew = "sum", mse = TRUE, shrink_fun = shrink_estim, eps_ridge = NULL){
  if(is.null(agg_order) || is.null(agg_mat)){
    cli_abort("Argument {.arg agg_mat} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  strc_mat <- cttools(agg_mat = agg_mat, agg_order = agg_order, tew = tew)$strc_mat

  if(NCOL(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
  }

  n <- sum(dim(agg_mat))
  na <- NROW(agg_mat)
  m <- max(kset)
  kid <- rep(rep(kset, m/kset), n)
  nid <- rep(1:n, each = sum(max(kset)/kset))
  res_mat <- mat2hmat(res, h = NCOL(res) / sum(max(kset)/kset), kset = kset, n = n)
  cov <- shrink_fun(res_mat[, nid > na & kid == 1], mse = mse)
  Omega <- strc_mat %*% cov %*% t(strc_mat)

  if(is.null(eps_ridge)){
    eigenvalues <- eigen(Omega)$values
    eps_ridge <- min(eigenvalues[eigenvalues > 1e-6])
  }
  print(eps_ridge)
  Omega + diag(eps_ridge, dim(Omega)[1])
}

#' @export
ctcov.hbsam <- function(comb = "hbsam", ..., agg_mat = NULL, agg_order = NULL, res = NULL,
                        tew = "sum", mse = TRUE, eps_ridge = NULL){
  if(is.null(agg_order) || is.null(agg_mat)){
    cli_abort("Argument {.arg agg_mat} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  strc_mat <- cttools(agg_mat = agg_mat, agg_order = agg_order, tew = tew)$strc_mat

  if(NCOL(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
  }

  n <- sum(dim(agg_mat))
  na <- NROW(agg_mat)
  m <- max(kset)
  N <- NCOL(res) / sum(max(kset)/kset)
  kid <- rep(rep(kset, m/kset), n)
  nid <- rep(1:n, each = sum(max(kset)/kset))
  res_mat <- mat2hmat(res, h = N, kset = kset, n = n)
  cov <- sample_estim(res_mat[, nid > na & kid == 1], mse = mse)
  Omega <- strc_mat %*% cov %*% t(strc_mat)

  if(is.null(eps_ridge)){
    eigenvalues <- eigen(Omega)$values
    eps_ridge <- min(eigenvalues[eigenvalues > 1e-6])
  }
  Omega + diag(eps_ridge, dim(Omega)[1])
}

#' @export
ctcov.hshr <- function(comb = "hshr", ..., agg_order = NULL, res = NULL, tew = "sum",
                       mse = TRUE, shrink_fun = shrink_estim, eps_ridge = NULL){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg n} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  n <- NROW(res)
  tmp <- tetools(agg_order = agg_order, tew = tew)
  kset <- tmp$set

  if(NCOL(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
  }

  Ste <- tmp$strc_mat
  m <- max(kset)
  N <- NCOL(res) / sum(max(kset)/kset)
  kid <- rep(rep(kset, m/kset), n)
  res_mat <- mat2hmat(res, h = N, kset = kset, n = n)
  cov <- shrink_fun(res_mat[, kid == 1], mse = mse)
  kSte <- kronecker(Diagonal(n), Ste)
  Omega <- kSte %*% cov %*% t(kSte)

  if(is.null(eps_ridge)){
    eigenvalues <- eigen(Omega)$values
    eps_ridge <- min(eigenvalues[eigenvalues > 1e-6])
  }
  Omega + diag(eps_ridge, dim(Omega)[1])
}

#' @export
ctcov.hsam <- function(comb = "hsam", ..., agg_order = NULL, res = NULL,
                       tew = "sum", mse = TRUE, eps_ridge = NULL){
  if(is.null(agg_order)){
    cli_abort("Argument {.arg n} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  n <- NROW(res)
  tmp <- tetools(agg_order = agg_order, tew = tew)
  kset <- tmp$set

  if(NCOL(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
  }

  Ste <- tmp$strc_mat
  m <- max(kset)
  N <- NCOL(res) / sum(max(kset)/kset)
  kid <- rep(rep(kset, m/kset), n)
  res_mat <- mat2hmat(res, h = N, kset = kset, n = n)
  cov <- sample_estim(res_mat[, kid == 1], mse = mse)
  kSte <- kronecker(Diagonal(n), Ste)
  Omega <- kSte %*% cov %*% t(kSte)

  if(is.null(eps_ridge)){
    eigenvalues <- eigen(Omega)$values
    eps_ridge <- min(eigenvalues[eigenvalues > 1e-6])
  }
  Omega + diag(eps_ridge, dim(Omega)[1])
}

#' @export
ctcov.bshr <- function(comb = "bshr", ..., agg_mat = NULL, agg_order = NULL, res = NULL,
                       mse = TRUE, shrink_fun = shrink_estim, eps_ridge = NULL){
  if(is.null(agg_order) || is.null(agg_mat)){
    cli_abort("Argument {.arg agg_mat} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
  }

  Scs <- cstools(agg_mat = agg_mat)$strc_mat
  n <- sum(dim(agg_mat))
  na <- NROW(agg_mat)
  N <- NCOL(res) / sum(max(kset)/kset)
  nid <- rep(1:n, each = sum(max(kset)/kset))
  res_mat <- mat2hmat(res, h = N, kset = kset, n = n)
  cov <- shrink_fun(res_mat[, nid > na], mse = mse)
  kScs <- kronecker(Scs, Diagonal(sum(max(kset)/kset)))
  Omega <- kScs %*% cov %*% t(kScs)
  if(is.null(eps_ridge)){
    eigenvalues <- eigen(Omega)$values
    eps_ridge <- min(eigenvalues[eigenvalues > 1e-6])
  }
  Omega + diag(eps_ridge, dim(Omega)[1])
}

#' @export
ctcov.bsam <- function(comb = "bsam", ..., agg_mat = NULL, agg_order = NULL, res = NULL,
                       mse = TRUE, eps_ridge = NULL){
  if(is.null(agg_order) || is.null(agg_mat)){
    cli_abort("Argument {.arg agg_mat} and/or {.arg agg_order} are NULL.", call = NULL)
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }

  if(NCOL(res) %% sum(max(kset)/kset) != 0){
    cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
  }

  Scs <- cstools(agg_mat = agg_mat)$strc_mat
  n <- sum(dim(agg_mat))
  na <- NROW(agg_mat)
  N <- NCOL(res) / sum(max(kset)/kset)
  nid <- rep(1:n, each = sum(max(kset)/kset))
  res_mat <- mat2hmat(res, h = N, kset = kset, n = n)
  cov <- sample_estim(res_mat[, nid > na], mse = mse)
  kScs <- kronecker(Scs, Diagonal(sum(max(kset)/kset)))
  Omega <- kScs %*% cov %*% t(kScs)
  if(is.null(eps_ridge)){
    eigenvalues <- eigen(Omega)$values
    eps_ridge <- min(eigenvalues[eigenvalues > 1e-6])
  }
  Omega + diag(eps_ridge, dim(Omega)[1])
}
