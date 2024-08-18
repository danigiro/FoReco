#' Projection matrix for optimal combination cross-sectional reconciliation
#'
#' This function computes the projection or the mapping matrix
#' \eqn{\mathbf{M}} and \eqn{\mathbf{G}}, respectively, such that
#' \eqn{\widetilde{\mathbf{y}} = \mathbf{M}\widehat{\mathbf{y}} = \mathbf{S}_{cs}\mathbf{G}\widehat{\mathbf{y}}},
#' where \eqn{\widetilde{\mathbf{y}}} is the vector of the reconciled forecasts,
#' \eqn{\widehat{\mathbf{y}}} is the vector of the base forecasts,
#' \eqn{\mathbf{S}_{cs}} is the cross-sectional structural matrix, and \eqn{\mathbf{M} = \mathbf{S}_{cs}\mathbf{G}}.
#' For further information regarding on the structure of these matrices,
#' refer to Girolimetto et al. (2023).
#'
#' @inheritParams csrec
#' @param mat A string specifying which matrix to return:
#' "\code{M}" (\emph{default}) for \eqn{\mathbf{M}} and "\code{G}" for \eqn{\mathbf{G}}.
#' @inheritDotParams cscov mse shrink_fun
#'
#' @return The projection matrix \eqn{\mathbf{M}} (\code{mat = "M"}) or
#' the mapping matrix \eqn{\mathbf{G}} (\code{mat = "G"}).
#'
#' @references
#' Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J. (2024),
#' Cross-temporal probabilistic forecast reconciliation: Methodological and
#' practical issues. \emph{International Journal of Forecasting}, 40, 3, 1134-1151.
#' \doi{10.1016/j.ijforecast.2023.10.003}
#'
#' @examples
#' # Cross-sectional framework
#' A <- t(c(1,1)) # Aggregation matrix for Z = X + Y
#' Mcs <- csprojmat(agg_mat = A, comb = "ols")
#' Gcs <- csprojmat(agg_mat = A, comb = "ols", mat = "G")
#'
#' @family Utilities
#' @export
csprojmat <- function(agg_mat, cons_mat, comb = "ols", res = NULL, mat = "M", ...){
  mat <- match.arg(mat, c("M", "G"))

  if(missing(agg_mat) && missing(cons_mat)){
    cli_abort("Argument {.arg agg_mat} (or {.arg cons_mat}) is missing,
              with no default.", call = NULL)
  } else if(!missing(agg_mat)){
    tmp <- cstools(agg_mat = agg_mat)
    n <- tmp$dim[["n"]]
    strc_mat <- tmp$strc_mat
    agg_mat <- tmp$agg_mat
    cons_mat <- tmp$cons_mat
  } else {
    n <- NCOL(cons_mat)
    strc_mat <- NULL
    agg_mat <- NULL
    id_nn <- NULL
  }

  cov_mat <- cscov(comb = comb, n = n, agg_mat = agg_mat, res = res, ...)
  if(NROW(cov_mat) != n | NCOL(cov_mat) != n){
    cli_abort(c("Incorrect covariance dimensions.",
                "i"="Check {.arg res} columns dimension."), call = NULL)
  }

  projmat(cons_mat = cons_mat, cov_mat = cov_mat, strc_mat = strc_mat, mat = mat)
}

#' Projection matrix for optimal combination temporal reconciliation
#'
#' This function computes the projection or the mapping matrix
#' \eqn{\mathbf{M}} and \eqn{\mathbf{G}}, respectively, such that
#' \eqn{\widetilde{\mathbf{y}} = \mathbf{M}\widehat{\mathbf{y}} = \mathbf{S}_{te}\mathbf{G}\widehat{\mathbf{y}}},
#' where \eqn{\widetilde{\mathbf{y}}} is the vector of the reconciled forecasts,
#' \eqn{\widehat{\mathbf{y}}} is the vector of the base forecasts,
#' \eqn{\mathbf{S}_{te}} is the temporal structural matrix, and \eqn{\mathbf{M} = \mathbf{S}_{te}\mathbf{G}}.
#' For further information regarding on the structure of these matrices,
#' refer to Girolimetto et al. (2023).
#'
#' @inheritParams terec
#' @inheritParams csprojmat
#' @inheritDotParams tecov mse shrink_fun
#'
#' @inherit csprojmat return references
#'
#' @examples
#' # Temporal framework (annual-quarterly)
#' Mte <- teprojmat(agg_order = 4, comb = "ols")
#' Gte <- teprojmat(agg_order = 4, comb = "ols", mat = "G")
#'
#' @family Utilities
#' @export
teprojmat <- function(agg_order, comb = "ols", res = NULL, mat = "M", tew = "sum", ...){
  mat <- match.arg(mat, c("M", "G"))

  if(missing(agg_order)){
    cli_abort("Argument {.arg agg_order} is missing, with no default.", call = NULL)
  }

  tmp <- tetools(agg_order = agg_order, tew = tew)
  kset <- tmp$set
  m <- tmp$dim[["m"]]
  kt <- tmp$dim[["kt"]]

  # matrix
  agg_mat <- tmp$agg_mat
  strc_mat <- tmp$strc_mat
  cons_mat <- tmp$cons_mat

  cov_mat <- tecov(comb = comb, res = res, agg_order = agg_order, tew = tew, ...)
  if(NROW(cov_mat) != kt | NCOL(cov_mat) != kt){
    cli_abort(c("Incorrect covariance dimensions.",
                "i"="Check {.arg res} length."), call = NULL)
  }

  projmat(cons_mat = cons_mat, cov_mat = cov_mat, strc_mat = strc_mat, mat = mat)
}


#' Projection matrix for optimal combination cross-temporal reconciliation
#'
#' This function computes the projection or the mapping matrix
#' \eqn{\mathbf{M}} and \eqn{\mathbf{G}}, respectively, such that
#' \eqn{\widetilde{\mathbf{y}} = \mathbf{M}\widehat{\mathbf{y}} = \mathbf{S}_{ct}\mathbf{G}\widehat{\mathbf{y}}},
#' where \eqn{\widetilde{\mathbf{y}}} is the vector of the reconciled forecasts,
#' \eqn{\widehat{\mathbf{y}}} is the vector of the base forecasts,
#' \eqn{\mathbf{S}_{ct}} is the cross-temporal structural matrix, and \eqn{\mathbf{M} = \mathbf{S}_{ct}\mathbf{G}}.
#' For further information regarding on the structure of these matrices,
#' refer to Girolimetto et al. (2023).
#'
#' @usage
#' ctprojmat(agg_mat, cons_mat, agg_order, comb = "ols", res = NULL,
#'           mat = "M", tew = "sum", ...)
#'
#' @inheritParams ctrec
#' @inheritParams csprojmat
#' @inheritDotParams ctcov mse shrink_fun
#'
#' @inherit csprojmat return references
#'
#' @examples
#' # Cross-temporal framework (Z = X + Y, annual-quarterly)
#' A <- t(c(1,1)) # Aggregation matrix for Z = X + Y
#' Mct <- ctprojmat(agg_mat = A, agg_order = 4, comb = "ols")
#' Gct <- ctprojmat(agg_mat = A, agg_order = 4, comb = "ols", mat = "G")
#'
#' @family Utilities
#' @export
ctprojmat <- function(agg_mat, cons_mat, agg_order, comb = "ols", res = NULL, mat = "M", tew = "sum", ...){
  mat <- match.arg(mat, c("M", "G"))

  # Check if 'agg_order' is provided
  if(missing(agg_order)){
    cli_abort("Argument {.arg agg_order} is missing, with no default.", call = NULL)
  }

  # Check if either 'agg_mat' or 'cons_mat' is specified
  if(missing(agg_mat) && missing(cons_mat)){
    cli_abort("Argument {.arg agg_mat} (or {.arg cons_mat}) is missing,
              with no default.", call = NULL)
  } else if(!missing(agg_mat)){
    tmp <- cttools(agg_mat = agg_mat, agg_order = agg_order, tew = tew)
    strc_mat <- tmp$strc_mat
    cons_mat <- tmp$cons_mat
  } else {
    tmp <- cttools(cons_mat = cons_mat, agg_order = agg_order, tew = tew)
    strc_mat <- tmp$strc_mat
    cons_mat <- tmp$cons_mat
  }

  # Compute covariance
  cov_mat <- ctcov(comb = comb, res = res, agg_order = agg_order, agg_mat = agg_mat,
                   n = tmp$dim[["n"]], tew = tew, ...)
  if(NROW(cov_mat) != prod(tmp$dim[c("kt", "n")]) | NCOL(cov_mat) != prod(tmp$dim[c("kt", "n")])){
    cli_abort(c("Incorrect covariance dimensions.",
                "i"="Check {.arg res} dimensions."), call = NULL)
  }

  projmat(cons_mat = cons_mat, cov_mat = cov_mat, strc_mat = strc_mat, mat = mat)
}

projmat <- function(cons_mat, cov_mat, strc_mat, mat = "M"){
  mat <- match.arg(mat, c("M", "G"))

  if(mat == "M"){
    slm <- lin_sys(cons_mat %*% cov_mat %*% t(cons_mat), cons_mat)
    M <- unname(.sparseDiagonal(NCOL(cov_mat)) - cov_mat %*% t(cons_mat) %*% slm)
    return(Matrix(M))
  }else{
    if(!is.null(strc_mat)){
      if(isDiagonal(cov_mat)){
        cov_mat_inv <- .sparseDiagonal(x = diag(cov_mat)^(-1))
        lm_sx1 <- t(strc_mat) %*% cov_mat_inv %*% strc_mat
        lm_dx1 <- methods::as(t(strc_mat) %*% cov_mat_inv, "CsparseMatrix")
        G <- unname(lin_sys(lm_sx1, lm_dx1))
      }else{
        Q <- lin_sys(cov_mat, strc_mat)
        lm_sx1 <- t(strc_mat) %*% Q
        G <- unname(lin_sys(lm_sx1, t(Q)))
      }
      return(Matrix(G))
    }else{
      cli_abort("Argument {.arg agg_mat} is missing, with no default.", call = NULL)
    }
  }
}
