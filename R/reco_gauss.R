#' Cross-sectional Gaussian probabilistic reconciliation
#'
#' This function performs cross-sectional probabilistic forecast reconciliation
#' assuming a multivariate normal base forecast distribution (Panagiotelis et
#' al., 2023; Girolimetto et al., 2024; Wickramasuriya, 2024) for linearly
#' constrained (e.g. hierarchical or grouped) multiple time series.
#'
#' @usage
#' csmvn(base, agg_mat, cons_mat, comb = "ols", comb_base = comb,
#'       res = NULL, approach = "proj", reduce_form = FALSE, ...)
#'
#' @inheritParams csrec
#' @param comb_base A string specifying the reconciliation method. For a
#' complete list, see [cscov].
#' @param approach A string specifying the approach used to compute the
#' reconciled mean and covariance matrix. Options include:
#'   \itemize{
#'   \item "\code{proj}" (\emph{default}): Projection approach according to
#'   Byron (1978, 1979).
#'   \item "\code{strc}": Structural approach as proposed by Hyndman et
#'   al. (2011).
#'   }
#' @param reduce_form A logical parameter indicating whether the function
#' should return the full distribution (\code{FALSE}, \emph{default}) or only
#' the distribution corresponding to the bottom-level time series (\code{TRUE}).
#' @inheritDotParams cscov mse shrink_fun
#'
#' @returns A [distributional::dist_multivariate_normal] object.
#'
#' @references
#' Byron, R.P. (1978), The estimation of large social account matrices,
#' \emph{Journal of the Royal Statistical Society, Series A}, 141, 3, 359-367.
#' \doi{10.2307/2344807}
#'
#' Byron, R.P. (1979), Corrigenda: The estimation of large social account
#' matrices, \emph{Journal of the Royal Statistical Society, Series A}, 142(3),
#' 405. \doi{10.2307/2982515}
#'
#' Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J. (2024),
#' Cross-temporal probabilistic forecast reconciliation: Methodological and
#' practical issues. \emph{International Journal of Forecasting}, 40, 3,
#' 1134-1151. \doi{10.1016/j.ijforecast.2023.10.003}
#'
#' Hyndman, R.J., Ahmed, R.A., Athanasopoulos, G. and Shang, H.L. (2011),
#' Optimal combination forecasts for hierarchical time series,
#' \emph{Computational Statistics & Data Analysis}, 55, 9, 2579-2589.
#' \doi{10.1016/j.csda.2011.03.006}
#'
#' Panagiotelis, A., Gamakumara, P., Athanasopoulos, G. and Hyndman, R.J.
#' (2023), Probabilistic forecast reconciliation: Properties, evaluation and
#' score optimisation, \emph{European Journal of Operational Research} 306(2),
#' 693–706. \doi{10.1016/j.ejor.2022.07.040}
#'
#' Wickramasuriya, S. L. (2024). Probabilistic Forecast Reconciliation under
#' the Gaussian Framework. \emph{Journal of Business & Economic Statistics},
#' 42(1), 272–285. \doi{10.1080/07350015.2023.2181176}
#'
#' @examples
#' set.seed(123)
#' # (2 x 3) base forecasts matrix (simulated), Z = X + Y
#' base <- matrix(rnorm(6, mean = c(20, 10, 10)), 2, byrow = TRUE)
#' # (10 x 3) in-sample residuals matrix (simulated)
#' res <- t(matrix(rnorm(n = 30), nrow = 3))
#'
#' # Aggregation matrix for Z = X + Y
#' A <- t(c(1,1))
#' reco_dist <- csmvn(base = base, agg_mat = A, comb = "shr", res = res)
#'
#' @family Reco: probabilistic
#' @family Framework: cross-sectional
#'
#' @export
csmvn <- function(
  base,
  agg_mat,
  cons_mat,
  comb = "ols",
  comb_base = comb,
  res = NULL,
  approach = "proj",
  reduce_form = FALSE,
  ...
) {
  # Check if either 'agg_mat' or 'cons_mat' is specified
  if (missing(agg_mat) && missing(cons_mat)) {
    cli_abort(
      "Argument {.arg agg_mat} (or {.arg cons_mat}) is missing,
              with no default.",
      call = NULL
    )
  } else if (!missing(agg_mat)) {
    tmp <- cstools(agg_mat = agg_mat)
  } else {
    tmp <- cstools(cons_mat = cons_mat)
  }

  n <- tmp$dim[["n"]]
  strc_mat <- tmp$strc_mat
  cons_mat <- tmp$cons_mat

  if (!is.null(tmp$agg_mat) && reduce_form) {
    idts <- c(rep(0, tmp$dim[["na"]]), rep(1, tmp$dim[["nb"]]))
  } else {
    idts <- NULL
  }

  # Check if 'base' is provided and its dimensions match with the data
  if (missing(base)) {
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  } else if (NCOL(base) == 1) {
    base <- t(base)
  }

  if (NCOL(base) != n) {
    cli_abort("Incorrect {.arg base} columns dimension.", call = NULL)
  }

  # Compute covariance for reconciliation
  if (is(comb, "Matrix") | is(comb, "matrix")) {
    cov_mat <- comb
  } else {
    cov_mat <- cscov(
      comb = comb,
      n = n,
      agg_mat = agg_mat,
      res = res,
      strc_mat = strc_mat,
      ...
    )
  }

  if (NROW(cov_mat) != n | NCOL(cov_mat) != n) {
    cli_abort(
      c(
        "Incorrect covariance dimensions.",
        "i" = "Check {.arg res} columns dimension."
      ),
      call = NULL
    )
  }

  # Compute covariance base forecasts
  if (is(comb_base, "Matrix") | is(comb_base, "matrix")) {
    cov_mat_base <- comb_base
  } else if (comb_base != comb) {
    cov_mat_base <- cscov(
      comb = comb_base,
      n = n,
      agg_mat = agg_mat,
      res = res,
      strc_mat = strc_mat,
      ...
    )
  } else {
    comb_base <- comb
    cov_mat_base <- NULL
  }

  if (!is.null(cov_mat_base)) {
    if (NROW(cov_mat_base) != n | NCOL(cov_mat_base) != n) {
      cli_abort(
        c(
          "Incorrect base covariance dimensions.",
          "i" = "Check {.arg res} columns dimension."
        ),
        call = NULL
      )
    }
  }

  reco_mat <- pgreco(
    base = base,
    cov_mat = cov_mat,
    cov_mat_base = cov_mat_base,
    strc_mat = strc_mat,
    cons_mat = cons_mat,
    approach = approach,
    idts = idts
  )

  if (missing(agg_mat)) {
    tmp_name <- namesCS(n = n, names_vec = colnames(base))
  } else {
    tmp_name <- namesCS(
      n = n,
      names_vec = colnames(base),
      names_list = dimnames(agg_mat)
    )
  }

  if (!is.null(idts)) {
    tmp_name <- tmp_name[idts == 1]
  }

  colnames(reco_mat$mu) <- tmp_name
  rownames(reco_mat$mu) <- paste0("h-", 1:NROW(reco_mat$mu))

  reco_dist <- distributional::dist_multivariate_normal(
    mu = split(reco_mat$mu, 1:NROW(reco_mat$mu)),
    sigma = list(as.matrix(reco_mat$sigma))
  )
  dimnames(reco_dist) <- colnames(reco_mat$mu)
  names(reco_dist) <- paste0("h-", 1:length(reco_dist))

  attr(reco_dist, "FoReco") <- new_foreco_info(list(
    framework = "Cross-sectional",
    forecast_horizon = NROW(reco_mat),
    comb = comb,
    comb_base = comb_base,
    cs_n = n,
    rfun = "csrec",
    pfun = "csmvn"
  ))
  return(reco_dist)
}

#' Temporal Gaussian probabilistic reconciliation
#'
#' This function performs temporal probabilistic forecast reconciliation
#' assuming a multivariate normal base forecast distribution (Girolimetto et
#' al., 2024) for a single time series using temporal hierarchies
#' (Athanasopoulos et al., 2017).
#'
#' @usage
#' temvn(base, agg_order, comb = "ols", comb_base = comb, res = NULL,
#'       tew = "sum", approach = "proj", reduce_form = FALSE, ...)
#'
#' @param comb_base A string specifying the reconciliation method. For a
#' complete list, see [tecov].
#' @param reduce_form A logical parameter indicating whether the function
#' should return the full distribution (\code{FALSE}, \emph{default}) or only
#' the distribution corresponding to the high-frequency time series
#' (\code{TRUE}).
#' @inheritParams terec
#' @inheritParams csmvn
#' @inheritDotParams tecov mse shrink_fun
#'
#' @returns A [distributional::dist_multivariate_normal] object.
#'
#' @references
#' Athanasopoulos, G., Hyndman, R.J., Kourentzes, N. and Petropoulos, F. (2017),
#' Forecasting with Temporal Hierarchies, \emph{European Journal of Operational
#' Research}, 262, 1, 60-74. \doi{10.1016/j.ejor.2017.02.046}
#'
#' Byron, R.P. (1978), The estimation of large social account matrices,
#' \emph{Journal of the Royal Statistical Society, Series A}, 141, 3, 359-367.
#' \doi{10.2307/2344807}
#'
#' Byron, R.P. (1979), Corrigenda: The estimation of large social account
#' matrices, \emph{Journal of the Royal Statistical Society, Series A}, 142(3),
#' 405. \doi{10.2307/2982515}
#'
#' Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J. (2024),
#' Cross-temporal probabilistic forecast reconciliation: Methodological and
#' practical issues. \emph{International Journal of Forecasting}, 40, 3,
#' 1134-1151. \doi{10.1016/j.ijforecast.2023.10.003}
#'
#' Hyndman, R.J., Ahmed, R.A., Athanasopoulos, G. and Shang, H.L. (2011),
#' Optimal combination forecasts for hierarchical time series,
#' \emph{Computational Statistics & Data Analysis}, 55, 9, 2579-2589.
#' \doi{10.1016/j.csda.2011.03.006}
#'
#' @examples
#' set.seed(123)
#' # (7 x 1) base forecasts vector (simulated), m = 4
#' base <- rnorm(7*2, rep(c(20, 10, 5), 2*c(1, 2, 4)))
#' # (70 x 1) in-sample residuals vector (simulated)
#' res <- rnorm(70)
#'
#' m <- 4 # from quarterly to annual temporal aggregation
#' reco_dist <- terec(base = base, agg_order = m, comb = "wlsv", res = res)
#'
#' @family Reco: probabilistic
#' @family Framework: temporal
#'
#' @export
temvn <- function(
  base,
  agg_order,
  comb = "ols",
  comb_base = comb,
  res = NULL,
  tew = "sum",
  approach = "proj",
  reduce_form = FALSE,
  ...
) {
  # Check if 'agg_order' is provided
  if (missing(agg_order)) {
    cli_abort(
      "Argument {.arg agg_order} is missing, with no default.",
      call = NULL
    )
  }

  tmp <- tetools(agg_order = agg_order, tew = tew)
  kset <- tmp$set
  m <- tmp$dim[["m"]]
  kt <- tmp$dim[["kt"]]

  if (reduce_form) {
    idts <- c(rep(0, tmp$dim[["ks"]]), rep(1, m))
  } else {
    idts <- NULL
  }

  strc_mat <- tmp$strc_mat
  cons_mat <- tmp$cons_mat

  # Check if 'base' is provided and its dimensions match with the data
  if (missing(base)) {
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  } else if (NCOL(base) != 1) {
    cli_abort("{.arg base} is not a vector.", call = NULL)
  }

  # Calculate 'h' and 'base_hmat'
  if (length(base) %% kt != 0) {
    cli_abort("Incorrect {.arg base} length.", call = NULL)
  } else {
    h <- length(base) / kt
    base <- vec2hmat(vec = base, h = h, kset = kset)
  }

  # Compute covariance for reconciliation
  if (is(comb, "Matrix") | is(comb, "matrix")) {
    cov_mat <- comb
  } else {
    cov_mat <- tecov(
      comb = comb,
      res = res,
      agg_order = kset,
      tew = tew,
      strc_mat = strc_mat,
      ...
    )
  }

  if (NROW(cov_mat) != kt | NCOL(cov_mat) != kt) {
    cli_abort(
      c("Incorrect covariance dimensions.", "i" = "Check {.arg res} length."),
      call = NULL
    )
  }

  # Compute covariance base forecasts
  if (is(comb_base, "Matrix") | is(comb_base, "matrix")) {
    cov_mat_base <- comb_base
  } else if (comb_base != comb) {
    cov_mat_base <- tecov(
      comb = comb_base,
      res = res,
      agg_order = kset,
      tew = tew,
      strc_mat = strc_mat,
      ...
    )
  } else {
    comb_base <- comb
    cov_mat_base <- NULL
  }

  if (!is.null(cov_mat_base)) {
    if (NROW(cov_mat_base) != kt | NCOL(cov_mat_base) != kt) {
      cli_abort(
        c(
          "Incorrect base covariance dimensions.",
          "i" = "Check {.arg res} columns dimension."
        ),
        call = NULL
      )
    }
  }

  reco_mat <- pgreco(
    base = base,
    cov_mat = cov_mat,
    cov_mat_base = cov_mat_base,
    strc_mat = strc_mat,
    cons_mat = cons_mat,
    approach = approach,
    idts = idts
  )

  tmp_name <- namesTE(kset = kset, h = 1)
  if (!is.null(idts)) {
    tmp_name <- tmp_name[idts == 1]
  }
  colnames(reco_mat$mu) <- tmp_name

  reco_dist <- distributional::dist_multivariate_normal(
    mu = split(reco_mat$mu, 1:NROW(reco_mat$mu)),
    sigma = list(as.matrix(reco_mat$sigma))
  )
  dimnames(reco_dist) <- colnames(reco_mat$mu)
  names(reco_dist) <- paste0("tao-", 1:length(reco_dist))

  attr(reco_dist, "FoReco") <- new_foreco_info(list(
    framework = "Temporal",
    forecast_horizon = NROW(reco_mat),
    comb = comb,
    comb_base = comb_base,
    te_set = tmp$set,
    rfun = "terec",
    pfun = "temvn"
  ))
  return(reco_dist)
}

#' Cross-temporal Gaussian probabilistic reconciliation
#'
#' This function performs cross-temporal probabilistic forecast reconciliation
#' assuming a multivariate normal base forecast distribution (Girolimetto et
#' al., 2024) for linearly constrained multiple time series observed across both
#' cross-sectional and temporal dimensions (Di Fonzo and Girolimetto, 2023).
#'
#' @usage
#' ctmvn(base, agg_order, agg_mat, cons_mat, comb = "ols", comb_base = comb,
#'       res = NULL, tew = "sum", approach = "proj", reduce_form = FALSE, ...)
#'
#' @param comb_base A string specifying the reconciliation method. For a complete list, see [ctcov].
#' @param reduce_form A logical parameter indicating whether the function should return the
#' full distribution (\code{FALSE}, \emph{default}) or only the distribution corresponding
#' to the high-frequency bottom time series (\code{TRUE}).
#' @inheritParams ctrec
#' @inheritParams csmvn
#' @inheritDotParams ctcov mse shrink_fun
#'
#' @returns A [distributional::dist_multivariate_normal] object.
#'
#' @references
#' Byron, R.P. (1978), The estimation of large social account matrices,
#' \emph{Journal of the Royal Statistical Society, Series A}, 141, 3, 359-367.
#' \doi{10.2307/2344807}
#'
#' Byron, R.P. (1979), Corrigenda: The estimation of large social account matrices,
#' \emph{Journal of the Royal Statistical Society, Series A}, 142(3), 405.
#' \doi{10.2307/2982515}
#'
#' Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J. (2024),
#' Cross-temporal probabilistic forecast reconciliation: Methodological and
#' practical issues. \emph{International Journal of Forecasting}, 40, 3, 1134-1151.
#' \doi{10.1016/j.ijforecast.2023.10.003}
#'
#' Hyndman, R.J., Ahmed, R.A., Athanasopoulos, G. and Shang, H.L. (2011),
#' Optimal combination forecasts for hierarchical time series,
#' \emph{Computational Statistics & Data Analysis}, 55, 9, 2579-2589.
#' \doi{10.1016/j.csda.2011.03.006}
#'
#' Panagiotelis, A., Gamakumara, P., Athanasopoulos, G. and Hyndman, R.J. (2023),
#' Probabilistic forecast reconciliation: Properties, evaluation and score optimisation,
#' \emph{European Journal of Operational Research} 306(2), 693–706.
#' \doi{http://dx.doi.org/10.1016/j.ejor.2022.07.040}
#'
#' @examples
#' set.seed(123)
#' # (3 x 7) base forecasts matrix (simulated), Z = X + Y and m = 4
#' base <- rbind(rnorm(7, rep(c(20, 10, 5), c(1, 2, 4))),
#'               rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))),
#'               rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))))
#' # (3 x 70) in-sample residuals matrix (simulated)
#' res <- rbind(rnorm(70), rnorm(70), rnorm(70))
#' A <- t(c(1,1))
#' reco_dist <- ctmvn(base = base, res = res, agg_mat = A, agg_order = 4)
#'
#' @family Reco: probabilistic
#' @family Framework: cross-temporal
#'
#' @export
ctmvn <- function(
  base,
  agg_mat,
  cons_mat,
  agg_order,
  comb = "ols",
  comb_base = comb,
  res = NULL,
  tew = "sum",
  approach = "proj",
  reduce_form = FALSE,
  ...
) {
  # Check if 'base' is provided and its dimensions match with the data
  if (missing(base)) {
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }

  # Check if 'agg_order' is provided
  if (missing(agg_order)) {
    cli_abort(
      "Argument {.arg agg_order} is missing, with no default.",
      call = NULL
    )
  }

  # Check if either 'agg_mat' or 'cons_mat' is specified
  if (missing(agg_mat) && missing(cons_mat)) {
    cli_abort(
      "Argument {.arg agg_mat} (or {.arg cons_mat}) is missing,
              with no default.",
      call = NULL
    )
  } else if (!missing(agg_mat)) {
    tmp <- cttools(agg_mat = agg_mat, agg_order = agg_order, tew = tew)
    strc_mat <- tmp$strc_mat
    cons_mat <- tmp$cons_mat
  } else {
    tmp <- cttools(cons_mat = cons_mat, agg_order = agg_order, tew = tew)
    strc_mat <- tmp$strc_mat
    cons_mat <- tmp$cons_mat
    agg_mat <- cstools(cons_mat = cons_mat)$agg_mat
  }

  if (!is.null(tmp$agg_mat) && reduce_form) {
    cs_nn <- c(rep(0, tmp$dim[["na"]]), rep(1, tmp$dim[["nb"]]))
    te_nn <- c(rep(0, tmp$dim[["ks"]]), rep(1, tmp$dim[["m"]]))
    idts <- as.numeric(kronecker(cs_nn, te_nn))
  } else {
    idts <- NULL
    reduce_form <- FALSE
  }

  if (NCOL(base) %% tmp$dim[["kt"]] != 0) {
    cli_abort("Incorrect {.arg base} columns dimension.", call = NULL)
  }

  if (NROW(base) != tmp$dim[["n"]]) {
    cli_abort("Incorrect {.arg base} rows dimension.", call = NULL)
  }

  # Calculate 'h' and 'base_hmat'
  h <- NCOL(base) / tmp$dim[["kt"]]
  base_hmat <- mat2hmat(base, h = h, kset = tmp$set, n = tmp$dim[["n"]])

  # Compute covariance for reconciliation
  if (is(comb, "Matrix") | is(comb, "matrix")) {
    cov_mat <- comb
  } else {
    cov_mat <- ctcov(
      comb = comb,
      res = res,
      agg_order = tmp$set,
      agg_mat = agg_mat,
      n = tmp$dim[["n"]],
      tew = tew,
      strc_mat = strc_mat,
      ...
    )
  }
  if (
    NROW(cov_mat) != prod(tmp$dim[c("kt", "n")]) |
      NCOL(cov_mat) != prod(tmp$dim[c("kt", "n")])
  ) {
    cli_abort(
      c(
        "Incorrect covariance dimensions.",
        "i" = "Check {.arg res} dimensions."
      ),
      call = NULL
    )
  }

  # Compute covariance base forecasts
  if (is(comb_base, "Matrix") | is(comb_base, "matrix")) {
    cov_mat_base <- comb_base
  } else if (comb_base != comb) {
    cov_mat_base <- ctcov(
      comb = comb_base,
      res = res,
      agg_order = tmp$set,
      agg_mat = agg_mat,
      n = tmp$dim[["n"]],
      tew = tew,
      strc_mat = strc_mat,
      ...
    )
  } else {
    comb_base <- comb
    cov_mat_base <- NULL
  }

  if (!is.null(cov_mat_base)) {
    if (
      NROW(cov_mat_base) != prod(tmp$dim[c("kt", "n")]) |
        NCOL(cov_mat_base) != prod(tmp$dim[c("kt", "n")])
    ) {
      cli_abort(
        c(
          "Incorrect covariance dimensions.",
          "i" = "Check {.arg res} dimensions."
        ),
        call = NULL
      )
    }
  }

  reco_mat <- pgreco(
    base = base_hmat,
    cov_mat = cov_mat,
    cov_mat_base = cov_mat_base,
    strc_mat = strc_mat,
    cons_mat = cons_mat,
    approach = approach,
    idts = idts
  )
  nsize <- ifelse(reduce_form, tmp$dim[["nb"]], tmp$dim[["n"]])
  names_cs <- paste0("s-", 1:nsize)
  names_te <- namesTE(kset = tmp$set, h = 1)
  names_ct <- paste(
    rep(names_cs, each = length(names_te)),
    rep(names_te, length(names_cs))
  )

  reco_dist <- distributional::dist_multivariate_normal(
    mu = split(reco_mat$mu, 1:NROW(reco_mat$mu)),
    sigma = list(as.matrix(reco_mat$sigma))
  )
  dimnames(reco_dist) <- names_ct
  names(reco_dist) <- paste0("tao-", 1:length(reco_dist))

  attr(reco_dist, "FoReco") <- new_foreco_info(list(
    info = attr(reco_mat, "info"),
    framework = "Cross-temporal",
    forecast_horizon = h,
    comb = comb,
    te_set = tmp$set,
    cs_n = tmp$dim[["n"]],
    rfun = "ctrec",
    pfun = "ctmvn"
  ))
  return(reco_dist)
}

pgreco <- function(
  approach,
  base,
  cons_mat,
  strc_mat,
  cov_mat,
  cov_mat_base,
  idts = idts,
  ...
) {
  if (approach == "proj") {
    # check input
    if (missing(base) | missing(cons_mat) | missing(cov_mat)) {
      cli_abort(
        "Mandatory arguments: {.arg base}, {.arg cons_mat} and {.arg cov_mat}.",
        call = NULL
      )
    }

    if (NCOL(cons_mat) != NROW(cov_mat) | NCOL(base) != NROW(cov_mat)) {
      cli_abort("The size of the matrices does not match.", call = NULL)
    }

    # Point reconciled forecasts
    lm_dx <- methods::as(Matrix::tcrossprod(cons_mat, base), "CsparseMatrix")
    lm_sx <- methods::as(
      Matrix::tcrossprod(cons_mat %*% cov_mat, cons_mat),
      "CsparseMatrix"
    )
    lm_sol <- Matrix::crossprod(cons_mat, lin_sys(lm_sx, cons_mat))

    if (!is.null(idts)) {
      idts <- idts == 1
      reco <- base[, idts, drop = FALSE] -
        t(cov_mat[idts, , drop = FALSE] %*% Matrix::tcrossprod(lm_sol, base))
      M <- unname(
        .sparseDiagonal(NCOL(cov_mat))[idts, , drop = FALSE] -
          cov_mat[idts, , drop = FALSE] %*% lm_sol
      )
    } else {
      reco <- base - t(cov_mat %*% Matrix::tcrossprod(lm_sol, base))
      M <- unname(.sparseDiagonal(NCOL(cov_mat)) - cov_mat %*% lm_sol)
    }

    if (is.null(cov_mat_base)) {
      if (!is.null(idts)) {
        covr <- M %*% cov_mat[, idts, drop = FALSE]
      } else {
        covr <- M %*% cov_mat
      }
    } else {
      covr <- M %*% Matrix::tcrossprod(cov_mat_base, M)
    }
  } else if (approach == "strc") {
    # check input
    if (missing(base) | missing(strc_mat) | missing(cov_mat)) {
      cli_abort(
        "Mandatory arguments: {.arg base}, {.arg strc_mat} and {.arg cov_mat}.",
        call = NULL
      )
    }

    if (is.null(strc_mat)) {
      cli_abort(
        "Please provide a valid {.arg agg_mat} for the structural approach.",
        call = NULL
      )
    }

    if (NROW(strc_mat) != NROW(cov_mat) | NCOL(base) != NROW(cov_mat)) {
      cli_abort("The size of the matrices does not match.", call = NULL)
    }

    # Point reconciled forecasts
    if (isDiagonal(cov_mat)) {
      cov_mat_inv <- .sparseDiagonal(x = diag(cov_mat)^(-1))
      StWm <- Matrix::crossprod(strc_mat, cov_mat_inv)
      lm_sx1 <- methods::as(StWm %*% strc_mat, "CsparseMatrix")
      lm_dx1 <- methods::as(Matrix::tcrossprod(StWm, base), "CsparseMatrix")

      if (!is.null(idts)) {
        idts <- idts == 1
        reco <- t(lin_sys(lm_sx1, lm_dx1))
        if (is.null(cov_mat_base)) {
          lm_dx1_cov <- methods::as(
            Matrix::tcrossprod(StWm, cov_mat[idts, , drop = FALSE]),
            "CsparseMatrix"
          )
          covr <- lin_sys(lm_sx1, lm_dx1_cov)
        } else {
          G <- lin_sys(lm_sx1, StWm)
          covr <- G %*% Matrix::tcrossprod(cov_mat_base, G)
        }
      } else {
        reco <- t(strc_mat %*% lin_sys(lm_sx1, lm_dx1))
        if (is.null(cov_mat_base)) {
          lm_dx1_cov <- methods::as(
            Matrix::tcrossprod(StWm, cov_mat),
            "CsparseMatrix"
          )
          covr <- strc_mat %*% lin_sys(lm_sx1, lm_dx1_cov)
        } else {
          SG <- strc_mat %*% lin_sys(lm_sx1, StWm)
          covr <- SG %*% Matrix::tcrossprod(cov_mat_base, SG)
        }
      }
    } else {
      Q <- lin_sys(cov_mat, strc_mat)
      lm_sx1 <- methods::as(crossprod(strc_mat, Q), "CsparseMatrix")
      lm_dx1 <- methods::as(t(base %*% Q), "CsparseMatrix")

      if (!is.null(idts)) {
        idts <- idts == 1
        reco <- t(lin_sys(lm_sx1, lm_dx1))

        if (is.null(cov_mat_base)) {
          lm_dx1_cov <- methods::as(
            t(cov_mat[idts, , drop = FALSE] %*% Q),
            "CsparseMatrix"
          )
          covr <- lin_sys(lm_sx1, lm_dx1_cov)
        } else {
          G <- lin_sys(lm_sx1, t(Q))
          covr <- G %*% Matrix::tcrossprod(cov_mat_base, G)
        }
      } else {
        reco <- t(strc_mat %*% lin_sys(lm_sx1, lm_dx1))

        if (is.null(cov_mat_base)) {
          lm_dx1_cov <- methods::as(t(cov_mat %*% Q), "CsparseMatrix")
          covr <- strc_mat %*% lin_sys(lm_sx1, lm_dx1_cov)
        } else {
          SG <- strc_mat %*% lin_sys(lm_sx1, t(Q))
          covr <- SG %*% Matrix::tcrossprod(cov_mat_base, SG)
        }
      }
    }
  } else if (approach == "bu") {
    # check input
    if (missing(base) | missing(strc_mat) | missing(cov_mat)) {
      cli_abort(
        "Mandatory arguments: {.arg base}, {.arg strc_mat} and {.arg cov_mat}.",
        call = NULL
      )
    }

    if (is.null(strc_mat)) {
      cli_abort(
        "Please provide a valid {.arg agg_mat} for the structural approach.",
        call = NULL
      )
    }

    reco <- tcrossprod(base, strc_mat)
    covr <- strc_mat %*% tcrossprod(cov_mat, strc_mat)
  } else {
    cli_abort("No support.", call = NULL)
  }

  return(list(mu = as.matrix(reco), sigma = covr))
}
