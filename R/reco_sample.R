#' Cross-sectional probabilistic reconciliation (sample approach)
#'
#' @param sample A (\eqn{h \times n \times L}) numeric array containing the base forecasts
#' samples to be reconciled; \eqn{h} is the forecast horizon, \eqn{n} is the total number
#' of time series (\eqn{n = n_a + n_b}), and \eqn{L} is the sample size.
#' @param fun A string specifying the reconciliation function to be used,
#' as implemented in \pkg{FoReco}.
#' @param ... Arguments passed on to \code{fun}
#'
#' @returns A [distributional::dist_sample] object.
#'
#' @references
#' Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J. (2024),
#' Cross-temporal probabilistic forecast reconciliation: Methodological and
#' practical issues. \emph{International Journal of Forecasting}, 40, 3, 1134-1151.
#' \doi{10.1016/j.ijforecast.2023.10.003}
#'
#' Panagiotelis, A., Gamakumara, P., Athanasopoulos, G. and Hyndman, R.J. (2023),
#' Probabilistic forecast reconciliation: Properties, evaluation and score optimisation,
#' \emph{European Journal of Operational Research} 306(2), 693–706.
#' \doi{http://dx.doi.org/10.1016/j.ejor.2022.07.040}
#'
#' @examples
#' set.seed(123)
#' A <- matrix(c(1,1,1,1,  # Z = X + Y
#'               1,1,0,0,  # X = XX + XY
#'               0,0,1,1), # Y = YX + YY
#'               nrow = 3, byrow = TRUE)
#' rownames(A) <- c("Z", "X", "Y")
#' colnames(A) <- c("XX", "XY", "YX", "YY")
#'
#' # (100 x 7) base forecasts sample (simulated) for h = 1
#' base_h1 <- matrix(rnorm(100*7, mean = c(20, rep(10, 2), rep(5, 4))),
#'                   100, byrow = TRUE)
#' # (100 x 7) base forecasts sample (simulated) for h = 2
#' base_h2 <- matrix(rnorm(100*7, mean = c(20, rep(10, 2), rep(5, 4))),
#'                   100, byrow = TRUE)
#' # (2 x 7 x 100) base forecasts sample array with
#' # 2 forecast horizons, 7 time series and 100 sample
#' base_sample <- aperm(simplify2array(list(base_h1, base_h2)), c(3,2,1))
#'
#' # Top-down probabilistic reconciliation
#' reco_dist_td <- cssmp(base_sample[, 1, , drop = FALSE], agg_mat = A,
#'                       fun = cstd, weights = c(0.3, 0.2, 0.1, 0.4))
#'
#' # Middle-out probabilistic reconciliation
#' reco_dist_mo <- cssmp(base_sample[, c(2,3), , drop = FALSE], agg_mat = A,
#'                       fun = csmo, weights = c(0.3, 0.7, 0.8, 0.2),
#'                       id_rows = 2:3)
#'
#' # Bottom-up probabilistic reconciliation
#' reco_dist_bu <- cssmp(base_sample[,-c(1:3),], agg_mat = A, fun = csbu)
#'
#' # Level conditional coherent probabilistic reconciliation
#' reco_dist_lcc <- cssmp(base_sample, agg_mat = A, fun = cslcc)
#'
#' # Optimal cross-sectional probabilistic reconciliation
#' reco_dist_opt <- cssmp(base_sample, agg_mat = A)
#'
#' @family Reco: probabilistic
#' @family Framework: cross-sectional
#'
#' @export
cssmp <- function(sample, fun = csrec, ...) {
  if (is.list(sample)) {
    sample <- simplify2array(sample)
  }

  if (length(dim(sample)) == 2) {
    # Matrix input: sample_size x variable
    tmp_dim <- c(1, rev(dim(sample))) # forecast_horizon x variable x sample_size
  } else if (length(dim(sample)) == 3) {
    # Array input: forecast_horizon x variable x sample_size
    tmp_dim <- dim(sample) # forecast_horizon x variable x sample_size

    # array to matrix form
    sample <- aperm(sample, c(2, 3, 1))
    dim(sample) <- c(tmp_dim[2], prod(tmp_dim[-2]))
    sample <- t(sample)
  } else {
    cli_abort("Incorrect {.arg sample} dimensions", call = NULL)
  }

  # reconciliation
  reco_mat <- fun(base = sample, ...)

  colnames_save <- colnames(reco_mat)
  obj <- recoinfo(reco_mat, verbose = FALSE)
  obj$sample_size <- tmp_dim[3]
  obj$forecast_horizon <- obj$forecast_horizon / obj$sample_size
  obj$pfun <- "cssmp"

  # matrix to array
  tmp_dim[2] <- NCOL(reco_mat)
  reco_mat <- t(reco_mat)
  dim(reco_mat) <- tmp_dim[c(2, 3, 1)]
  reco_mat <- aperm(reco_mat, c(2, 1, 3))
  dimnames(reco_mat) <- list(
    paste0("l-", 1:tmp_dim[3]),
    colnames_save,
    paste0("h-", 1:tmp_dim[1])
  )

  reco_dist <- do.call(
    c,
    apply(reco_mat, 3, function(x) {
      dist <- distributional::dist_sample(list(x))
      dist
    })
  )
  dimnames(reco_dist) <- colnames_save
  names(reco_dist) <- paste0("h-", 1:tmp_dim[1])

  attr(reco_dist, "FoReco") <- new_foreco_info(obj)
  return(reco_dist)
}


#' Temporal probabilistic reconciliation (sample approach)
#'
#' @param sample A (\eqn{L \times h(k^\ast + m)}) numeric matrix containing the base forecasts
#' samples to be reconciled; \eqn{m} is the max aggregation order, \eqn{k^\ast} is the sum of
#' (a subset of) (\eqn{p-1}) factors of \eqn{m}, excluding \eqn{m}, \eqn{h} is the forecast
#' horizon for the lowest frequency time series, and \eqn{L} is the sample size.
#' @inheritParams cssmp
#' @inheritParams terec
#'
#' @returns A [distributional::dist_sample] object.
#'
#' @references
#' Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J. (2024),
#' Cross-temporal probabilistic forecast reconciliation: Methodological and
#' practical issues. \emph{International Journal of Forecasting}, 40, 3, 1134-1151.
#' \doi{10.1016/j.ijforecast.2023.10.003}
#'
#' Panagiotelis, A., Gamakumara, P., Athanasopoulos, G. and Hyndman, R.J. (2023),
#' Probabilistic forecast reconciliation: Properties, evaluation and score optimisation,
#' \emph{European Journal of Operational Research} 306(2), 693–706.
#' \doi{http://dx.doi.org/10.1016/j.ejor.2022.07.040}
#'
#' @examples
#' set.seed(123)
#' m <- 4 # from quarterly to annual temporal aggregation
#'
#' # (100 x 14) base forecasts sample matrix (simulated), m = 4, h = 2
#' sample <- t(sapply(1:100, function(x) rnorm(14, rep(c(20, 10, 5), 2*c(1, 2, 4)))))
#' # (70 x 1) in-sample residuals vector (simulated)
#' res <- rnorm(70)
#'
#' # Top-down probabilistic reconciliation
#' reco_dist_td <- tesmp(sample[,c(1:2), drop = FALSE], agg_order = m,
#'                       fun = tetd, weights = c(0.2, 0.5, 0.3, 0.3))
#'
#' # Middle-out probabilistic reconciliation
#' tesmp(sample[,c(3:6), drop = FALSE], agg_order = m,
#'       fun = temo, weights = c(0.2, 0.5, 0.3, 0.3), order = 2)
#'
#' # Bottom-up probabilistic reconciliation
#' reco_dist_bu <- tesmp(sample[,-c(1:6)], agg_order = m, fun = tebu)
#'
#' # Level conditional coherent probabilistic reconciliation
#' reco_dist_lcc <- tesmp(sample, agg_order = m, fun = telcc)
#'
#' # Optimal cross-sectional probabilistic reconciliation
#' reco_dist_opt <- tesmp(sample, agg_order = m, res = res, comb = "wlsv")
#'
#' @family Reco: probabilistic
#' @family Framework: temporal
#'
#' @export
tesmp <- function(sample, agg_order, fun = terec, ...) {
  if (length(dim(sample)) == 2) {
    # Matrix input: sample_size x (variable*h)
    sample_size <- NROW(sample)
    if (as.character(substitute(fun)) %in% c("temo", "tetd")) {
      sample <- as.vector(t(sample))
    } else {
      sample <- unname(unlist(lapply(
        FoReco2matrix(sample, agg_order),
        as.vector
      )))
    }
  } else if (is.vector(sample)) {
    # vector input: variable * sample_size x 1
    sample_size <- NULL
  } else {
    cli_abort("Incorrect {.arg sample} dimensions", call = NULL)
  }

  # reconciliation
  reco_mat <- fun(base = sample, agg_order = agg_order, ...)

  # info
  obj <- recoinfo(reco_mat, verbose = FALSE)
  if (is.null(sample_size)) {
    obj$sample_size <- obj$forecast_horizon
  } else {
    obj$sample_size <- sample_size
  }
  obj$forecast_horizon <- obj$forecast_horizon / obj$sample_size
  obj$pfun <- "tesmp"

  # Names
  names_te <- namesTE(kset = obj$te_set, h = 1)

  # Distributional obj
  reco_mat <- vec2hmat(
    reco_mat,
    obj$forecast_horizon * obj$sample_size,
    obj$te_set
  )
  reco_mat <- lapply(
    split(reco_mat, rep(1:obj$forecast_horizon, obj$sample_size)),
    matrix,
    nrow = obj$sample_size
  )
  reco_mat <- lapply(
    reco_mat,
    `dimnames<-`,
    list(paste0("l-", 1:obj$sample_size), names_te)
  )
  reco_dist <- dist_sample(reco_mat)
  dimnames(reco_dist) <- names_te
  names(reco_dist) <- paste0("tao-", 1:length(reco_dist))

  attr(reco_dist, "FoReco") <- new_foreco_info(obj)
  return(reco_dist)
}

#' Cross-temporal probabilistic reconciliation (sample approach)
#'
#' @param sample A (\eqn{n \times h(k^\ast + m) \times L}) numeric array containing the
#' base forecasts samples to be reconciled; \eqn{n} is the total number of variables,
#' \eqn{m} is the max. order of temporal aggregation, \eqn{k^\ast} is the sum of (a subset
#' of) (\eqn{p-1}) factors of \eqn{m},  excluding \eqn{m}, \eqn{h} is the forecast horizon
#' for the lowest frequency time series, and \eqn{L} is the sample size.
#' The row identifies a time series, and the forecasts in each row are ordered from the
#' lowest frequency (most temporally aggregated) to the highest frequency.
#' @inheritParams cssmp
#' @inheritParams ctrec
#'
#' @returns A [distributional::dist_sample] object.
#'
#' @references
#' Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J. (2024),
#' Cross-temporal probabilistic forecast reconciliation: Methodological and
#' practical issues. \emph{International Journal of Forecasting}, 40, 3, 1134-1151.
#' \doi{10.1016/j.ijforecast.2023.10.003}
#'
#' Panagiotelis, A., Gamakumara, P., Athanasopoulos, G. and Hyndman, R.J. (2023),
#' Probabilistic forecast reconciliation: Properties, evaluation and score optimisation,
#' \emph{European Journal of Operational Research} 306(2), 693–706.
#' \doi{http://dx.doi.org/10.1016/j.ejor.2022.07.040}
#'
#' @examples
#' set.seed(123)
#' A <- t(c(1,1)) # Aggregation matrix for Z = X + Y
#' m <- 4 # from quarterly to annual temporal aggregation
#'
#' # (100 x 14) base forecasts sample matrix (simulated), m = 4, h = 2, n = 3
#' sample <- simplify2array(lapply(1:100, function(x){
#'   rbind(rnorm(14, rep(c(20, 10, 5), 2*c(1, 2, 4))),
#'         rnorm(14, rep(c(10, 5, 2.5), 2*c(1, 2, 4))),
#'         rnorm(14, rep(c(10, 5, 2.5), 2*c(1, 2, 4))))
#' }))
#' # (3 x 70) in-sample residuals matrix (simulated)
#' res <- rbind(rnorm(70), rnorm(70), rnorm(70))
#'

#' # Top-down probabilistic reconciliation
#' reco_dist_td <- ctsmp(sample[1, 1:2, , drop = FALSE], agg_order = m,
#'                       agg_mat = A, fun = cttd, weights = matrix(runif(8), 2))
#'
#' # Middle-out probabilistic reconciliation
#' reco_dist_mo <- ctsmp(sample[1, 3:6, , drop = FALSE], agg_order = m,
#'                       agg_mat = A, fun = ctmo, weights = matrix(runif(8), 2),
#'                       id_rows = 1, order = 2)
#' # Bottom-up probabilistic reconciliation
#' reco_dist_bu <- ctsmp(sample[-1,-c(1:6), ], agg_order = m, agg_mat = A, fun = ctbu)
#'
#' # Level conditional coherent probabilistic reconciliation
#' reco_dist_lcc <- ctsmp(sample, agg_order = m, agg_mat = A, fun = ctlcc)
#'
#' # Optimal cross-sectional probabilistic reconciliation
#' reco_dist_opt <- ctsmp(sample, agg_order = m, agg_mat = A, res = res, comb = "bdshr")
#'
#' @family Reco: probabilistic
#' @family Framework: cross-temporal
#'
#' @export
ctsmp <- function(sample, agg_order, fun = ctrec, ...) {
  if (length(dim(sample)) == 3) {
    tmp_dim <- dim(sample)
    if (as.character(substitute(fun)) %in% c("ctmo", "cttd")) {
      dim(sample) <- c(tmp_dim[1], tmp_dim[2] * tmp_dim[3])
      if (tmp_dim[1] == 1) {
        sample <- as.vector(sample)
      }
    } else {
      rownames_save <- rownames(sample)
      sample <- apply(sample, 1, function(x) {
        (unlist(lapply(FoReco::FoReco2matrix(t(x), agg_order), as.vector)))
      })
      sample <- t(unname(sample))
      rownames(sample) <- rownames_save
    }
  } else {
    cli_abort("Incorrect {.arg sample} dimensions", call = NULL)
  }

  # reconciliation
  if ("agg_order" %in% names(formals(fun))) {
    reco_mat <- fun(base = sample, agg_order = agg_order, ...)
  } else {
    reco_mat <- fun(base = sample, ...)
  }
  rownames_save <- colnames(reco_mat)

  # Reco info
  obj <- recoinfo(reco_mat, verbose = FALSE)
  obj$sample_size <- tmp_dim[3]
  obj$forecast_horizon <- obj$forecast_horizon / obj$sample_size
  obj$pfun <- "ctsmp"

  # Names
  names_cs <- paste0("s-", 1:obj$cs_n)
  names_te <- namesTE(kset = obj$te_set, h = 1)
  names_ct <- paste(
    rep(names_cs, each = length(names_te)),
    rep(names_te, length(names_cs))
  )

  # Distribution obj
  reco_mat <- mat2hmat(
    reco_mat,
    obj$forecast_horizon * obj$sample_size,
    obj$te_set,
    obj$cs_n
  )
  reco_mat <- lapply(
    split(reco_mat, rep(1:obj$forecast_horizon, obj$sample_size)),
    matrix,
    nrow = obj$sample_size
  )
  reco_mat <- lapply(
    reco_mat,
    `dimnames<-`,
    list(paste0("l-", 1:obj$sample_size), names_ct)
  )
  reco_dist <- dist_sample(reco_mat)
  dimnames(reco_dist) <- names_ct
  names(reco_dist) <- paste0("tao-", 1:length(reco_dist))

  attr(reco_dist, "FoReco") <- new_foreco_info(obj)
  return(reco_dist)
}
