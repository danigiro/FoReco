#' Cross-sectional top-down reconciliation
#'
#' Top-down forecast reconciliation for genuine hierarchical/grouped time series
#' (Gross and Sohl, 1990), where the forecast of a `Total' (top-level series,
#' expected to be positive) is disaggregated according to a proportional scheme
#' (weights). Besides fulfilling any aggregation constraint, the top-down
#' reconciled forecasts should respect two main properties:
#' \itemize{
#'   \item the top-level value remains unchanged;
#'   \item all the bottom time series reconciled forecasts are non-negative.
#' }
#'
#' @param base A (\eqn{h \times 1}) numeric vector containing the top-level base forecast;
#' \eqn{h} is the forecast horizon.
#' @inheritParams csrec
#' @param weights A (\eqn{h \times n_b}) numeric matrix containing the proportions for the
#' bottom time series; \eqn{h} is the forecast horizon, and \eqn{n_b} is the total number
#' of bottom variables.
#' @param normalize If \code{TRUE} (\emph{default}), the \code{weights} will sum to 1.
#'
#' @references
#' Gross, C.W. and Sohl, J.E. (1990), Disaggregation methods to expedite product
#' line forecasting. \emph{Journal of Forecasting} 9(3), 233–254. \doi{10.1002/for.3980090304}
#'
#' @examples
#' set.seed(123)
#' # Aggregation matrix for Z = X + Y, X = XX + XY and Y = YX + YY
#' A <- matrix(c(1,1,1,1,1,1,0,0,0,0,1,1), 3, byrow = TRUE)
#' # (3 x 1) top base forecasts vector (simulated), forecast horizon = 3
#' topf <- rnorm(3, 10)
#' # Same weights for different forecast horizons
#' fix_weights <- runif(4)
#' reco <- cstd(base = topf, agg_mat = A, weights = fix_weights)
#'
#' # Different weights for different forecast horizons
#' h_weights <- matrix(runif(4*3), 3, 4)
#' recoh <- cstd(base = topf, agg_mat = A, weights = h_weights)
#'
#' @inherit csrec return
#'
#' @family Reco: top-down
#' @family Framework: cross-sectional
#' @export
#'
cstd <- function(base, agg_mat, weights, normalize = TRUE){
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }else if(NCOL(base) != 1){
    cli_abort("Argument {.arg base} is not a vector.", call = NULL)
  }

  if(missing(agg_mat)){
    cli_abort("Argument {.arg agg_mat} is missing, with no default.", call = NULL)
  }else{
    tmp <- cstools(agg_mat = agg_mat, sparse = TRUE)
  }

  if(missing(weights)){
    cli_abort("Argument {.arg weights} is missing, with no default.", call = NULL)
  }

  if(NCOL(weights) == 1){
    weights <- t(weights)
  }

  if(NCOL(weights)!=tmp$dim[["nb"]]){
    cli_abort("Incorrect {.arg agg_mat} or {.arg weights} dimensions.", call = NULL)
  }

  if(NROW(weights)!=length(base)){
    weights <- matrix(rep(as.vector(weights), length(base)), length(base), byrow = TRUE)
  }

  wsum <- apply(weights, 1, sum)

  if(!normalize){
    if(any(wsum != 1)){
      cli_warn("The {.arg weights} do not add up to 1", call = NULL)
    }
    wsum <- rep(1, length(base))
  }

  bbf <- weights*matrix(base/wsum, nrow = NROW(weights), ncol = NCOL(weights))
  reco_mat <- csbu(bbf, agg_mat = tmp$agg_mat)
  attr(reco_mat, "FoReco") <- list2env(list(framework = "Cross-sectional",
                                            forecast_horizon = NROW(reco_mat),
                                            cs_n = NCOL(reco_mat),
                                            rfun = "cstd"))
  return(reco_mat)
}

#' Temporal top-down reconciliation
#'
#' Top-down forecast reconciliation for a univariate time series, where the forecast
#' of the most aggregated temporal level is disaggregated according to a proportional
#' scheme (weights). Besides fulfilling any aggregation constraint, the
#' top-down reconciled forecasts should respect two main properties:
#' \itemize{
#'   \item the top-level value remains unchanged;
#'   \item all the bottom time series reconciled forecasts are non-negative.
#' }
#'
#' @param base A (\eqn{hm \times 1}) numeric vector containing the temporal
#' aggregated base forecasts of order \eqn{m}; \eqn{m} is the max aggregation
#' order, and \eqn{h} is the forecast horizon for the lowest frequency
#' time series.
#' @inheritParams terec
#' @param weights A (\eqn{hm \times 1}) numeric vector containing the proportions for the
#' high-frequency time series; \eqn{m} is the max aggregation order, and \eqn{h} is the
#' forecast horizon for the lowest frequency time series.
#' @inheritParams cstd
#'
#' @inherit terec return
#'
#' @examples
#' set.seed(123)
#' # (2 x 1) top base forecasts vector (simulated), forecast horizon = 2
#' topf <- rnorm(2, 10)
#' # Same weights for different forecast horizons
#' fix_weights <- runif(4)
#' reco <- tetd(base = topf, agg_order = 4, weights = fix_weights)
#'
#' # Different weights for different forecast horizons
#' h_weights <- runif(4*2)
#' recoh <- tetd(base = topf, agg_order = 4, weights = h_weights)
#'
#' @family Reco: top-down
#' @family Framework: temporal
#' @export
#'
tetd <- function(base, agg_order, weights, tew = "sum", normalize = TRUE){
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }else if(NCOL(base) != 1){
    cli_abort("Argument {.arg base} is not a vector.", call = NULL)
  }

  if(missing(agg_order)){
    cli_abort("Argument {.arg agg_order} is missing, with no default.", call = NULL)
  }else{
    tmp <- tetools(agg_order = agg_order, sparse = TRUE)
  }

  if(missing(weights)){
    cli_abort("Argument {.arg weights} is missing, with no default.", call = NULL)
  }else if(NCOL(weights) != 1){
    cli_abort("Argument {.arg weights} is not a vector.", call = NULL)
  }

  if(length(weights) != length(base)*tmp$dim[["m"]]){
    if(length(weights) == tmp$dim[["m"]]){
      weights <- rep(weights, length(base))
    }else{
      cli_abort("Incorrect {.arg weights} length.", call = NULL)
    }
  }

  wsum <- unname(tapply(weights, rep(1:length(base), each = tmp$dim[["m"]]), sum))

  if(!normalize){
    if(any(wsum != 1)){
      cli_warn("The {.arg weights} do not add up to 1", call = NULL)
    }
    wsum <- rep(1, length(base))
  }

  hfbf <- rep(base/wsum, each = tmp$dim[["m"]])*weights
  reco_vec <- tebu(hfbf, agg_order = tmp$set, tew = tew)
  attr(reco_vec, "FoReco") <- list2env(list(framework = "Temporal",
                                            forecast_horizon = length(base),
                                            te_set = tmp$set,
                                            rfun = "tetd"))
  return(reco_vec)
}

#' Cross-temporal top-down reconciliation
#'
#' Top-down forecast reconciliation for cross-temporal hierarchical/grouped
#' time series, where the forecast of a `Total' (top-level series, expected
#' to be positive) is disaggregated according to a proportional scheme (weights).
#' Besides fulfilling any aggregation constraint, the top-down reconciled
#' forecasts should respect two main properties:
#' \itemize{
#'   \item the top-level value remains unchanged;
#'   \item all the bottom time series reconciled forecasts are non-negative.
#' }
#'
#' @param base A (\eqn{hm \times 1}) numeric vector containing top- and \eqn{m} temporal
#' aggregated level base forecasts; \eqn{m} is the max aggregation order, and \eqn{h} is
#' the forecast horizon for the lowest frequency time series.
#' @inheritParams ctrec
#' @param weights A (\eqn{n_b \times hm}) numeric matrix containing the proportions for each
#' high-frequency bottom time series; \eqn{n_b} is the total number of high-frequency
#' bottom variables, \eqn{m} is the max aggregation order, and \eqn{h} is the forecast horizon
#' for the lowest frequency time series.
#' @inheritParams cstd
#'
#' @inherit ctrec return
#'
#' @examples
#' set.seed(123)
#' # (3 x 1) top base forecasts vector (simulated), forecast horizon = 3
#' topf <- rnorm(3, 10)
#' A <- t(c(1,1)) # Aggregation matrix for Z = X + Y
#'
#' # Same weights for different forecast horizons, agg_order = 4
#' fix_weights <- matrix(runif(4*2), 2, 4)
#' reco <- cttd(base = topf, agg_mat = A, agg_order = 4, weights = fix_weights)
#'
#' # Different weights for different forecast horizons
#' h_weights <- matrix(runif(4*2*3), 2, 3*4)
#' recoh <- cttd(base = topf, agg_mat = A, agg_order = 4, weights = h_weights)
#'
#' @family Reco: top-down
#' @family Framework: cross-temporal
#' @export
#'
cttd <- function(base, agg_mat, agg_order, weights, tew = "sum", normalize = TRUE){
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }else if(NCOL(base) != 1){
    cli_abort("Argument {.arg base} is not a vector.", call = NULL)
  }

  if(missing(agg_mat)){
    cli_abort("Argument {.arg agg_mat} is missing, with no default.", call = NULL)
  }else{
    cstmp <- cstools(agg_mat = agg_mat, sparse = TRUE)
  }

  if(missing(agg_order)){
    cli_abort("Argument {.arg agg_order} is missing, with no default.", call = NULL)
  }else{
    tetmp <- tetools(agg_order = agg_order, sparse = TRUE)
  }

  if(missing(weights)){
    cli_abort("Argument {.arg weights} is missing, with no default.", call = NULL)
  }else if(NROW(weights)!=cstmp$dim[["nb"]]){
    cli_abort("Incorrect {.arg agg_mat} or {.arg weights} dimensions.", call = NULL)
  }

  if(NCOL(weights) != length(base)*tetmp$dim[["m"]]){
    if(NCOL(weights) == tetmp$dim[["m"]]){
      weights <- do.call(cbind, rep(list(weights), length(base)))
    }else{
      cli_abort("Incorrect {.arg weights} length.", call = NULL)
    }
  }

  wsum <- sapply(1:length(base), function(h){
    sum(weights[, rep(1:length(base), each = tetmp$dim[["m"]]) == h])
  })

  if(!normalize){
    if(any(wsum != 1)){
      cli_warn("The {.arg weights} do not add up to 1", call = NULL)
    }
    wsum <- rep(1, length(base))
  }

  hfbbf <- matrix(rep(rep(base/wsum, each = tetmp$dim[["m"]]), cstmp$dim[["nb"]]),
                  cstmp$dim[["nb"]], byrow = TRUE)
  hfbbf <- hfbbf*weights
  reco_mat <- ctbu(hfbbf, agg_mat = cstmp$agg_mat, agg_order = tetmp$set, tew = tew)
  attr(reco_mat, "FoReco") <- list2env(list(framework = "Cross-temporal",
                                            forecast_horizon = length(base),
                                            te_set = tetmp$set,
                                            cs_n = cstmp$dim[["n"]],
                                            rfun = "cttd"))
  return(reco_mat)
}

