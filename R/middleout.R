#' Cross-sectional middle-out reconciliation
#'
#' The middle-out forecast reconciliation (Athanasopoulos et al., 2009) combines
#' top-down ([cstd]) and bottom-up ([csbu]) for genuine hierarchical/grouped
#' time series. Given the base forecasts of variables at an intermediate
#' level \eqn{l}, it performs
#' \itemize{
#'   \item a top-down approach for the levels \eqn{<l};
#'   \item a bottom-up approach for the levels \eqn{>l}.
#' }
#'
#' @param base A (\eqn{h \times n_l}) numeric matrix containing the \eqn{l}-level
#' base forecast; \eqn{n_l} is the number of variables at level \eqn{l}, and
#' \eqn{h} is the forecast horizon.
#' @inheritParams cstd
#' @inheritParams ctmo
#'
#' @inherit csrec return
#'
#' @references
#' Athanasopoulos, G., Ahmed, R. A. and Hyndman, R.J. (2009) Hierarchical forecasts
#' for Australian domestic tourism. \emph{International Journal of Forecasting} 25(1),
#' 146–166. \doi{10.1016/j.ijforecast.2008.07.004}
#'
#' @examples
#' set.seed(123)
#' # Aggregation matrix for Z = X + Y, X = XX + XY and Y = YX + YY
#' A <- matrix(c(1,1,1,1,1,1,0,0,0,0,1,1), 3, byrow = TRUE)
#' # (3 x 2) top base forecasts vector (simulated), forecast horizon = 3
#' baseL2 <- matrix(rnorm(2*3, 5), 3, 2)
#' # Same weights for different forecast horizons
#' fix_weights <- runif(4)
#' reco <- csmo(base = baseL2, agg_mat = A, id_rows = 2:3, weights = fix_weights)
#'
#' # Different weights for different forecast horizons
#' h_weights <- matrix(runif(4*3), 3, 4)
#' recoh <- csmo(base = baseL2, agg_mat = A, id_rows = 2:3, weights = h_weights)
#'
#' @family Reco: middle-out
#' @family Framework: cross-sectional
#' @export
#'
csmo <- function(base, agg_mat, id_rows = 1, weights, normalize = TRUE){
  if(missing(agg_mat)){
    cli_abort("Argument {.arg agg_mat} is missing, with no default.", call = NULL)
  }else{
    tmp <- cstools(agg_mat = agg_mat, sparse = TRUE)
    agg_mat <- tmp$agg_mat
    agg_mat_split <- agg_mat[id_rows, , drop = FALSE]

    # Check all bts in agg_mat_split
    if(sum(agg_mat_split) != tmp$dim[["nb"]] | any(colSums(agg_mat_split) != 1)){
      cli_abort("The matrix selected by {id_rows} rows of {.arg agg_mat},
                must contain only one values equal to 1 in each column.", call = NULL)
    }
  }

  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }else if(is.vector(base)){
    base <- t(base)
  }

  if(NCOL(base) != length(id_rows)){
    cli_abort("Incorrect {.arg id_rows} or {.arg base} dimensions.", call = NULL)
  }

  if(missing(weights)){
    cli_abort("Argument {.arg weights} is missing, with no default.", call = NULL)
  }

  if(NCOL(weights) == 1){
    weights <- t(weights)
  }

  if(NCOL(weights)!=tmp$dim[["nb"]]){
    cli_abort("Incorrect {.arg weights} dimensions.", call = NULL)
  }

  if(NROW(weights)!=NROW(base)){
    weights <- matrix(rep(as.vector(weights), NROW(base)), NROW(base), byrow = TRUE)
  }

  wsum <- t(rbind(apply(weights, 1, function(x)
    as.vector(agg_mat_split%*%x))))

  if(!normalize){
    if(any(wsum != 1)){
      cli_warn("The {.arg weights} do not add up to 1", call = NULL)
    }
    wsum <- matrix(1, NROW(base), NCOL(base))
  }

  bbf <- ((base/wsum) %*% agg_mat_split)*weights
  reco_mat <- csbu(bbf, agg_mat = agg_mat)
  attr(reco_mat, "FoReco") <- list2env(list(framework = "Cross-sectional",
                                            forecast_horizon = NROW(reco_mat),
                                            cs_n = NCOL(reco_mat),
                                            rfun = "csmid"))
  return(reco_mat)
}

#' Temporal middle-out reconciliation
#'
#' The middle-out forecast reconciliation for temporal hierarchies
#' combines top-down ([tetd]) and bottom-up ([tebu]) methods. Given
#' the base forecasts of an intermediate temporal aggregation order \eqn{k}, it performs
#' \itemize{
#'   \item a top-down approach for the aggregation orders \eqn{<k};
#'   \item a bottom-up approach for the aggregation orders \eqn{>k}.
#' }
#'
#' @usage
#' temo(base, agg_order, order = max(agg_order), weights, tew = "sum",
#'      normalize = TRUE)
#'
#' @param base A (\eqn{hk \times 1}) numeric vector containing the temporal
#' aggregated base forecasts of order \eqn{k}; \eqn{k} is an aggregation
#' order (a factor of \eqn{m}, and \eqn{1<k<m}), \eqn{m} is the max aggregation
#' order, and \eqn{h} is the forecast horizon for the lowest frequency time series.
#' @inheritParams tetd
#' @inheritParams ctmo
#'
#' @inherit terec return
#'
#' @examples
#' set.seed(123)
#' # (6 x 1) base forecasts vector (simulated), forecast horizon = 3
#' # and intermediate aggregation order k = 2 (max agg order = 4)
#' basek2 <- rnorm(3*2, 5)
#' # Same weights for different forecast horizons
#' fix_weights <- runif(4)
#' reco <- temo(base = basek2, order = 2, agg_order = 4, weights = fix_weights)
#'
#' # Different weights for different forecast horizons
#' h_weights <- runif(4*3)
#' recoh <- temo(base = basek2, order = 2, agg_order = 4, weights = h_weights)
#'
#' @family Reco: middle-out
#' @family Framework: temporal
#' @export
#'
temo <- function(base, agg_order, order = max(agg_order), weights, tew = "sum", normalize = TRUE){
  if(missing(agg_order)){
    cli_abort("Argument {.arg agg_order} is missing, with no default.", call = NULL)
  }else{
    tmp <- tetools(agg_order = agg_order, sparse = TRUE)
  }

  if(!(order %in% tmp$set)){
    cli_abort("Incorrect {.arg order}
              (different from {tmp$set[-length(tmp$set)]}).", call = NULL)
  }

  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }else if(NCOL(base) != 1){
    cli_abort("Argument {.arg base} is not a vector.", call = NULL)
  }else if(length(base) %% (tmp$dim[["m"]]/order) != 0){
    cli_abort("Incorrect {.arg base} length.", call = NULL)
  }

  if(missing(weights)){
    cli_abort("Argument {.arg weights} is missing, with no default.", call = NULL)
  }else if(NCOL(weights) != 1){
    cli_abort("Argument {.arg weights} is not a vector.", call = NULL)
  }

  h <- length(base)*order/tmp$dim[["m"]]
  if(length(weights) != tmp$dim[["m"]]*h){
    if(length(weights) == tmp$dim[["m"]]){
      weights <- rep(weights, h)
    }else{
      cli_abort("Incorrect {.arg weights} length.", call = NULL)
    }
  }

  wsum <- unname(tapply(weights, rep(1:length(base), each = order), sum))

  if(!normalize){
    if(any(wsum != 1)){
      cli_warn("The {.arg weights} do not add up to 1", call = NULL)
    }
    wsum <- rep(1, length(base))
  }

  hfbf <- rep(base/wsum, each = order)*weights
  reco_vec <- tebu(hfbf, agg_order = tmp$set, tew = tew)
  attr(reco_vec, "FoReco") <- list2env(list(framework = "Temporal",
                                            forecast_horizon = length(base),
                                            te_set = tmp$set,
                                            rfun = "temo"))
  return(reco_vec)
}

#' Cross-temporal middle-out reconciliation
#'
#' The cross-temporal middle-out forecast reconciliation combines top-down
#' ([cttd]) and bottom-up ([ctbu]) methods in the cross-temporal framework for
#' genuine hierarchical/grouped time series. Given the base forecasts of an
#' intermediate cross-sectional level \eqn{l} and aggregation order \eqn{k},
#' it performs
#' \itemize{
#'   \item a top-down approach for the aggregation orders \eqn{\geq k} and
#'   cross-sectional levels \eqn{\geq l};
#'   \item a bottom-up approach, otherwise.
#' }
#'
#' @usage
#' ctmo(base, agg_mat, agg_order, id_rows = 1, order = max(agg_order),
#'      weights, tew = "sum", normalize = TRUE)
#'
#' @param base A (\eqn{n_l \times hk}) numeric matrix containing the \eqn{l}-level
#' base forecasts of temporal aggregation order \eqn{k}; \eqn{n_l} is the number of variables at
#' level \eqn{l}, \eqn{k} is an aggregation order (a factor of \eqn{m}, and \eqn{1<k<m}),
#' \eqn{m} is the max aggregation order, and \eqn{h} is the forecast horizon for the
#' lowest frequency time series.
#' @param id_rows A numeric vector indicating the \eqn{l}-level rows of \code{agg_mat}.
#' @param order The intermediate fixed aggregation order \eqn{k}.
#' @inheritParams cttd
#'
#' @inherit ctrec return
#'
#' @examples
#' set.seed(123)
#' # Aggregation matrix for Z = X + Y, X = XX + XY and Y = YX + YY
#' A <- matrix(c(1,1,1,1,1,1,0,0,0,0,1,1), 3, byrow = TRUE)
#' # (2 x 6) base forecasts matrix (simulated), forecast horizon = 3
#' # and intermediate aggregation order k = 2 (max agg order = 4)
#' baseL2k2 <- rbind(rnorm(3*2, 5), rnorm(3*2, 5))
#'
#' # Same weights for different forecast horizons, agg_order = 4
#' fix_weights <- matrix(runif(4*4), 4, 4)
#' reco <- ctmo(base = baseL2k2, id_rows = 2:3, agg_mat = A,
#'              order = 2, agg_order = 4, weights = fix_weights)
#'
#' # Different weights for different forecast horizons
#' h_weights <- matrix(runif(4*4*3), 4, 3*4)
#' recoh <- ctmo(base = baseL2k2, id_rows = 2:3, agg_mat = A,
#'              order = 2, agg_order = 4, weights = h_weights)
#'
#' @family Reco: middle-out
#' @family Framework: cross-temporal
#' @export
#'
ctmo <- function(base, agg_mat, agg_order, id_rows = 1,
                 order = max(agg_order), weights, tew = "sum", normalize = TRUE){
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }else if(is.vector(base)){
    base <- t(base)
  }

  if(missing(agg_mat)){
    cli_abort("Argument {.arg agg_mat} is missing, with no default.", call = NULL)
  }else{
    cstmp <- cstools(agg_mat = agg_mat, sparse = TRUE)
    strc_mat <- cstmp$strc_mat
    agg_mat_split <- strc_mat[id_rows, , drop = FALSE]
    if(NROW(base) != NROW(agg_mat_split)){
      cli_abort("Incorrect {.arg base} rows dimension.", call = NULL)
    }
  }

  if(missing(agg_order)){
    cli_abort("Argument {.arg agg_order} is missing, with no default.", call = NULL)
  }else{
    tetmp <- tetools(agg_order = agg_order)
    if(!(order %in% tetmp$set)){
      cli_abort("Incorrect {.arg order}
                (different from {tmp$set[-length(tetmp$set)]}).", call = NULL)
    }

    if(NCOL(base) %% (tetmp$dim[["m"]]/order) != 0){
      cli_abort("Incorrect {.arg base} columns dimension.", call = NULL)
    }

    h <- NCOL(base)*order/tetmp$dim[["m"]]
  }

  if(NCOL(weights) != tetmp$dim[["m"]]*h){
    if(NCOL(weights) == tetmp$dim[["m"]]){
      weights <- do.call(cbind, rep(list(weights), h))
    }else{
      cli_abort("Incorrect {.arg weights} length.", call = NULL)
    }
  }

  if(NROW(weights)!=cstmp$dim[["nb"]]){
    cli_abort("Incorrect {.arg weights} row dimensions.", call = NULL)
  }

  idmat <- matrix(c(1:prod(dim(base))), nrow = NROW(base), byrow = TRUE)

  idmat <- do.call(cbind, apply(t(agg_mat_split)%*%idmat, 2, function(x){
    matrix(rep(x, each = order),
           cstmp$dim[["nb"]], byrow = TRUE)
  }, simplify = FALSE))

  base <- do.call(cbind, apply(t(agg_mat_split)%*%base, 2, function(x){
    matrix(rep(x, each = order),
           cstmp$dim[["nb"]], byrow = TRUE)
  }, simplify = FALSE))

  warn_check <- FALSE
  reco_mat <- Reduce("+", lapply(1:max(idmat), function(x){
    id1 <- 1*(idmat==x)
    w <- weights*id1
    if(normalize){
      w <- w/sum(w)
    }else if(sum(w) != 1){
      warn_check <- TRUE
    }
    base*w
  }))

  if(warn_check){
    cli_warn("The {.arg weights} do not add up to 1", call = NULL)
  }

  reco_mat <- ctbu(reco_mat, agg_mat = agg_mat, agg_order = agg_order, tew = tew)

  attr(reco_mat, "FoReco") <- list2env(list(framework = "Cross-temporal",
                                            forecast_horizon = length(base),
                                            te_set = tetmp$set,
                                            cs_n = cstmp$dim[["n"]],
                                            rfun = "ctmo"))

  return(reco_mat)
}
