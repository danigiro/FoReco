#' Cross-sectional bottom-up reconciliation
#'
#' This function computes the cross-sectional bottom-up reconciled forecasts
#' (Dunn et al., 1976) for all series by appropriate summation of the bottom
#' base forecasts \eqn{\widehat{\mathbf{b}}}:
#' \deqn{\widetilde{\mathbf{y}} = \mathbf{S}_{cs}\widehat{\mathbf{b}},}
#' where \eqn{\mathbf{S}_{cs}} is the cross-sectional structural matrix.
#'
#' @param base A (\eqn{h \times n_b}) numeric matrix or multivariate time series
#' (\code{mts} class) containing bottom base forecasts; \eqn{h} is the
#' forecast horizon, and \eqn{n_b} is the total number of bottom variables.
#' @inheritParams csrec
#' @param sntz If \code{TRUE}, the negative base forecasts are
#' set to zero before applying bottom-up.
#'
#' @family Reco: bottom-up
#' @family Framework: cross-sectional
#'
#' @inherit csrec return
#'
#' @references
#' Dunn, D. M., Williams, W. H. and Dechaine, T. L. (1976), Aggregate versus subaggregate
#' models in local area forecasting, \emph{Journal of the American Statistical Association}
#' 71(353), 68–71. \doi{10.1080/01621459.1976.10481478}
#'
#' @examples
#' set.seed(123)
#' # (3 x 2) bottom base forecasts matrix (simulated), Z = X + Y
#' bts <- matrix(rnorm(6, mean = c(10, 10)), 3, byrow = TRUE)
#'
#' # Aggregation matrix for Z = X + Y
#' A <- t(c(1,1))
#' reco <- csbu(base = bts, agg_mat = A)
#'
#' # Non negative reconciliation
#' bts[2,2] <- -bts[2,2] # Making negative one of the base forecasts for variable Y
#' nnreco <- csbu(base = bts, agg_mat = A, sntz = TRUE)
#'
#' @export
#'
csbu <- function(base, agg_mat, sntz = FALSE){
  # base forecasts condition
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }

  if(NCOL(base) == 1){
    base <- t(base)
  }

  if(missing(agg_mat)){
    cli_abort("Argument {.arg agg_mat} is missing, with no default.", call = NULL)
  }else{
    tmp <- cstools(agg_mat = agg_mat, sparse = TRUE)
  }

  csdim <- tmp$dim
  agg_mat <- tmp$agg_mat
  strc_mat <- tmp$strc_mat
  nn_cons_var <- c(rep(0, csdim[2]), rep(1, csdim[3]))

  if(NCOL(base) != csdim[["nb"]]){
    cli_abort("Incorrect {.arg agg_mat} or {.arg base} dimensions.", call = NULL)
  }

  if(sntz){
    base[base<0] <- 0
  }

  reco_mat <- as.matrix(base%*%t(strc_mat))
  rownames(reco_mat) <- paste0("h-", 1:NROW(reco_mat))
  if(!is.null(colnames(agg_mat)) & !is.null(rownames(agg_mat))){
    colnames(reco_mat) <- unlist(dimnames(agg_mat))
  }else{
    colnames(reco_mat) <- paste0("s-", 1:NCOL(reco_mat))
  }
  attr(reco_mat, "FoReco") <- list2env(list(framework = "Cross-sectional",
                                            forecast_horizon = NROW(reco_mat),
                                            cs_n = csdim[["n"]],
                                            rfun = "csbu"))
  return(reco_mat)
}

#' Temporal bottom-up reconciliation
#'
#' Temporal bottom-up reconciled forecasts at any temporal aggregation level are computed by
#' appropriate aggregation of the high-frequency base forecasts, \eqn{\widehat{\mathbf{x}}^{[1]}}:
#' \deqn{\widetilde{\mathbf{x}} = \mathbf{S}_{te}\widehat{\mathbf{x}}^{[1]},}
#' where \eqn{\mathbf{S}_{te}} is the temporal structural matrix.
#'
#' @param base A (\eqn{hm \times 1}) numeric vector containing the high-frequency base forecasts;
#' \eqn{m} is the max. temporal aggregation order, and \eqn{h} is the forecast horizon for the
#' lowest frequency time series.
#' @inheritParams terec
#' @inheritParams csbu
#'
#' @inherit terec return
#'
#' @examples
#' set.seed(123)
#' # (4 x 1) high frequency base forecasts vector (simulated),
#' # agg_order = 4 (annual-quarterly)
#' hfts <- rnorm(4, 5)
#'
#' reco <- tebu(base = hfts, agg_order = 4)
#'
#' # Non negative reconciliation
#' hfts[4] <- -hfts[4] # Making negative one of the quarterly base forecasts
#' nnreco <- tebu(base = hfts, agg_order = 4, sntz = TRUE)
#'
#' @family Reco: bottom-up
#' @family Framework: temporal
#' @export
#'
tebu <- function(base, agg_order, tew = "sum", sntz = FALSE){

  if(missing(agg_order)){
    cli_abort("Argument {.arg agg_order} is missing, with no default.", call = NULL)
  }

  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }else if(NCOL(base) != 1){
    cli_abort("{.arg base} is not a vector.", call = NULL)
  }else if(length(base) %% max(agg_order) != 0){
    cli_abort("Incorrect {.arg base} length.", call = NULL)
  }

  h <- length(base)/max(agg_order)

  tmp <- tetools(agg_order = agg_order, tew = tew)
  kset <- tmp$set
  m <- tmp$dim[["m"]]
  strc_mat <- tmp$strc_mat

  base <- matrix(base, nrow = h, ncol = m, byrow = TRUE)
  if(sntz){
    base[base<0] <- 0
  }

  reco_vec <- base%*%t(strc_mat)
  reco_vec <- hmat2vec(reco_vec, h = h, kset = kset)
  attr(reco_vec, "FoReco") <- list2env(list(framework = "Temporal",
                                            forecast_horizon = h,
                                            te_set = tmp$set,
                                            rfun = "tebu"))
  return(reco_vec)
}

#' Cross-temporal bottom-up reconciliation
#'
#' Cross-temporal bottom-up reconciled forecasts for all series at any temporal
#' aggregation level are computed by appropriate summation of the high-frequency
#' bottom base forecasts \eqn{\widehat{\mathbf{B}^{[1]}}}:
#' \deqn{\widetilde{\mathbf{X}} = \mathbf{S}_{cs}\widehat{\mathbf{B}^{[1]}}\mathbf{S}'_{te},}
#' where \eqn{\mathbf{S}_{cs}} and \eqn{\mathbf{S}_{te}} are the cross-sectional and
#' temporal structural matrices, respectively.
#'
#' @param base A (\eqn{n_b \times hm}) numeric matrix containing high-frequency bottom base
#' forecasts; \eqn{n_b} is the total number of high-frequency bottom variables, \eqn{m} is
#' the max aggregation order, and \eqn{h} is the forecast horizon for the lowest frequency
#' time series.
#' @inheritParams ctrec
#' @inheritParams csbu
#'
#' @inherit ctrec return
#'
#' @examples
#' set.seed(123)
#' # Aggregation matrix for Z = X + Y
#' A <- t(c(1,1))
#' # (2 x 4) high frequency bottom base forecasts matrix (simulated),
#' # agg_order = 4 (annual-quarterly)
#' hfbts <- matrix(rnorm(4*2, 2.5), 2, 4)
#'
#' reco <- ctbu(base = hfbts, agg_mat = A, agg_order = 4)
#'
#' # Non negative reconciliation
#' hfbts[1,4] <- -hfbts[1,4] # Making negative one of the quarterly base forecasts for variable X
#' nnreco <- ctbu(base = hfbts, agg_mat = A, agg_order = 4, sntz = TRUE)
#'
#' @family Reco: bottom-up
#' @family Framework: cross-temporal
#' @export
#'
ctbu <- function(base, agg_mat, agg_order, tew = "sum", sntz = FALSE){
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }

  if(missing(agg_mat)){
    cli_abort("Argument {.arg agg_mat} is missing, with no default.", call = NULL)
  }else{
    cstmp <- cstools(agg_mat = agg_mat, sparse = TRUE)
    cs_strc_mat <- cstmp$strc_mat
    if(NROW(base) != cstmp$dim[["nb"]]){
      cli_abort("Incorrect {.arg base} rows dimension.", call = NULL)
    }
  }

  if(missing(agg_order)){
    cli_abort("Argument {.arg agg_order} is missing, with no default.", call = NULL)
  }else{
    if(NCOL(base) %% max(agg_order) != 0){
      cli_abort("Incorrect {.arg base} columns dimension.", call = NULL)
    }
    h <- NCOL(base)/max(agg_order)
    tetmp <- tetools(agg_order = agg_order, tew = tew, fh = h)
    te_strc_mat <- tetmp$strc_mat
  }

  if(sntz){
    base[base<0] <- 0
  }

  reco_mat <- as.matrix(cs_strc_mat %*% base %*% t(te_strc_mat))

  rownames(reco_mat) <- namesCS(n = NROW(reco_mat), names_list = dimnames(agg_mat))
  colnames(reco_mat) <- namesTE(kset = tetmp$set, h = h)
  attr(reco_mat, "FoReco") <- list2env(list(framework = "Cross-temporal",
                                            forecast_horizon = h,
                                            te_set = tetmp$set,
                                            cs_n = cstmp$dim[["n"]],
                                            rfun = "ctbu"))
  return(reco_mat)
}
