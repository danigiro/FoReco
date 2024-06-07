#' Level conditional coherent reconciliation for genuine hierarchical/grouped time series
#'
#' @description
#' This function implements the cross-sectional forecast reconciliation procedure that
#' extends the original proposal by Hollyman et al. (2021). Level conditional coherent
#' reconciled forecasts are conditional on (i.e., constrained by) the base forecasts
#' of a specific upper level in the hierarchy (exogenous constraints). It also allows
#' handling the linear constraints linking the variables endogenously (Di Fonzo and
#' Girolimetto, 2022). The function can calculate Combined Conditional Coherent (CCC)
#' forecasts as simple averages of Level-Conditional Coherent (LCC) and bottom-up
#' reconciled forecasts, with either endogenous or exogenous constraints.
#'
#' @usage
#' cslcc(base, agg_mat, nodes = "auto", comb = "ols", res = NULL, CCC = TRUE,
#'       const = "exogenous", bts = NULL, approach = "proj", nn = NULL,
#'       settings = NULL, ...)
#'
#' @inheritParams csrec
#' @param nodes A (\eqn{L \times 1}) numeric vector indicating the number of variables
#' in each of the upper \eqn{L} levels of the hierarchy. The \emph{default}
#' value is the string "\code{auto}" which calculates the number of variables in each level.
#' @param const A string specifying the reconciliation constraints:
#' \itemize{
#'   \item "\code{exogenous}" (\emph{default}): Fixes the top level of each sub-hierarchy.
#'   \item "\code{endogenous}": Coherently revises both the top and bottom levels.
#' }
#' @param bts A (\eqn{h \times n_b}) numeric matrix or multivariate time series (\code{mts} class)
#' containing bottom base forecasts defined by the user (e.g., seasonal averages, as in Hollyman et al., 2021).
#' This parameter can be omitted if only base forecasts are used
#' (see Di Fonzo and Girolimetto, 2024).
#' @param CCC A logical value indicating whether the Combined Conditional Coherent reconciled
#' forecasts reconciliation should include bottom-up forecasts (\code{TRUE}, \emph{default}), or not.
#' @inheritDotParams cscov mse shrink_fun
#'
#' @inherit csrec return
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
#' Di Fonzo, T. and Girolimetto, D. (2024), Forecast combination-based forecast reconciliation:
#' Insights and extensions, \emph{International Journal of Forecasting}, 40(2), 490–514.
#' \doi{10.1016/j.ijforecast.2022.07.001}
#'
#' Di Fonzo, T. and Girolimetto, D. (2023b) Spatio-temporal reconciliation of solar forecasts.
#' \emph{Solar Energy} 251, 13–29. \doi{10.1016/j.solener.2023.01.003}
#'
#' Hyndman, R.J., Ahmed, R.A., Athanasopoulos, G. and Shang, H.L. (2011),
#' Optimal combination forecasts for hierarchical time series,
#' \emph{Computational Statistics & Data Analysis}, 55, 9, 2579-2589.
#' \doi{10.1016/j.csda.2011.03.006}
#'
#' Hollyman, R., Petropoulos, F. and Tipping, M.E. (2021), Understanding forecast reconciliation.
#' \emph{European Journal of Operational Research}, 294, 149–160. \doi{10.1016/j.ejor.2021.01.017}
#'
#' Stellato, B., Banjac, G., Goulart, P., Bemporad, A. and Boyd, S. (2020), OSQP:
#' An Operator Splitting solver for Quadratic Programs,
#' \emph{Mathematical Programming Computation}, 12, 4, 637-672.
#' \doi{10.1007/s12532-020-00179-2}
#'
#' @examples
#' set.seed(123)
#' # Aggregation matrix for Z = X + Y, X = XX + XY and Y = YX + YY
#' A <- matrix(c(1,1,1,1,1,1,0,0,0,0,1,1), 3, byrow = TRUE)
#' # (2 x 7) base forecasts matrix (simulated)
#' base <- matrix(rnorm(7*2, mean = c(40, 20, 20, 10, 10, 10, 10)), 2, byrow = TRUE)
#' # (10 x 7) in-sample residuals matrix (simulated)
#' res <- matrix(rnorm(n = 7*10), ncol = 7)
#' # (2 x 7) Naive bottom base forecasts matrix: all forecasts are set equal to 10
#' naive <- matrix(10, 2, 4)
#'
#' ## EXOGENOUS CONSTRAINTS (Hollyman et al., 2021)
#' # Level Conditional Coherent (LCC) reconciled forecasts
#' exo_LC <- cslcc(base = base, agg_mat = A, comb = "wls", bts = naive,
#'                 res = res, nodes = "auto", CCC = FALSE)
#'
#' # Combined Conditional Coherent (CCC) reconciled forecasts
#' exo_CCC <- cslcc(base = base, agg_mat = A, comb = "wls", bts = naive,
#'                  res = res, nodes = "auto", CCC = TRUE)
#'
#' # Results detailed by level:
#' # L-1: Level 1 immutable reconciled forecasts for the whole hierarchy
#' # L-2: Middle-Out reconciled forecasts
#' # L-3: Bottom-Up reconciled forecasts
#' info_exo <- recoinfo(exo_CCC, verbose = FALSE)
#' info_exo$lcc
#'
#' ## ENDOGENOUS CONSTRAINTS (Di Fonzo and Girolimetto, 2024)
#' # Level Conditional Coherent (LCC) reconciled forecasts
#' endo_LC <- cslcc(base = base, agg_mat = A, comb = "wls",
#'                  res = res, nodes = "auto", CCC = FALSE,
#'                  const = "endogenous")
#'
#' # Combined Conditional Coherent (CCC) reconciled forecasts
#' endo_CCC <- cslcc(base = base, agg_mat = A, comb = "wls",
#'                   res = res, nodes = "auto", CCC = TRUE,
#'                   const = "endogenous")
#'
#' # Results detailed by level:
#' # L-1: Level 1 reconciled forecasts for L1 + L3 (bottom level)
#' # L-2: Level 2 reconciled forecasts for L2 + L3 (bottom level)
#' # L-3: Bottom-Up reconciled forecasts
#' info_endo <- recoinfo(endo_CCC, verbose = FALSE)
#' info_endo$lcc
#'
#' @family Reco: level conditional coherent
#' @family Framework: cross-sectional
#' @export
cslcc <- function(base, agg_mat, nodes = "auto", comb = "ols", res = NULL,
                  CCC = TRUE, const = "exogenous", bts = NULL,
                  approach = "proj", nn = NULL, settings = NULL, ...){

  const <- match.arg(const, c("exogenous", "endogenous"))

  # Check if 'agg_mat' is specified
  if(missing(agg_mat)){
    utd <- TRUE
    agg_mat <- Matrix(1, ncol = NCOL(base)-1)
    nodes <- 1
  } else {
    utd <- FALSE
  }

  # Check if 'base' is provided and its dimensions match with the data
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  } else if(NCOL(base) == 1){
    base <- t(base)
  }

  if(NCOL(base) != sum(dim(agg_mat))){
    cli_abort("Incorrect {.arg base} columns dimension.", call = NULL)
  }

  # Balanced hierarchy
  balh <- balance_hierarchy(agg_mat = agg_mat, nodes = nodes)
  names_base <- colnames(base)
  base <- base[, balh$id, drop = FALSE]
  nodes <- balh$nodes
  agg_mat <- balh$bam
  ubam <- balh$agg_mat

  # Cross-sectional tools
  n <- sum(dim(agg_mat))
  na <- NROW(agg_mat)
  nb <- NCOL(agg_mat)
  strc_mat <- rbind(agg_mat, .sparseDiagonal(nb))

  # Check if 'res' is provided and its dimensions match with the data - TODO
  if(!is.null(res)){
    res <- res[, balh$id]
  }

  lev_num <- 1:(length(nodes)+1)
  lev_bts <- max(lev_num)
  lev_id <- rep(lev_num , c(nodes, nb))

  # check bts
  if(!is.null(bts)){
    if(NCOL(bts) != nb | NROW(bts) != NROW(base)){
      cli_abort("Incorrect {.arg bts} dimensions ({NROW(base)} x {nb}).", call = NULL)
    }
  }else{
    bts <- base[, sort(which(lev_id == lev_bts)), drop = FALSE]
  }

  # Compute covariance
  cov_mat <- cscov(comb = comb, n = n, agg_mat = agg_mat, res = res, ...)
  if(NROW(cov_mat) != n | NCOL(cov_mat) != n){
    cli_abort(c("Incorrect covariance dimensions.",
                "i"="Check {.arg res} dimensions."), call = NULL)
  }

  lccmat <- lapply(lev_num, function(x){
    idx <- sort(which(lev_id == x))

    if(x != lev_bts){
      if(const == "exogenous"){
        imx <- 1:length(idx)
      }else{
        imx <- NULL
      }

      basex <- cbind(base[, idx, drop = FALSE], bts)
      idxb <- sort(which(lev_id %in% c(x, lev_bts)))
      strc_matx <- strc_mat[idxb, , drop = FALSE]
      cons_matx <- cbind(.sparseDiagonal(length(idx)), -agg_mat[idx, , drop = FALSE])
      cov_matx <- cov_mat[idxb, idxb, drop = FALSE]

      rmat <- reco(base = basex,
                   cov_mat = cov_matx,
                   strc_mat = strc_matx,
                   cons_mat = cons_matx,
                   id_nn = as.numeric(lev_id == lev_bts)[idxb],
                   approach = approach,
                   nn = nn,
                   immutable = imx,
                   settings = settings)

      rbts <- rmat[, -c(1:nodes[x]), drop = FALSE]
      rmat <- csbu(rbts, agg_mat = ubam)
    }else{
      rmat <- csbu(base[, idx, drop = FALSE], agg_mat = ubam, sntz = !is.null(nn))
    }
    colnames(rmat) <- namesCS(n = NCOL(rmat), names_vec = names_base,
                              names_list = dimnames(ubam))
    rmat
  })

  if(utd){ # Unbiased top-down forecasts
    out <- lccmat[[1]]
    attr(out, "FoReco") <- list2env(list(info = attr(lccmat[[1]], "info"),
                                         framework = "Cross-sectional",
                                         forecast_horizon = NROW(out),
                                         comb = comb,
                                         cs_n = n,
                                         rfun = "cslcc"))
    return(out)
  }else{ # Mean
    if(CCC){ # Mean LCC + BU
      out <- Reduce("+", lccmat)/length(lccmat)
    }else{   # Mean LCC
      out <- Reduce("+", lccmat[-length(lccmat)])/(length(lccmat)-1)
    }
    names(lccmat) <- paste0("L-", lev_num)
    attr(out, "FoReco") <- list2env(list(lcc = lccmat,
                                         framework = "Cross-sectional",
                                         forecast_horizon = NROW(out),
                                         comb = comb,
                                         cs_n = n,
                                         rfun = "cslcc"))
    return(out)
  }
}

#' Level conditional coherent reconciliation for temporal hierarchies
#'
#' @description
#' This function implements a forecast reconciliation procedure inspired by the original proposal
#' by Hollyman et al. (2021) for temporal hierarchies. Level conditional coherent
#' reconciled forecasts are conditional on (i.e., constrained by) the base forecasts
#' of a specific upper level in the hierarchy (exogenous constraints). It also allows
#' handling the linear constraints linking the variables endogenously (Di Fonzo and
#' Girolimetto, 2022). The function can calculate Combined Conditional Coherent (CCC)
#' forecasts as simple averages of Level-Conditional Coherent (LCC) and bottom-up
#' reconciled forecasts, with either endogenous or exogenous constraints.
#'
#' @usage
#' telcc(base, agg_order, comb = "ols", res = NULL, CCC = TRUE,
#'       const = "exogenous", hfts = NULL, tew = "sum",
#'       approach = "proj", nn = NULL, settings = NULL, ...)
#'
#' @inheritParams terec
#' @inheritParams cslcc
#' @param hfts A (\eqn{mh \times 1}) numeric vector containing high frequency base forecasts defined
#' by the user. This parameter can be omitted if only base forecasts in
#' \code{base} are used (see Di Fonzo and Girolimetto, 2024).
#' @inheritDotParams tecov mse shrink_fun
#'
#' @inherit terec return
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
#' Di Fonzo, T. and Girolimetto, D. (2024), Forecast combination-based forecast reconciliation:
#' Insights and extensions, \emph{International Journal of Forecasting}, 40(2), 490–514.
#' \doi{10.1016/j.ijforecast.2022.07.001}
#'
#' Di Fonzo, T. and Girolimetto, D. (2023b) Spatio-temporal reconciliation of solar forecasts.
#' \emph{Solar Energy} 251, 13–29. \doi{10.1016/j.solener.2023.01.003}
#'
#' Hyndman, R.J., Ahmed, R.A., Athanasopoulos, G. and Shang, H.L. (2011),
#' Optimal combination forecasts for hierarchical time series,
#' \emph{Computational Statistics & Data Analysis}, 55, 9, 2579-2589.
#' \doi{10.1016/j.csda.2011.03.006}
#'
#' Hollyman, R., Petropoulos, F. and Tipping, M.E. (2021), Understanding forecast reconciliation.
#' \emph{European Journal of Operational Research}, 294, 149–160. \doi{10.1016/j.ejor.2021.01.017}
#'
#' Stellato, B., Banjac, G., Goulart, P., Bemporad, A. and Boyd, S. (2020), OSQP:
#' An Operator Splitting solver for Quadratic Programs,
#' \emph{Mathematical Programming Computation}, 12, 4, 637-672.
#' \doi{10.1007/s12532-020-00179-2}
#'
#' @examples
#' set.seed(123)
#' # (7 x 1) base forecasts vector (simulated), agg_order = 4
#' base <- rnorm(7, rep(c(20, 10, 5), c(1, 2, 4)))
#' # (70 x 1) in-sample residuals vector (simulated)
#' res <- rnorm(70)
#' # (4 x 1) Naive high frequency base forecasts vector: all forecasts are set equal to 2.5
#' naive <- rep(2.5, 4)
#'
#' ## EXOGENOUS CONSTRAINTS
#' # Level Conditional Coherent (LCC) reconciled forecasts
#' exo_LC <- telcc(base = base, agg_order = 4, comb = "wlsh", hfts = naive,
#'                 res = res, nodes = "auto", CCC = FALSE)
#'
#' # Combined Conditional Coherent (CCC) reconciled forecasts
#' exo_CCC <- telcc(base = base, agg_order = 4, comb = "wlsh", hfts = naive,
#'                  res = res, nodes = "auto", CCC = TRUE)
#'
#' # Results detailed by level:
#' info_exo <- recoinfo(exo_CCC, verbose = FALSE)
#' # info_exo$lcc
#'
#' ## ENDOGENOUS CONSTRAINTS
#' # Level Conditional Coherent (LCC) reconciled forecasts
#' endo_LC <- telcc(base = base, agg_order = 4, comb = "wlsh", res = res,
#'                  nodes = "auto", CCC = FALSE, const = "endogenous")
#'
#' # Combined Conditional Coherent (CCC) reconciled forecasts
#' endo_CCC <- telcc(base = base, agg_order = 4, comb = "wlsh", res = res,
#'                   nodes = "auto", CCC = TRUE, const = "endogenous")
#'
#' # Results detailed by level:
#' info_endo <- recoinfo(endo_CCC, verbose = FALSE)
#' # info_endo$lcc
#'
#' @family Reco: level conditional coherent
#' @family Framework: temporal
#' @export
telcc <- function(base, agg_order, comb = "ols", res = NULL, CCC = TRUE,
                  const = "exogenous", hfts = NULL, tew = "sum",
                  approach = "proj", nn = NULL, settings = NULL, ...){

  # Check if 'agg_order' is provided
  if(missing(agg_order)){
    cli_abort("Argument {.arg agg_order} is missing, with no default.", call = NULL)
  }

  tmp <- tetools(agg_order = agg_order, tew = tew)
  kset <- tmp$set
  kt <- tmp$dim[["kt"]]
  strc_mat <- tmp$strc_mat
  agg_mat <- tmp$agg_mat
  nodes <- tmp$dim[["m"]]/kset[-length(kset)]

  # Check if 'base' is provided and its dimensions match with the data
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  } else if(NCOL(base) != 1){
    cli_abort("{.arg base} is not a vector.", call = NULL)
  }

  # Calculate 'h' and 'base_hmat'
  if(length(base) %% kt != 0){
    cli_abort("Incorrect {.arg base} length.", call = NULL)
  } else {
    h <- length(base) / kt
    base <- vec2hmat(vec = base, h = h, kset = kset)
  }

  lev_num <- 1:(length(nodes)+1)
  lev_hfts <- max(lev_num)
  lev_id <- rep(lev_num , c(nodes, tmp$dim[["m"]]))

  # check bts
  if(!is.null(hfts)){
    if(NCOL(hfts) != 1 | NROW(hfts) != h*tmp$dim[["m"]]){
      cli_abort('Incorrect {.arg hfts} lenght ({h*tmp$dim[["m"]]}).', call = NULL)
    }
    hfts <- matrix(hfts, nrow = h, byrow = TRUE)
  }else{
    hfts <- base[, sort(which(lev_id == lev_hfts)), drop = FALSE]
  }

  # Compute covariance
  cov_mat <- tecov(comb = comb, res = res, agg_order = agg_order, tew = tew, ...)
  if(NROW(cov_mat) != kt | NCOL(cov_mat) != kt){
    cli_abort(c("Incorrect covariance dimensions.",
                "i"="Check {.arg res} length."), call = NULL)
  }

  lccmat <- lapply(lev_num, function(x){
    idx <- sort(which(lev_id == x))

    if(x != lev_hfts){
      if(const == "exogenous"){
        imx <- idx
      }else{
        imx <- NULL
      }

      basex <- cbind(base[, idx, drop = FALSE], hfts)
      idxb <- sort(which(lev_id %in% c(x, lev_hfts)))
      strc_matx <- strc_mat[idxb, , drop = FALSE]
      cons_matx <- cbind(.sparseDiagonal(length(idx)), -agg_mat[idx, , drop = FALSE])
      cov_matx <- cov_mat[idxb, idxb, drop = FALSE]

      rmat <- reco(base = basex,
                   cov_mat = cov_matx,
                   strc_mat = strc_matx,
                   cons_mat = cons_matx,
                   id_nn = as.numeric(lev_id == lev_hfts)[idxb],
                   approach = approach,
                   nn = nn,
                   immutable = 1:length(idx),
                   settings = settings)

      rhfts <- rmat[, -c(1:nodes[x]), drop = FALSE]

      tebu(as.vector(t(rhfts)), agg_order = agg_order, tew = tew)
    }else{
      tebu(as.vector(t(base[, idx, drop = FALSE])), agg_order = agg_order,
           tew = tew, sntz = !is.null(nn))
    }
  })


  if(CCC){ # Mean LCC + BU
    out <- Reduce("+", lccmat)/length(lccmat)
  }else{   # Mean LCC
    out <- Reduce("+", lccmat[-length(lccmat)])/(length(lccmat)-1)
  }
  names(lccmat) <- paste0("k-", kset)
  attr(out, "FoReco") <- list2env(list(lcc = lccmat,
                                       framework = "Temporal",
                                       forecast_horizon = h,
                                       comb = comb,
                                       te_set = tmp$set,
                                       rfun = "telcc"))
  return(out)
}


#' Level conditional coherent reconciliation for cross-temporal hierarchies
#'
#' @description
#' This function implements a forecast reconciliation procedure inspired by the original
#' proposal by Hollyman et al. (2021) for cross-temporal hierarchies. Level conditional
#' coherent reconciled forecasts are conditional on (i.e., constrained by) the base
#' forecasts of a specific upper level in the hierarchy (exogenous constraints). It also
#' allows handling the linear constraints linking the variables endogenously (Di Fonzo
#' and Girolimetto, 2022). The function can calculate Combined Conditional Coherent (CCC)
#' forecasts as simple averages of Level-Conditional Coherent (LCC) and bottom-up
#' reconciled forecasts, with either endogenous or exogenous constraints.
#'
#' @usage
#' ctlcc(base, agg_mat, nodes = "auto", agg_order, comb = "ols", res = NULL,
#'       CCC = TRUE, const = "exogenous", hfbts = NULL, tew = "sum",
#'       approach = "proj", nn = NULL, settings = NULL, ...)
#'
#' @inheritParams ctrec
#' @inheritParams cslcc
#' @param hfbts A (\eqn{n \times mh}) numeric matrix containing high frequency bottom base
#' forecasts defined by the user. This parameter can be omitted
#' if only base forecasts are used (see Di Fonzo and Girolimetto, 2024).
#' @inheritDotParams ctcov mse shrink_fun
#'
#' @inherit ctrec return
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
#' Di Fonzo, T. and Girolimetto, D. (2024), Forecast combination-based forecast reconciliation:
#' Insights and extensions, \emph{International Journal of Forecasting}, 40(2), 490–514.
#' \doi{10.1016/j.ijforecast.2022.07.001}
#'
#' Di Fonzo, T. and Girolimetto, D. (2023b) Spatio-temporal reconciliation of solar forecasts.
#' \emph{Solar Energy} 251, 13–29. \doi{10.1016/j.solener.2023.01.003}
#'
#' Hyndman, R.J., Ahmed, R.A., Athanasopoulos, G. and Shang, H.L. (2011),
#' Optimal combination forecasts for hierarchical time series,
#' \emph{Computational Statistics & Data Analysis}, 55, 9, 2579-2589.
#' \doi{10.1016/j.csda.2011.03.006}
#'
#' Hollyman, R., Petropoulos, F. and Tipping, M.E. (2021), Understanding forecast reconciliation.
#' \emph{European Journal of Operational Research}, 294, 149–160. \doi{10.1016/j.ejor.2021.01.017}
#'
#' Stellato, B., Banjac, G., Goulart, P., Bemporad, A. and Boyd, S. (2020), OSQP:
#' An Operator Splitting solver for Quadratic Programs,
#' \emph{Mathematical Programming Computation}, 12, 4, 637-672.
#' \doi{10.1007/s12532-020-00179-2}
#'
#' @examples
#' set.seed(123)
#' # Aggregation matrix for Z = X + Y, X = XX + XY and Y = YX + YY
#' A <- matrix(c(1,1,1,1,1,1,0,0,0,0,1,1), 3, byrow = TRUE)
#' # (7 x 7) base forecasts matrix (simulated), agg_order = 4
#' base <- rbind(rnorm(7, rep(c(40, 20, 10), c(1, 2, 4))),
#'               rnorm(7, rep(c(20, 10, 5), c(1, 2, 4))),
#'               rnorm(7, rep(c(20, 10, 5), c(1, 2, 4))),
#'               rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))),
#'               rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))),
#'               rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))),
#'               rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))))
#' # (7 x 70) in-sample residuals matrix (simulated)
#' res <- matrix(rnorm(70*7), nrow = 7)
#' # (4 x 4) Naive high frequency bottom base forecasts vector:
#' # all forecasts are set equal to 2.5
#' naive <- matrix(2.5, 4, 4)
#'
#' ## EXOGENOUS CONSTRAINTS (Hollyman et al., 2021)
#' # Level Conditional Coherent (LCC) reconciled forecasts
#' exo_LC <- ctlcc(base = base, agg_mat = A, agg_order = 4, comb = "wlsh", nn = "osqp",
#'                 hfbts = naive, res = res, nodes = "auto", CCC = FALSE)
#'
#' # Combined Conditional Coherent (CCC) reconciled forecasts
#' exo_CCC <- ctlcc(base = base, agg_mat = A, agg_order = 4, comb = "wlsh",
#'                 hfbts = naive, res = res, nodes = "auto", CCC = TRUE)
#'
#' # Results detailed by level:
#' info_exo <- recoinfo(exo_CCC, verbose = FALSE)
#' # info_exo$lcc
#'
#' ## ENDOGENOUS CONSTRAINTS (Di Fonzo and Girolimetto, 2024)
#' # Level Conditional Coherent (LCC) reconciled forecasts
#' endo_LC <- ctlcc(base = base, agg_mat = A, agg_order = 4, comb = "wlsh",
#'                  res = res, nodes = "auto", CCC = FALSE,
#'                  const = "endogenous")
#'
#' # Combined Conditional Coherent (CCC) reconciled forecasts
#' endo_CCC <- ctlcc(base = base, agg_mat = A, agg_order = 4, comb = "wlsh",
#'                   res = res, nodes = "auto", CCC = TRUE,
#'                   const = "endogenous")
#'
#' # Results detailed by level:
#' info_endo <- recoinfo(endo_CCC, verbose = FALSE)
#' # info_endo$lcc
#'
#' @family Reco: level conditional coherent
#' @family Framework: cross-temporal
#' @export
ctlcc <- function(base, agg_mat, nodes = "auto", agg_order, comb = "ols", res = NULL,
                  CCC = TRUE, const = "exogenous", hfbts = NULL, tew = "sum",
                  approach = "proj", nn = NULL, settings = NULL, ...){

  # Check if 'agg_order' is provided
  if(missing(agg_order)){
    cli_abort("Argument {.arg agg_order} is missing, with no default.", call = NULL)
  }

  # Check if 'base' is provided and its dimensions match with the data
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }

  if(NROW(base) != sum(dim(agg_mat))){
    cli_abort("Incorrect {.arg base} rows dimension.", call = NULL)
  }

  # Check if 'agg_mat' is provided and balance the hierarchy
  if(missing(agg_mat)){
    cli_abort("Argument {.arg agg_mat} is missing, with no default.", call = NULL)
  }else{
    balh <- balance_hierarchy(agg_mat = agg_mat, nodes = nodes)
    names_base <- rownames(base)
    base <- base[balh$id,]
    nodes <- balh$nodes
    agg_mat <- balh$bam
    ubam <- balh$agg_mat
    n <- sum(dim(agg_mat))
    na <- NROW(agg_mat)
    nb <- NCOL(agg_mat)
    tmp <- cttools(agg_mat = agg_mat, agg_order = agg_order, tew = tew)
    m <- tmp$dim[["m"]]
    kset <- tmp$set
    strc_mat <- tmp$strc_mat
  }

  # Calculate 'h' and 'base_hmat'
  if(NCOL(base) %% tmp$dim[["kt"]] != 0){
    cli_abort("Incorrect {.arg base} columns dimension.", call = NULL)
  }else{
    h <- NCOL(base) / tmp$dim[["kt"]]
    base_hmat <- mat2hmat(base, h = h, kset = tmp$set, n = tmp$dim[["n"]])
  }

  # Balance the residuals
  if(!is.null(res)){
    res <- res[balh$id, ]
  }

  # Compute covariance
  cov_mat <- ctcov(comb = comb, res = res, agg_order = agg_order, agg_mat = agg_mat,
                    n = tmp$dim[["n"]], tew = tew, ...)
  if(NROW(cov_mat) != prod(tmp$dim[c("kt", "n")]) | NCOL(cov_mat) != prod(tmp$dim[c("kt", "n")])){
    cli_abort(c("Incorrect covariance dimensions.",
                "i"="Check {.arg res} dimensions."), call = NULL)
  }

  lev_num <- 1:((length(nodes)+1)*length(kset))
  MLk <- matrix(1:length(lev_num), ncol = length(kset),  byrow = TRUE)
  lev_id <- as.vector(sapply(rep(split(MLk, 1:NROW(MLk)), c(nodes, nb)),
                             function(x) rep(x, m/kset)))
  lev_bts <- max(lev_num)

  # check bts
  if(!is.null(hfbts)){
    if(NCOL(hfbts) != h*m | NROW(hfbts) != nb){
      cli_abort("Incorrect {.arg hfbts} dimensions ({nb} x {h*m}).", call = NULL)
    }
    hfbts <- matrix(as.vector(hfbts), h, m*nb)
  }else{
    hfbts <- base_hmat[, sort(which(lev_id == lev_bts)), drop = FALSE]
  }

  lccmat <- lapply(lev_num, function(x){
    idx <- sort(which(lev_id == x))

    if(x != lev_bts){
      if(const == "exogenous"){
        imx <- idx
      }else{
        imx <- NULL
      }

      basex <- cbind(base_hmat[, idx, drop = FALSE], hfbts)
      idxb <- sort(which(lev_id %in% c(x, lev_bts)))
      strc_matx <- strc_mat[idxb, , drop = FALSE]
      cons_matx <- cbind(.sparseDiagonal(length(idx)), -strc_mat[idx, , drop = FALSE])
      cov_matx <- cov_mat[idxb, idxb, drop = FALSE]

      rmat <- reco(base = basex,
                   cov_mat = cov_matx,
                   strc_mat = strc_matx,
                   cons_mat = cons_matx,
                   id_nn = as.numeric(lev_id == lev_bts)[idxb],
                   approach = approach,
                   nn = nn,
                   immutable = 1:length(idx),
                   settings = settings)

      rhfbts <- rmat[, -c(1:length(idx)), drop = FALSE]

      rmat <- as.matrix(rhfbts%*%t(strc_mat))
      remid <- !duplicated(balh$id, fromLast = TRUE)
      rmat <- hmat2mat(rmat, h = h, kset = tmp$set, n = tmp$dim[["n"]])[remid,,drop = FALSE]
    }else{
      rmat <- ctbu(base[-c(1:sum(nodes)), -c(1:(h*tmp$dim[["ks"]])), drop = FALSE],
                   agg_order = agg_order,
                   agg_mat = ubam, sntz = !is.null(nn), tew = tew)
    }
    rownames(rmat) <- namesCS(n = NROW(rmat), names_vec = names_base,
                              names_list = dimnames(ubam))
    rmat
  })

  if(CCC){ # Mean LCC + BU
    out <- Reduce("+", lccmat)/length(lccmat)
  }else{   # Mean LCC
    out <- Reduce("+", lccmat[-length(lccmat)])/(length(lccmat)-1)
  }
  names(lccmat) <- paste0("L-", rep(1:(length(nodes)+1), each = length(kset)),
                          "_k-", rep(kset, (length(nodes)+1)))
  attr(out, "FoReco") <- list2env(list(lcc = lccmat,
                                       framework = "Cross-temporal",
                                       forecast_horizon = h,
                                       comb = comb,
                                       te_set = tmp$set,
                                       cs_n = unname(sum(dim(ubam))),
                                       balanced = sum(dim(ubam))==sum(dim(agg_mat)),
                                       rfun = "cslcc"))
  return(out)
}
