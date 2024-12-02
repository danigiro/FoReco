#' Optimal combination cross-sectional reconciliation
#'
#' @description
#' This function performs optimal (in least squares sense) combination cross-sectional
#' forecast reconciliation for a linearly constrained (e.g., hierarchical/grouped)
#' multiple time series (Wickramasuriya et al., 2019, Panagiotelis et al., 2022,
#' Girolimetto and Di Fonzo, 2023). The reconciled forecasts are calculated
#' using either a projection approach (Byron, 1978, 1979) or the equivalent structural
#' approach by Hyndman et al. (2011). Non-negative (Di Fonzo and Girolimetto, 2023) and
#' immutable (including Zhang et al., 2023) reconciled forecasts can be considered.
#'
#' @usage
#' csrec(base, agg_mat, cons_mat, comb = "ols", res = NULL, approach = "proj",
#'       nn = NULL, settings = NULL, bounds = NULL, immutable = NULL, ...)
#'
#' @param base A (\eqn{h \times n}) numeric matrix or multivariate time series (\code{mts} class)
#' containing the base forecasts to be reconciled; \eqn{h} is the forecast horizon, and \eqn{n} is
#' the total number of time series (\eqn{n = n_a + n_b}).
#' @inheritParams ctrec
#' @param res An (\eqn{N \times n}) optional numeric matrix containing the in-sample
#' residuals. This matrix is used to compute some covariance matrices.
#' @param comb A string specifying the reconciliation method. For a complete list, see [cscov].
#' @param bounds A matrix (see [set_bounds]) with 3 columns (\eqn{i,lower,upper}), such that
#' \itemize{
#'   \item Column 1 represents the cross-sectional series (\eqn{i = 1, \dots, n}).
#'   \item Columns 2 and 3 indicates the \emph{lower} and \emph{lower} bounds, respectively.
#' }
#' @param immutable A numeric vector containing the column indices of the base forecasts
#' (\code{base} parameter) that should be fixed.
#' @inheritDotParams cscov mse shrink_fun
#'
#' @returns A (\eqn{h \times n}) numeric matrix of cross-sectional reconciled forecasts.
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
#' Di Fonzo, T. and Girolimetto, D. (2023), Spatio-temporal reconciliation of solar
#' forecasts, \emph{Solar Energy}, 251, 13–29. \doi{10.1016/j.solener.2023.01.003}
#'
#' Girolimetto, D. and Di Fonzo, T. (2023), Point and probabilistic forecast reconciliation
#' for general linearly constrained multiple time series,
#' \emph{Statistical Methods & Applications}, 33, 581-607. \doi{10.1007/s10260-023-00738-6}.
#'
#' Hyndman, R.J., Ahmed, R.A., Athanasopoulos, G. and Shang, H.L. (2011),
#' Optimal combination forecasts for hierarchical time series,
#' \emph{Computational Statistics & Data Analysis}, 55, 9, 2579-2589.
#' \doi{10.1016/j.csda.2011.03.006}
#'
#' Panagiotelis, A., Athanasopoulos, G., Gamakumara, P. and Hyndman, R.J. (2021), Forecast
#' reconciliation: A geometric view with new insights on bias correction,
#' \emph{International Journal of Forecasting}, 37, 1, 343–359.
#' \doi{10.1016/j.ijforecast.2020.06.004}
#'
#' Stellato, B., Banjac, G., Goulart, P., Bemporad, A. and Boyd, S. (2020), OSQP:
#' An Operator Splitting solver for Quadratic Programs,
#' \emph{Mathematical Programming Computation}, 12, 4, 637-672.
#' \doi{10.1007/s12532-020-00179-2}
#'
#' Wickramasuriya, S.L., Athanasopoulos, G. and Hyndman, R.J. (2019), Optimal forecast
#' reconciliation for hierarchical and grouped time series through trace minimization,
#' \emph{Journal of the American Statistical Association}, 114, 526, 804-819.
#' \doi{10.1080/01621459.2018.1448825}
#'
#' Zhang, B., Kang, Y., Panagiotelis, A. and Li, F. (2023), Optimal reconciliation with
#' immutable forecasts, \emph{European Journal of Operational Research}, 308(2), 650–660.
#' \doi{10.1016/j.ejor.2022.11.035}
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
#' reco <- csrec(base = base, agg_mat = A, comb = "wls", res = res)
#'
#' # Zero constraints matrix for Z - X - Y = 0
#' C <- t(c(1, -1, -1))
#' reco <- csrec(base = base, cons_mat = C, comb = "wls", res = res) # same results
#'
#' # Non negative reconciliation
#' base[1,3] <- -base[1,3] # Making negative one of the base forecasts for variable Y
#' nnreco <- csrec(base = base, agg_mat = A, comb = "wls", res = res, nn = "osqp")
#' recoinfo(nnreco, verbose = FALSE)$info
#'
#' @family Reco: regression-based
#' @family Framework: cross-sectional
#' @export
#'
csrec <- function(base, agg_mat, cons_mat,
                  comb = "ols", res = NULL, approach = "proj",
                  nn = NULL, settings = NULL, bounds = NULL,
                  immutable = NULL, ...){

  # Check if either 'agg_mat' or 'cons_mat' is specified
  if(missing(agg_mat) && missing(cons_mat)){
    cli_abort("Argument {.arg agg_mat} (or {.arg cons_mat}) is missing,
              with no default.", call = NULL)
  } else if(!missing(agg_mat)){
    tmp <- cstools(agg_mat = agg_mat)
  } else {
    tmp <- cstools(cons_mat = cons_mat)
  }

  n <- tmp$dim[["n"]]
  strc_mat <- tmp$strc_mat
  cons_mat <- tmp$cons_mat

  if(is.null(tmp$agg_mat)){
    id_nn <- NULL
  }else{
    id_nn <- c(rep(0, tmp$dim[["na"]]), rep(1, tmp$dim[["nb"]]))
  }

  # Check if 'base' is provided and its dimensions match with the data
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  } else if(NCOL(base) == 1){
    base <- t(base)
  }

  if(NCOL(base) != n){
    cli_abort("Incorrect {.arg base} columns dimension.", call = NULL)
  }

  # Check immutable
  if(!is.null(immutable)){
    if(NCOL(immutable) != 1){
      cli_abort("{.arg immutable} is not a vector.", call = NULL)
    }

    if(max(immutable) > n){
      cli_abort("{.code max(immutable)} must be less or equal to {n}", call = NULL)
    }

    if(length(immutable) >= n){
      # Answer issue: https://github.com/danigiro/FoReco/issues/6#issue-2397642027 (@AngelPone)
      cli_abort("{.code length(immutable)} must be less than {n}", call = NULL)
    }
  }

  # Check bounds cs
  if(!is.null(bounds)){
    if(is.vector(bounds)){
      bounds <- matrix(bounds, ncol = length(bounds))
    }
    if(NCOL(bounds) != 3){
      cli_abort("{.arg bounds} is not a matrix with 3 columns.", call = NULL)
    }
    bounds_approach <- attr(bounds, "approach")

    bounds <- bounds[bounds[,1] <= n, , drop = FALSE]

    if(is.null(bounds)){
      cli_warn("No valid bounds", call = NULL)
    }else{
      attr(bounds, "approach") <- bounds_approach
    }
  }

  # Compute covariance
  if(is(comb, "Matrix") | is(comb, "matrix")){
    cov_mat <- comb
  }else{
    cov_mat <- cscov(comb = comb, n = n, agg_mat = agg_mat, res = res, strc_mat = strc_mat, ...)
  }
  if(NROW(cov_mat) != n | NCOL(cov_mat) != n){
    cli_abort(c("Incorrect covariance dimensions.",
                "i"="Check {.arg res} columns dimension."), call = NULL)
  }

  reco_mat <- reco(base = base,
                   cov_mat = cov_mat,
                   strc_mat = strc_mat,
                   cons_mat = cons_mat,
                   approach = approach,
                   nn = nn,
                   immutable = immutable,
                   id_nn = id_nn,
                   bounds = bounds,
                   settings = settings)

  if(missing(agg_mat)){
    colnames(reco_mat) <- namesCS(n = NCOL(reco_mat), names_vec = colnames(base))
  }else{
    colnames(reco_mat) <- namesCS(n = NCOL(reco_mat), names_vec = colnames(base),
                                  names_list = dimnames(agg_mat))
  }

  rownames(reco_mat) <- paste0("h-", 1:NROW(reco_mat))
  attr(reco_mat, "FoReco") <- list2env(list(info = attr(reco_mat, "info"),
                                            framework = "Cross-sectional",
                                            forecast_horizon = NROW(reco_mat),
                                            comb = comb,
                                            cs_n = n,
                                            rfun = "csrec"))
  attr(reco_mat, "info") <- NULL
  return(reco_mat)
}

#' Optimal combination temporal reconciliation
#'
#' @description
#' This function performs forecast reconciliation for a single time series using temporal
#' hierarchies (Athanasopoulos et al., 2017, Nystrup et al., 2020). The reconciled forecasts can be computed
#' using either a projection approach (Byron, 1978, 1979) or the equivalent structural
#' approach by Hyndman et al. (2011). Non-negative (Di Fonzo and Girolimetto, 2023)
#' and immutable reconciled forecasts can be considered.
#'
#' @usage
#' terec(base, agg_order, comb = "ols", res = NULL, tew = "sum",
#'       approach = "proj", nn = NULL, settings = NULL, bounds = NULL,
#'       immutable = NULL, ...)
#'
#' @param base A (\eqn{h(k^\ast + m) \times 1}) numeric vector containing base forecasts
#' to be reconciled ordered from the lowest frequency to the highest frequency; \eqn{m}
#' is the max aggregation order, \eqn{k^\ast} is the sum of (a subset of) (\eqn{p-1})
#' factors of \eqn{m}, excluding \eqn{m}, and \eqn{h} is the forecast horizon for the
#' lowest frequency time series.
#' @param comb A string specifying the reconciliation method. For a complete list, see [tecov].
#' @param res A (\eqn{N(k^\ast+m) \times 1}) optional numeric vector containing the
#' in-sample residuals at all the temporal frequencies ordered from the lowest frequency
#' to the highest frequency. This vector is used to compute come covariance matrices.
#' @inheritParams ctrec
#' @param bounds A matrix (see [set_bounds]) with 4 columns (\eqn{k,j,lower,upper}), such that
#' \itemize{
#'   \item Column 1 represents the temporal aggregation order (\eqn{k = m,\dots,1}).
#'   \item Column 2 represents the temporal forecast horizon (\eqn{j = 1,\dots,m/k}).
#'   \item Columns 3 and 4 indicates the \emph{lower} and \emph{lower} bounds, respectively.
#' }
#' @param immutable A matrix with 2 columns (\eqn{k,j}), such that
#' \itemize{
#'   \item Column 1 represents the temporal aggregation order (\eqn{k = m,\dots,1}).
#'   \item Column 2 represents the temporal forecast horizon (\eqn{j = 1,\dots,m/k}).
#' }
#' For example, when working with a quarterly time series:
#' \itemize{
#'   \item \code{t(c(4, 1))} - Fix the one step ahead annual forecast.
#'   \item \code{t(c(1, 2))} - Fix the two step ahead quarterly forecast.
#' }
#' @inheritDotParams tecov mse shrink_fun
#'
#' @return A (\eqn{h(k^\ast+m) \times 1}) numeric vector of temporal reconciled forecasts.
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
#' Byron, R.P. (1979), Corrigenda: The estimation of large social account matrices,
#' \emph{Journal of the Royal Statistical Society, Series A}, 142(3), 405.
#' \doi{10.2307/2982515}
#'
#' Di Fonzo, T. and Girolimetto, D. (2023), Spatio-temporal reconciliation of solar
#' forecasts, \emph{Solar Energy}, 251, 13–29. \doi{10.1016/j.solener.2023.01.003}
#'
#' Hyndman, R.J., Ahmed, R.A., Athanasopoulos, G. and Shang, H.L. (2011),
#' Optimal combination forecasts for hierarchical time series,
#' \emph{Computational Statistics & Data Analysis}, 55, 9, 2579-2589.
#' \doi{10.1016/j.csda.2011.03.006}
#'
#' Nystrup, P.,  \enc{Lindström}{Lindstrom}, E., Pinson, P. and Madsen, H. (2020),
#' Temporal hierarchies with autocorrelation for load forecasting,
#' \emph{European Journal of Operational Research}, 280, 1, 876-888.
#' \doi{10.1016/j.ejor.2019.07.061}
#'
#' Stellato, B., Banjac, G., Goulart, P., Bemporad, A. and Boyd, S. (2020), OSQP:
#' An Operator Splitting solver for Quadratic Programs,
#' \emph{Mathematical Programming Computation}, 12, 4, 637-672.
#' \doi{10.1007/s12532-020-00179-2}
#'
#' @examples
#' set.seed(123)
#' # (7 x 1) base forecasts vector (simulated), m = 4
#' base <- rnorm(7, rep(c(20, 10, 5), c(1, 2, 4)))
#' # (70 x 1) in-sample residuals vector (simulated)
#' res <- rnorm(70)
#'
#' m <- 4 # from quarterly to annual temporal aggregation
#' reco <- terec(base = base, agg_order = m, comb = "wlsv", res = res)
#'
#' # Immutable reconciled forecast
#' # E.g. fix all the quarterly forecasts
#' imm_q <- expand.grid(k = 1, j = 1:4)
#' immreco <- terec(base = base, agg_order = m, comb = "wlsv",
#'                  res = res, immutable = imm_q)
#'
#' # Non negative reconciliation
#' base[7] <- -base[7] # Making negative one of the quarterly base forecasts
#' nnreco <- terec(base = base, agg_order = m, comb = "wlsv",
#'                 res = res, nn = "osqp")
#' recoinfo(nnreco, verbose = FALSE)$info
#'
#' @family Reco: regression-based
#' @family Framework: temporal
#' @export
#'
terec <- function(base, agg_order, comb = "ols", res = NULL, tew = "sum",
                  approach = "proj", nn = NULL, settings = NULL, bounds = NULL,
                  immutable = NULL, ...){

  # Check if 'agg_order' is provided
  if(missing(agg_order)){
    cli_abort("Argument {.arg agg_order} is missing, with no default.", call = NULL)
  }

  tmp <- tetools(agg_order = agg_order, tew = tew)
  kset <- tmp$set
  m <- tmp$dim[["m"]]
  kt <- tmp$dim[["kt"]]
  id_nn <- c(rep(0, tmp$dim[["ks"]]), rep(1, m))
  strc_mat <- tmp$strc_mat
  cons_mat <- tmp$cons_mat

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

  # Check immutable
  if(!is.null(immutable)){
    if(is.vector(immutable)){
      immutable <- matrix(immutable, ncol = length(immutable))
    }

    if(NCOL(immutable) != 2){
      cli_abort("{.arg immutable} is not a matrix with 2 columns.", call = NULL)
    }

    immutable <- immutable[immutable[,1] %in% kset, , drop = FALSE]
    immutable <- immutable[immutable[,2] <= m/immutable[,1], , drop = FALSE]
    immutable <- apply(immutable, 1,  function(x){
      if(x[1] %in% kset & x[2] <= m/x[1]){
        which(rep(kset, m/kset) == x[1] &
                do.call(c, sapply(m/kset, seq.int)) == x[2])
      }
    })

    if(is.null(immutable)){
      cli_warn("No immutable forecasts", call = NULL)
    }else if(length(immutable) >= kt){
      # Answer issue: https://github.com/danigiro/FoReco/issues/6#issue-2397642027 (@AngelPone)
      cli_abort("The number of immutable constraints must be less than {kt}", call = NULL)
    }
  }

  # Check bounds te
  if(!is.null(bounds)){
    if(is.vector(bounds)){
      bounds <- matrix(bounds, ncol = length(bounds))
    }
    if(NCOL(bounds) != 4){
      cli_abort("{.arg bounds} is not a matrix with 4 columns.", call = NULL)
    }
    bounds_approach <- attr(bounds, "approach")

    bounds <- bounds[bounds[,1] %in% kset, , drop = FALSE]
    bounds <- bounds[bounds[,2] <= m/bounds[,1], , drop = FALSE]
    bounds <- t(apply(bounds, 1,  function(x){
      if(x[1] %in% kset & x[2] <= m/x[1]){
        c(which(rep(kset, m/kset) == x[1] &
                  do.call(c, sapply(m/kset, seq.int)) == x[2]), x[-c(1:2)])
      }
    }))

    if(is.null(bounds)){
      cli_warn("No valid bounds", call = NULL)
    }else{
      attr(bounds, "approach") <- bounds_approach
    }
  }

  # Compute covariance
  if(is(comb, "Matrix") | is(comb, "matrix")){
    cov_mat <- comb
  }else{
    cov_mat <- tecov(comb = comb, res = res, agg_order = kset, tew = tew, strc_mat = strc_mat, ...)
  }

  if(NROW(cov_mat) != kt | NCOL(cov_mat) != kt){
    cli_abort(c("Incorrect covariance dimensions.",
                "i"="Check {.arg res} length."), call = NULL)
  }

  reco_mat <- reco(base = base,
                   cov_mat = cov_mat,
                   strc_mat = strc_mat,
                   cons_mat = cons_mat,
                   approach = approach,
                   nn = nn,
                   immutable = immutable,
                   id_nn = id_nn,
                   bounds = bounds,
                   settings = settings)

  # Convert 'reco_mat' back to vector
  out <- hmat2vec(reco_mat, h = h, kset = kset)
  attr(out, "FoReco") <- list2env(list(info = attr(reco_mat, "info"),
                                       framework = "Temporal",
                                       forecast_horizon = h,
                                       comb = comb,
                                       te_set = tmp$set,
                                       rfun = "terec"))
  return(out)
}

#' Optimal combination cross-temporal reconciliation
#'
#' @description
#' This function performs optimal (in least squares sense) combination cross-temporal forecast
#' reconciliation (Di Fonzo and Girolimetto 2023a, Girolimetto et al. 2023). The reconciled
#' forecasts are calculated using either a projection approach (Byron, 1978, 1979) or the
#' equivalent structural approach by Hyndman et al. (2011). Non-negative (Di Fonzo and
#' Girolimetto, 2023) and immutable reconciled forecasts can be considered.
#'
#' @usage
#' ctrec(base, agg_mat, cons_mat, agg_order, comb = "ols", res = NULL,
#'       tew = "sum", approach = "proj", nn = NULL, settings = NULL,
#'       bounds = NULL, immutable = NULL, ...)
#'
#' @param base A (\eqn{n \times h(k^\ast+m)}) numeric matrix containing the base forecasts to
#' be reconciled; \eqn{n} is the total number of variables, \eqn{m} is the max. order of temporal
#' aggregation, \eqn{k^\ast} is the sum of (a subset of) (\eqn{p-1}) factors of \eqn{m},
#' excluding \eqn{m}, and \eqn{h} is the forecast horizon for the lowest frequency time series.
#' The row identifies a time series, and the forecasts in each row are ordered from the
#' lowest frequency (most temporally aggregated) to the highest frequency.
#' @param agg_mat A (\eqn{n_a \times n_b}) numeric matrix representing the cross-sectional
#' aggregation matrix. It maps the \eqn{n_b} bottom-level (free)
#' variables into the \eqn{n_a} upper (constrained) variables.
#' @param cons_mat A (\eqn{n_a \times n}) numeric matrix representing the cross-sectional
#' zero constraints: each row represents a constraint equation, and each column represents
#' a variable. The matrix can be of full rank, meaning the rows are linearly independent,
#' but this is not a strict requirement, as the function allows for redundancy in the
#' constraints.
#' @param agg_order Highest available sampling frequency per seasonal cycle (max. order
#' of temporal aggregation, \eqn{m}), or a vector representing a subset of \eqn{p} factors
#' of \eqn{m}.
#' @param comb A string specifying the reconciliation method. For a complete list, see [ctcov].
#' @param res A (\eqn{n \times N(k^\ast+m)}) optional numeric matrix containing the
#' in-sample residuals at all the temporal frequencies ordered from the lowest frequency
#' to the highest frequency (columns) for each variable (rows). This matrix is used
#' to compute some covariance matrices.
#' @param tew A string specifying the type of temporal aggregation. Options include:
#' "\code{sum}" (simple summation, \emph{default}), "\code{avg}" (average),
#' "\code{first}" (first value of the period), and "\code{last}"
#' (last value of the period).
#' @param approach A string specifying the approach used to compute the reconciled
#' forecasts. Options include:
#'   \itemize{
#'   \item "\code{proj}" (\emph{default}): Projection approach according to Byron (1978, 1979).
#'   \item "\code{strc}": Structural approach as proposed by Hyndman et al. (2011).
#'   \item "\code{proj_osqp}": Numerical solution using \href{https://osqp.org/}{\pkg{osqp}}
#'   for projection approach.
#'   \item "\code{strc_osqp}": Numerical solution using \href{https://osqp.org/}{\pkg{osqp}}
#'   for structural approach.
#'   }
#' @param nn A string specifying the algorithm to compute non-negative forecasts:
#'   \itemize{
#'   \item "\code{osqp}": quadratic programming optimization
#'   (\href{https://osqp.org/}{\pkg{osqp}} solver).
#'   \item "\code{bpv}": block principal pivoting algorithm.
#'   \item "\code{sntz}": heuristic "set-negative-to-zero" (Di Fonzo and Girolimetto, 2023).
#'   }
#' @param settings A list of control parameters.
#'   \itemize{
#'   \item \code{nn = "osqp"} An object of class \code{osqpSettings} specifying settings
#'   for the \href{https://osqp.org/}{\pkg{osqp}} solver. For details, refer to the
#'   \href{https://osqp.org/}{\pkg{osqp} documentation} (Stellato et al., 2020).
#'   \item \code{nn = "bpv"} It includes: \code{ptype} for permutation method ("\code{random}"
#'   or "\code{fixed}", \emph{default}), \code{par} for the number of full exchange rules that
#'   may be attempted (\code{10}, \emph{default}), \code{tol} for the tolerance in convergence
#'   criteria (\code{sqrt(.Machine$double.eps)}, \emph{default}), \code{gtol} for the gradient
#'   tolerance in convergence criteria (\code{sqrt(.Machine$double.eps)}, \emph{default}),
#'   \code{itmax} for the maximum number of algorithm iterations (\code{100}, \emph{default})
#'   }
#' @param bounds A matrix (see [set_bounds]) with 5 columns (\eqn{i,k,j,lower,upper}), such that
#' \itemize{
#'   \item Column 1 represents the cross-sectional series (\eqn{i = 1, \dots, n}).
#'   \item Column 2 represents the temporal aggregation order (\eqn{k = m,\dots,1}).
#'   \item Column 3 represents the temporal forecast horizon (\eqn{j = 1,\dots,m/k}).
#'   \item Columns 4 and 5 indicates the \emph{lower} and \emph{lower} bounds, respectively.
#' }
#' @param immutable A matrix with three columns (\eqn{i,k,j}), such that
#' \itemize{
#'   \item Column 1 represents the cross-sectional series (\eqn{i = 1, \dots, n}).
#'   \item Column 2 represents the temporal aggregation order (\eqn{k = m,\dots,1}).
#'   \item Column 3 represents the temporal forecast horizon (\eqn{j = 1,\dots,m/k}).
#' }
#' For example, when working with a quarterly multivariate time series (\eqn{n = 3}):
#' \itemize{
#'   \item \code{t(c(1, 4, 1))} - Fix the one step ahead annual forecast of the first time series.
#'   \item \code{t(c(2, 1, 2))} - Fix the two step ahead quarterly forecast of the second time series.
#' }
#' @inheritDotParams ctcov mse shrink_fun
#'
#' @return A (\eqn{n \times h(k^\ast+m)}) numeric matrix of cross-temporal reconciled forecasts.
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
#' Di Fonzo, T. and Girolimetto, D. (2023a), Cross-temporal forecast reconciliation:
#' Optimal combination method and heuristic alternatives, \emph{International Journal
#' of Forecasting}, 39, 1, 39-57. \doi{10.1016/j.ijforecast.2021.08.004}
#'
#' Di Fonzo, T. and Girolimetto, D. (2023), Spatio-temporal reconciliation of solar
#' forecasts, \emph{Solar Energy}, 251, 13–29. \doi{10.1016/j.solener.2023.01.003}
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
#' Stellato, B., Banjac, G., Goulart, P., Bemporad, A. and Boyd, S. (2020), OSQP:
#' An Operator Splitting solver for Quadratic Programs,
#' \emph{Mathematical Programming Computation}, 12, 4, 637-672.
#' \doi{10.1007/s12532-020-00179-2}
#'
#' @examples
#' set.seed(123)
#' # (3 x 7) base forecasts matrix (simulated), Z = X + Y and m = 4
#' base <- rbind(rnorm(7, rep(c(20, 10, 5), c(1, 2, 4))),
#'               rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))),
#'               rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))))
#' # (3 x 70) in-sample residuals matrix (simulated)
#' res <- rbind(rnorm(70), rnorm(70), rnorm(70))
#'
#' A <- t(c(1,1)) # Aggregation matrix for Z = X + Y
#' m <- 4 # from quarterly to annual temporal aggregation
#' reco <- ctrec(base = base, agg_mat = A, agg_order = m, comb = "wlsv", res = res)
#'
#' C <- t(c(1, -1, -1)) # Zero constraints matrix for Z - X - Y = 0
#' reco <- ctrec(base = base, cons_mat = C, agg_order = m, comb = "wlsv", res = res)
#'
#' # Immutable reconciled forecasts
#' # Fix all the quarterly forecasts of the second variable.
#' imm_mat <- expand.grid(i = 2, k = 1, j = 1:4)
#' immreco <- ctrec(base = base, cons_mat = C, agg_order = m, comb = "wlsv",
#'                  res = res, immutable = imm_mat)
#'
#' # Non negative reconciliation
#' base[2,7] <- -2*base[2,7] # Making negative one of the quarterly base forecasts for variable X
#' nnreco <- ctrec(base = base, cons_mat = C, agg_order = m, comb = "wlsv",
#'                 res = res, nn = "osqp")
#' recoinfo(nnreco, verbose = FALSE)$info
#'
#' @family Reco: regression-based
#' @family Framework: cross-temporal
#' @export
#'
ctrec <- function(base, agg_mat, cons_mat, agg_order, comb = "ols", res = NULL,
                  tew = "sum", approach = "proj", nn = NULL, settings = NULL,
                  bounds = NULL, immutable = NULL, ...){

  # Check if 'base' is provided and its dimensions match with the data
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }

  # Check if 'agg_order' is provided
  if(missing(agg_order)){
    cli_abort("Argument {.arg agg_order} is missing, with no default.", call = NULL)
  }

  # Check if either 'agg_mat' or 'cons_mat' is specified
  if(missing(agg_mat) && missing(cons_mat)){
    cli_abort("Argument {.arg agg_mat} (or {.arg cons_mat}) is missing,
              with no default.", call = NULL)
  }else if(!missing(agg_mat)){
    tmp <- cttools(agg_mat = agg_mat, agg_order = agg_order, tew = tew)
    strc_mat <- tmp$strc_mat
    cons_mat <- tmp$cons_mat
  }else{
    tmp <- cttools(cons_mat = cons_mat, agg_order = agg_order, tew = tew)
    strc_mat <- tmp$strc_mat
    cons_mat <- tmp$cons_mat
    agg_mat <- cstools(cons_mat = cons_mat)$agg_mat
  }

  if(is.null(strc_mat)){
    id_nn <- NULL
  }else{
    cs_nn <- c(rep(0, tmp$dim[["na"]]), rep(1, tmp$dim[["nb"]]))
    te_nn <- c(rep(0, tmp$dim[["ks"]]), rep(1, tmp$dim[["m"]]))
    id_nn <- as.numeric(kronecker(cs_nn, te_nn))
  }

  if(NCOL(base) %% tmp$dim[["kt"]] != 0){
    cli_abort("Incorrect {.arg base} columns dimension.", call = NULL)
  }

  if(NROW(base) != tmp$dim[["n"]]){
    cli_abort("Incorrect {.arg base} rows dimension.", call = NULL)
  }

  # Check immutable
  if(!is.null(immutable)){
    if(is.vector(immutable)){
      immutable <- matrix(immutable, ncol = length(immutable))
    }
    if(NCOL(immutable) != 3){
      cli_abort("{.arg immutable} is not a matrix with 3 columns.", call = NULL)
    }

    immutable <- immutable[immutable[,1] <= tmp$dim[["n"]], , drop = FALSE]
    immutable <- immutable[immutable[,2] %in% tmp$set, , drop = FALSE]
    immutable <- immutable[immutable[,3] <= tmp$dim[["m"]]/immutable[,2], , drop = FALSE]

    if(NROW(immutable) == 0){
      cli_warn("No immutable forecasts.", call = NULL)
      immutable <- NULL
    }else{
      imm_mat <- sparseMatrix(integer(), integer(), x = numeric(),
                              dims = c(tmp$dim[["n"]], tmp$dim[["kt"]]))
      col_id <- apply(immutable, 1, function(x){
        which(rep(tmp$set, tmp$dim[["m"]]/tmp$set) == x[2] &
                do.call(c, sapply(tmp$dim[["m"]]/tmp$set, seq.int)) == x[3])
      })
      imm_mat[immutable[,1], col_id] <- 1
      imm_vec <- as(t(imm_mat), "sparseVector")
      immutable <- imm_vec@i

      if(length(immutable) >= tmp$dim[["n"]]*tmp$dim[["kt"]]){
        # Answer issue: https://github.com/danigiro/FoReco/issues/6#issue-2397642027 (@AngelPone)
        cli_abort("The number of immutable constraints must be less than {tmp$dim[['n']]*tmp$dim[['kt']]}",
                  call = NULL)
      }
    }
  }

  # Check bounds ct
  if(!is.null(bounds)){
    if(is.vector(bounds)){
      bounds <- matrix(bounds, ncol = length(bounds))
    }
    if(NCOL(bounds) != 5){
      cli_abort("{.arg bounds} is not a matrix with 5 columns.", call = NULL)
    }
    bounds_approach <- attr(bounds, "approach")

    bounds <- bounds[bounds[,1] <= tmp$dim[["n"]], , drop = FALSE]
    bounds <- bounds[bounds[,2] %in% tmp$set, , drop = FALSE]
    bounds <- bounds[bounds[,3] <= tmp$dim[["m"]]/bounds[,2], , drop = FALSE]

    if(NROW(bounds) == 0){
      cli_warn("No valid bounds.", call = NULL)
      bounds <- NULL
    }else{
      bounds_id <- bounds[, 1:3]
      bounds_mat <- sparseMatrix(integer(), integer(), x = numeric(),
                                 dims = c(tmp$dim[["n"]], tmp$dim[["kt"]]))
      col_id <- apply(bounds, 1, function(x){
        which(rep(tmp$set, tmp$dim[["m"]]/tmp$set) == x[2] &
                do.call(c, sapply(tmp$dim[["m"]]/tmp$set, seq.int)) == x[3])
      })
      bounds_mat[bounds[,1], col_id] <- 1
      bounds_vec <- as(t(bounds_mat), "sparseVector")
      bounds_id <- bounds_vec@i

      bounds <- cbind(bounds_id, bounds[, -c(1:3), drop = FALSE])
      attr(bounds, "approach") <- bounds_approach
    }
  }

  # Calculate 'h' and 'base_hmat'
  h <- NCOL(base) / tmp$dim[["kt"]]
  base_hmat <- mat2hmat(base, h = h, kset = tmp$set, n = tmp$dim[["n"]])

  # Compute covariance
  if(is(comb, "Matrix") | is(comb, "matrix")){
    cov_mat <- comb
  }else{
    cov_mat <- ctcov(comb = comb, res = res, agg_order = tmp$set, agg_mat = agg_mat,
                     n = tmp$dim[["n"]], tew = tew, strc_mat = strc_mat, ...)
  }
  if(NROW(cov_mat) != prod(tmp$dim[c("kt", "n")]) | NCOL(cov_mat) != prod(tmp$dim[c("kt", "n")])){
    cli_abort(c("Incorrect covariance dimensions.",
                "i"="Check {.arg res} dimensions."), call = NULL)
  }

  reco_mat <- reco(base = base_hmat,
                   cov_mat = cov_mat,
                   strc_mat = strc_mat,
                   cons_mat = cons_mat,
                   approach = approach,
                   nn = nn,
                   immutable = immutable,
                   id_nn = id_nn,
                   bounds = bounds,
                   settings = settings)

  # Convert 'reco_mat' back to matrix
  out <- hmat2mat(reco_mat, h = h, kset = tmp$set, n = tmp$dim[["n"]])
  rownames(out) <- namesCS(n = NROW(out), names_vec = rownames(base),
                           names_list = dimnames(agg_mat))

  attr(out, "FoReco") <- list2env(list(info = attr(reco_mat, "info"),
                                       framework = "Cross-temporal",
                                       forecast_horizon = h,
                                       comb = comb,
                                       te_set = tmp$set,
                                       cs_n = tmp$dim[["n"]],
                                       rfun = "ctrec"))
  return(out)
}
