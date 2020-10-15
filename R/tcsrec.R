#' @title Heuristic first-temporal-then-cross-sectional cross-temporal forecast reconciliation
#'
#' @description The cross-temporal forecast reconciliation procedure by
#' Kourentzes and Athanasopoulos (2019) can be viewed as an ensemble forecasting
#' procedure which exploits the simple averaging of different forecasts.
#' First, for each time series the forecasts at any temporal aggregation order are
#' reconciled using temporal hierarchies (\code{\link[FoReco]{thfrec}}), then
#' time-by-time cross-sectional reconciliation is performed (\code{\link[FoReco]{htsrec}}).
#' The projection matrices obtained at this step are then averaged and used to
#' cross-sectionally reconcile the forecasts obtained at step 1, by this way
#' fulfilling both cross-sectional and temporal constraints.
#'
#' @param basef  (\code{n x h(k* + m)}) matrix of base forecasts to be reconciled;
#' \code{n} is the total number of variables, \code{m} is the highest frequency,
#' \code{k*} is the sum of (\code{p-1}) factors of \code{m}, excluding \code{m},
#' and \code{h} is the forecast horizon. Each row identifies, a time series, and the forecasts
#' are ordered as [lowest_freq' ...  highest_freq']'.
#' @param m Highest available sampling frequency per seasonal cycle (max. order of temporal aggregation).
#' @param C (\code{na x nb}) cross-sectional (contemporaneous) matrix mapping the bottom
#' level series into the higher level ones.
#' @param Ut Zero constraints cross-sectional (contemporaneous) kernel matrix
#' \eqn{\textbf{U}'\textbf{Y} = \mathbf{0}}{} spanning the null space valid for the reconciled
#' forecasts. It can be used instead of parameter \code{C}, but needs \code{nb} (n = na + nb). If
#' the hierarchy admits a structural representation, \code{Ut} has dimension (\code{na x n}).
#' @param nb Number of bottom time series; if \code{C} is present, \code{nb} is not used.
#' @param thf_comb Type of the (\code{(k* + m) x (k* + m)}) covariance matrix to be used in
#' the temporal reconciliation, see more in \code{comb} param of \code{\link[FoReco]{thfrec}}.
#' @param hts_comb Type of the (\code{n x n}) covariance matrix to be used in the
#' cross-sectional reconciliation, see more in \code{comb} param of \code{\link[FoReco]{htsrec}}.
#' @param Omega This option permits to directly enter the covariance matrix in the
#' reconciliation through temporal hierarchies, see more in \code{Omega} param of \code{\link[FoReco]{thfrec}}.
#' @param W This option permits to directly enter the covariance matrix in the
#' cross-sectional reconciliation, see more in \code{W} param of \code{\link[FoReco]{htsrec}}.
#' @param res (\code{n x N(k* + m)}) matrix containing the residuals at all the
#' temporal frequencies ordered [lowest_freq' ...  highest_freq']' (columns) for
#' each variable (row), needed to estimate the covariance matrix when \code{hts_comb =}
#' \code{\{"wls",} \code{"shr",} \code{"sam"\}} and/or \code{hts_comb =} \code{\{"wlsv",}
#' \code{"wlsh",} \code{"acov",} \code{"strar1",} \code{"sar1",} \code{"har1",}
#' \code{"shr",} \code{"sam"\}}. The row must be in the same order as \code{basef}.
#' @param mse Logical value: \code{TRUE} (\emph{default}) calculates the
#' covariance matrix of the in-sample residuals (when necessary) according to the original
#' \pkg{hts} and \pkg{thief} formulation: no mean correction, T as denominator.
#' @param corpcor Logical value: \code{TRUE} if \pkg{corpcor} (\enc{Schäfer}{Schafer} et
#' al., 2017) must be used to shrink the sample covariance matrix according to
#' \enc{Schäfer}{Schafer} and Strimmer (2005), otherwise the function uses the same
#' implementation as package \pkg{hts}.
#' @param avg If \code{avg = "KA"} (\emph{default}), the final projection matrix \code{M} is the one proposed
#' by Kourentzes and Athanasopoulos (2019), otherwise it is calculated as simple average of
#' all the involved projection matrices at step 2 of th procedure (see Di Fonzo and
#' Girolimetto, 2020).
#' @param nn Logical value, \code{TRUE} if non-negative reconciled forecasts are wished. \strong{Warning},
#' the two-step heuristic reconciliation allows non negativity constraints only in the first step.
#' This means that non-negativity is not guaranteed in the final reconciled values.
#' @param settings Settings for \pkg{osqp} (object \code{\link[osqp]{osqpSettings}}). The default options
#' are: \code{verbose = FALSE}, \code{eps_abs = 1e-5}, \code{eps_rel = 1e-5},
#' \code{polish_refine_iter = 100} and \code{polish = TRUE}. For details, see the
#' \href{https://osqp.org/}{\pkg{osqp} documentation} (Stellato et al., 2019).
#'
#' @details
#' This function performs a two-step cross-temporal forecast reconciliation using
#' the covariance matrices chosen by the user. If the combinations used by Kourentzes and Athanasopoulos (2019) are
#' wished, \code{thf_comb} must be set equal to either \code{"struc"} or \code{"wlsv"},
#' and \code{hts_comb} equal to either \code{"shr"} or \code{"wls"}.
#'
#' @return
#' The function returns a list with two elements:
#' \item{\code{recf}}{(\code{n x h(k* + m)}) reconciled forecasts matrix.}
#' \item{\code{M}}{Projection matrix (projection approach).}
#'
#' @references
#' Di Fonzo, T., Girolimetto, D. (2020), Cross-Temporal Forecast Reconciliation:
#' Optimal Combination Method and Heuristic Alternatives, Department of Statistical
#' Sciences, University of Padua, \href{https://arxiv.org/abs/2006.08570}{arXiv:2006.08570}.
#'
#' Kourentzes, N., Athanasopoulos, G. (2019), Cross-temporal coherent forecasts
#' for Australian tourism, \emph{Annals of Tourism Research}, 75, 393-409.
#'
#' \enc{Schäfer}{Schafer}, J.L., Opgen-Rhein, R., Zuber, V., Ahdesmaki, M.,
#' Duarte Silva, A.P., Strimmer, K. (2017), \emph{Package `corpcor'}, R
#' package version 1.6.9 (April 1, 2017), \href{https://CRAN.R-project.org/package=corpcor}{https://CRAN.R-project.org/package= corpcor}.
#'
#' \enc{Schäfer}{Schafer}, J.L., Strimmer, K. (2005), A Shrinkage Approach to Large-Scale Covariance
#' Matrix Estimation and Implications for Functional Genomics, \emph{Statistical
#' Applications in Genetics and Molecular Biology}, 4, 1.
#'
#' Stellato, B., Banjac, G., Goulart, P., Bemporad, A., Boyd, S. (2018). OSQP:
#' An Operator Splitting Solver for Quadratic Programs, \href{https://arxiv.org/abs/1711.08013}{arXiv:1711.08013}.
#'
#' Stellato, B., Banjac, G., Goulart, P., Boyd, S., Anderson, E. (2019), OSQP:
#' Quadratic Programming Solver using the 'OSQP' Library, R package version 0.6.0.3
#' (October 10, 2019), \href{https://CRAN.R-project.org/package=osqp}{https://CRAN.R-project.org/package=osqp}.
#'
#' @keywords reconciliation heuristic
#' @examples
#' data(FoReco_data)
#' obj <- tcsrec(FoReco_data$base, m = 12, C = FoReco_data$C, thf_comb = "acov",
#'               hts_comb = "shr", res = FoReco_data$res)
#'
#' @export
#'
#' @usage tcsrec(basef, m, C, thf_comb, hts_comb, Ut, nb, res, W, Omega,
#'        mse = TRUE, corpcor = FALSE, avg="KA", nn= FALSE,
#'        settings = osqpSettings(verbose = FALSE, eps_abs = 1e-5,
#'        eps_rel = 1e-5, polish_refine_iter = 100, polish = TRUE))
tcsrec <- function(basef, m, C, thf_comb, hts_comb, Ut, nb, res, W, Omega,
                   mse = TRUE, corpcor = FALSE, avg = "KA", nn = FALSE,
                   settings = osqpSettings(
                     verbose = FALSE, eps_abs = 1e-5,
                     eps_rel = 1e-5, polish_refine_iter = 100,
                     polish = TRUE
                   )) {
  if (missing(basef)) {
    stop("the argument basef is not specified")
  }
  if (missing(m)) {
    stop("the argument m is not specified")
  }
  if (missing(thf_comb)) {
    stop("the argument thf_comb is not specified")
  }
  if (missing(hts_comb)) {
    stop("the argument hts_comb is not specified")
  }

  if (missing(W)) {
    W <- NULL
  }

  if (missing(Omega)) {
    Omega <- NULL
  }

  # k <- divisors(m) # set of (p-1) factors of m (exclude m as factors)
  kset <- rev(divisors(m))
  kt <- sum(kset) # sum of p (include m) factors of m

  ## Step 1: compute the temporally reconciled forecasts for each individual variable
  # (basef -> Y1)
  if (missing(res)) {
    Y1 <- t(apply(basef, 1, function(x) {
      thfrec(x,
        m = m, comb = thf_comb, mse = mse,
        corpcor = corpcor, type = "M", Omega = Omega,
        nn = nn, settings = settings
      )$recf
    }))
  } else {
    Y1 <- t(mapply(function(Y, X) {
      thfrec(Y,
        m = m, comb = thf_comb, res = X, mse = mse,
        corpcor = corpcor, type = "M", Omega = Omega,
        nn = nn, settings = settings
      )$recf
    },
    Y = split(basef, row(basef)), X = split(res, row(res))
    ))
  }

  ## Step 2: compute time by time cross sectional M matrix
  if (missing(C)) {
    if (missing(Ut) | missing(nb)) {
      stop("Please, give C (or Ut AND nb)", call. = FALSE)
    }
    hts_mod <- function(...) {
      htsrec(
        Ut = Ut, nb = nb, mse = mse, corpcor = corpcor,
        type = "M", W = W, ...
      )
    }
  } else {
    hts_mod <- function(...) {
      htsrec(
        C = C, mse = mse, corpcor = corpcor,
        type = "M", W = W, ...
      )
    }
  }

  # Create list with lenght p, with time by time temporally reconciled forecasts matrices
  h <- NCOL(Y1) / kt
  Y <- lapply(kset, function(x) Y1[, rep(kset, rev(kset) * h) == x, drop = FALSE])


  if (missing(res)) {
    M <- lapply(Y, function(x) hts_mod(basef = t(x), comb = hts_comb)$M)

    if (avg == "KA") {
      meanM <- apply(simplify2array(lapply(M, as.matrix)), c(1, 2), sum) / length(M)
    } else {
      Mw <- mapply(function(a, A) a * A, A = M, a = split(kset, 1:length(kset)), SIMPLIFY = FALSE)
      meanM <- apply(simplify2array(lapply(Mw, as.matrix)), c(1, 2), sum) / sum(kset)
    }
  } else {
    # Create list with lenght p, with time by time temporally reconciled residuals matrices
    r <- NCOL(res) / kt
    E <- lapply(kset, function(x) res[, rep(kset, rev(kset) * r) == x, drop = FALSE])

    ## list of time by time cross sectional M matrix
    M <- mapply(function(Y, E) hts_mod(basef = t(Y), comb = hts_comb, res = t(E))$M,
      Y = Y, E = E, SIMPLIFY = FALSE
    )

    if (avg == "KA") {
      meanM <- apply(simplify2array(lapply(M, as.matrix)), c(1, 2), sum) / length(M)
    } else {
      Mw <- mapply(function(a, A) a * A, A = M, a = split(kset, 1:length(kset)), SIMPLIFY = FALSE)
      meanM <- apply(simplify2array(lapply(Mw, as.matrix)), c(1, 2), sum) / sum(kset)
    }
  }

  ## Step 3: Cross-Temporal reconciled forecasts with heuristic
  Y3 <- meanM %*% Y1

  rec_sol <- list()
  rec_sol$recf <- Y3
  rownames(rec_sol$recf) <- if (is.null(rownames(basef))) paste("serie", 1:NROW(rec_sol$recf), sep = "") else rownames(basef)
  colnames(rec_sol$recf) <- paste("k", rep(kset, h * rev(kset)), "h",
    do.call("c", as.list(sapply(
      rev(kset) * h,
      function(x) seq(1:x)
    ))),
    sep = ""
  )
  rec_sol$M <- meanM
  rec_sol$nn_check <- sum(rec_sol$recf < 0)
  return(rec_sol)
}
