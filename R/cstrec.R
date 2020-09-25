#' @title Heuristic first-cross-sectional-then-temporal cross-temporal forecast reconciliation
#'
#' @description
#' The order of application of the two reconciliation steps proposed by Kourentzes and Athanasopoulos (2019),
#' implemented in the function \code{\link[FoReco]{tcsrec}}, may be inverted. The function
#' \code{\link[FoReco]{cstrec}} performs cross-sectional reconciliation (\code{\link[FoReco]{htsrec}})
#' first, then temporal reconciliation (\code{\link[FoReco]{thfrec}}), and finally applies the
#' average of the projection matrices obtained in the second step to the one dimensional reconciled
#' values obtained in the first step.
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
#' \eqn{(\textbf{U}'\textbf{Y} = \mathbf{0})}{} spanning the null space valid for the reconciled
#' forecasts. It can be used instead of parameter \code{C}, but in this case \code{nb} (n = na + nb) is needed. If
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
#' temporal frequencies, ordered [lowest_freq' ...  highest_freq']' (columns) for
#' each variable (row), needed to estimate the covariance matrix when \code{hts_comb =}
#' \code{\{"wls",} \code{"shr",} \code{"sam"\}} and/or \code{hts_comb =} \code{\{"wlsv",}
#' \code{"wlsh",} \code{"acov",} \code{"strar1",} \code{"sar1",} \code{"har1",}
#' \code{"shr",} \code{"sam"\}}. The rows must be in the same order as \code{basef}.
#' @param mse Logical value: \code{TRUE} (\emph{default}) calculates the
#' covariance matrix of the in-sample residuals (when necessary) according to the original
#' \pkg{hts} and \pkg{thief} formulation: no mean correction, T as denominator.
#' @param corpcor Logical value: \code{TRUE} if \pkg{corpcor} (\enc{Schäfer}{Schafer} et
#' al., 2017) must be used to shrink the sample covariance matrix according to
#' \enc{Schäfer}{Schafer} and Strimmer (2005), otherwise the function uses the same
#' implementation as package \pkg{hts}.
#' @param nn Logical value, \code{TRUE} if non-negative reconciled forecasts are wished. \strong{Warning},
#' the two-step heuristic reconciliation allows non negativity constraints only in the first step.
#' This means that non-negativity is not guaranteed in the final reconciled values.
#' @param settings Settings for \pkg{osqp} (object \code{\link[osqp]{osqpSettings}}). The default options
#' are: \code{verbose = FALSE}, \code{eps_abs = 1e-5}, \code{eps_rel = 1e-5},
#' \code{polish_refine_iter = 100} and \code{polish = TRUE}. For details, see the
#' \href{https://osqp.org/}{\pkg{osqp} documentation} (Stellato et al., 2019).
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
#' obj <- cstrec(FoReco_data$base, m = 12, C = FoReco_data$C, thf_comb = "acov",
#'               hts_comb = "shr", res = FoReco_data$res)
#'
#' @usage cstrec(basef, m, C, thf_comb, hts_comb, Ut, nb, res, W, Omega,
#'        mse = TRUE, corpcor = FALSE, nn = FALSE,
#'        settings = osqpSettings(verbose = FALSE, eps_abs = 1e-5,
#'        eps_rel = 1e-5, polish_refine_iter = 100, polish = TRUE))
#'
#' @export
cstrec <- function(basef, m, C, thf_comb, hts_comb, Ut, nb, res,
                   W, Omega, mse = TRUE, corpcor = FALSE, nn = FALSE,
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

  kset <- rev(divisors(m))
  kt <- sum(kset) # sum of p (include m) factors of m

  # Step 1
  if (missing(C)) {
    if (missing(Ut) | missing(nb)) {
      stop("Please, give C (or Ut AND nb)", call. = FALSE)
    }
    hts_mod <- function(...) {
      htsrec(
        Ut = Ut, nb = nb, mse = mse, corpcor = corpcor,
        type = "M", W = W, nn = nn,
        settings = settings, ...
      )
    }
  } else {
    hts_mod <- function(...) {
      htsrec(
        C = C, mse = mse, corpcor = corpcor,
        type = "M", W = W, nn = nn,
        settings = settings, ...
      )
    }
  }

  h <- NCOL(basef) / kt
  Y <- lapply(kset, function(x) basef[, rep(kset, rev(kset) * h) == x, drop = FALSE])

  if (missing(res)) {
    Y1 <- lapply(Y, function(x) t(hts_mod(basef = t(x), comb = hts_comb)$recf))
  } else {
    # Create list with lenght p, with time by time temporally reconcilked residuals matrices
    r <- NCOL(res) / kt
    E <- lapply(kset, function(x) res[, rep(kset, rev(kset) * r) == x, drop = FALSE])

    ## list of time by time cross sectional M matrix
    Y1 <- mapply(function(Y, E) t(hts_mod(basef = t(Y), comb = hts_comb, res = t(E))$recf),
      Y = Y, E = E
    )
  }
  Y1 <- do.call("cbind", Y1)

  # Step 2
  if (missing(res)) {
    M <- lapply(split(Y1, row(Y1)), function(x) {
      thfrec(x,
        m = m, comb = thf_comb, mse = mse,
        corpcor = corpcor, type = "M", Omega = Omega
      )$M
    })
  } else {
    M <- mapply(function(Y, X) {
      thfrec(Y,
        m = m, comb = thf_comb, res = X, mse = mse,
        corpcor = corpcor, type = "M", Omega = Omega
      )$M
    },
    Y = split(Y1, row(Y1)), X = split(res, row(res)), SIMPLIFY = FALSE
    )
  }

  ## Step 3: Cross-Temporal reconciled forecasts with heuristic
  meanM <- apply(simplify2array(lapply(M, as.matrix)), c(1, 2), sum) / length(M)

  if (h == 1) {
    Y3 <- Y1 %*% t(meanM)
  } else {
    h_pos <- rep(rep(1:h, length(kset)), rep(rev(kset), each = h))
    D <- Dmat(h = h, kset = kset, n = 1)

    Y3 <- list()
    for (i in 1:h) {
      Y3[[i]] <- Y1[, which(h_pos == i)] %*% t(meanM)
    }
    Y3 <- do.call("cbind", Y3)
    Y3 <- Y3 %*% D
  }

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
