#' @title Iterative heuristic cross-temporal forecast reconciliation
#'
#' @description
#' Iterative procedure which produces cross-temporally reconciled
#' forecasts by alternating forecast reconciliation along one single dimension
#' (either cross-sectional or temporal) at each iteration step. \strong{Each iteration}
#' consists in the first two steps of the heuristic KA procedure, so the forecasts are
#' reconciled by alternating cross-sectional (contemporaneous) reconciliation, and
#' reconciliation through temporal hierarchies in a cyclic fashion.
#' The choice of the dimension along with the first reconciliation step in each
#' iteration is performed is up to the user (param \code{start_rec}), and there is
#' no particular reason why one should perform the temporal reconciliation first, and
#' the cross-sectional reconciliation then.
#' The iterative procedure allows the user to get non-negative reconciled forecasts.
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
#' @param type Approach used to compute the reconciled forecasts: \code{"M"} for
#' the projection approach with matrix M (\emph{default}), or \code{"S"} for the
#' structural approach with summing matrix S.
#' @param sol Solution technique for the reconciliation problem: either \code{"direct"} for the direct
#' solution or \code{"osqp"} for the numerical solution (solving a linearly constrained quadratic
#' program using \code{\link[osqp]{solve_osqp}}).
#' @param nn Logical value: \code{TRUE} if non-negative reconciled forecasts are wished.
#' @param settings Settings for \pkg{osqp} (object \code{\link[osqp]{osqpSettings}}). The default options
#' are: \code{verbose = FALSE}, \code{eps_abs = 1e-5}, \code{eps_rel = 1e-5},
#' \code{polish_refine_iter = 100} and \code{polish = TRUE}. For details, see the
#' \href{https://osqp.org/}{\pkg{osqp} documentation} (Stellato et al., 2019).
#' @param tol Convergence tolerance (\code{1e-5}, \emph{default}).
#' @param maxit Max number of iteration (\code{100}, \emph{default}).
#' @param start_rec Dimension along with the first reconciliation step in each
#' iteration is performed: it start from temporal reconciliation with "\code{thf}" (\emph{default}),
#' from cross-sectional with "\code{hts}" and it does both reconciliation with "\code{auto}".
#' @param note If \code{note = TRUE} (\emph{default}) the function writes some notes to the console,
#' otherwise no note is produced (also no plot).
#' @param plot Some useful plots: \code{"mti"} (\emph{default}) marginal trend inconsistencies,
#' \code{"pat"} step by step inconsistency pattern for each iteration, \code{"distf"} distance
#' forecasts iteration i and i-1, \code{"all"} all the plots.
#'
#' @details
#' The above procedure can be seen as an extension of the well known iterative proportional
#' fitting procedure (Deming and Stephan, 1940, Johnston and Pattie, 1993), also known as
#' RAS method (Miller and Blair, 2009), to adjust the internal cell values of a
#' two-dimensional matrix until they sum to some predetermined row and column
#' totals. In that case the adjustment follows a proportional adjustment scheme, whereas
#' in the cross-temporal reconciliation framework each adjustment step is made according to
#' the penalty function associated to the single-dimension reconciliation procedure adopted.
#'
#' Control status of iterative reconciliation:
#' \describe{
#' \item{\code{-2}}{Temporal/Cross-sectional reconciliation does not work.}
#' \item{\code{-1}}{Convergence not achieved (maximum iteration limit reached).}
#' \item{\code{0}}{Convergence achieved.}
#' \item{\code{+1}}{Convergence achieved: incoherence has increased in the next
#' iteration (at least one time).}
#' \item{\code{+2}}{Convergence achieved: incoherence has increased in the next
#' two or more iteration (at least one time).}
#' \item{\code{+3}}{The forecasts are already reconciled.}
#' }
#'
#' @return \code{iterec} returns a list with:
#' \item{\code{recf}}{(\code{n x h(k* + m)}) matrix of reconciled forecasts.}
#' \item{\code{d_cs}}{Cross-sectional incoherence at each iteration.}
#' \item{\code{d_te}}{Temporal incoherence at each iteration.}
#' \item{\code{start_rec}}{Starting coherence dimension (thf or hts).}
#' \item{\code{tol}}{Tolerance.}
#' \item{\code{flag}}{Control code (see \emph{details}).}
#' \item{\code{time}}{Elapsed time.}
#' \item{\code{dist}}{If \code{start_rec = "auto"}, matrix of distances of the forecasts reconciled
#' from the base.}
#'
#' @references
#' Deming, E., Stephan, F.F. (1940), On a least squares adjustment of a sampled frequency
#' table when the expected marginal totals are known, \emph{The Annals of Mathematical
#' Statistics}, 11, 4, 427–444.
#'
#' Di Fonzo, T., Girolimetto, D. (2020), Cross-Temporal Forecast Reconciliation:
#' Optimal Combination Method and Heuristic Alternatives, Department of Statistical
#' Sciences, University of Padua, \href{https://arxiv.org/abs/2006.08570}{arXiv:2006.08570}.
#'
#' Johnston, R.J., Pattie, C.J. (1993), Entropy-Maximizing and the Iterative Proportional
#' Fitting Procedure, \emph{The Professional Geographer}, 45, 3, 317–322.
#'
#' Kourentzes, N., Athanasopoulos, G. (2019), Cross-temporal coherent forecasts
#' for Australian tourism, \emph{Annals of Tourism Research}, 75, 393-409.
#'
#' Miller, R.E., Blair, P.D. (2009), Input-output analysis: foundations and extensions,
#' 2nd edition, New York, Cambridge University Press.
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
#' \donttest{
#' data(FoReco_data)
#' obj <- iterec(FoReco_data$base,
#'   m = 12, C = FoReco_data$C, thf_comb = "acov",
#'   hts_comb = "shr", res = FoReco_data$res, start_rec = "thf")
#' }
#'
#' @usage iterec(basef, m, C, thf_comb, hts_comb, Ut, nb, res, W,
#'        Omega, mse = TRUE, corpcor = FALSE,type = "M",
#'        sol = "direct", nn = FALSE, maxit = 100, tol = 1e-5,
#'        start_rec = "thf", note = TRUE, plot = "mti",
#'        settings = osqpSettings(verbose = FALSE, eps_abs = 1e-5,
#'        eps_rel = 1e-5, polish_refine_iter = 100, polish = TRUE))
#'
#' @export
iterec <- function(basef, m, C, thf_comb, hts_comb, Ut, nb, res, W, Omega, mse = TRUE,
                   corpcor = FALSE, type = "M", sol = "direct", nn = FALSE, maxit = 100,
                   tol = 1e-5, start_rec = "thf", note = TRUE, plot = "mti",
                   settings = osqpSettings(
                     verbose = FALSE, eps_abs = 1e-5, eps_rel = 1e-5,
                     polish_refine_iter = 100, polish = TRUE
                   )) {
  start_rec <- match.arg(start_rec, c("thf", "hts", "auto"))

  if (missing(res)) {
    res <- NULL
  }

  if (missing(W)) {
    W <- NULL
  }

  if (missing(Omega)) {
    Omega <- NULL
  }

  # Tools for cross-sectional reconciliation
  if (missing(C)) {
    if (missing(Ut) | missing(nb)) {
      stop("Please, give C (or Ut AND nb)", call. = FALSE)
    }
    r_u <- NROW(Ut)

    hts_mod <- function(...) {
      out <- htsrec(
        Ut = Ut, nb = nb, comb = hts_comb, mse = mse,
        corpcor = corpcor, W = W, type = type, keep = "recf",
        nn = nn, settings = settings, sol = sol, ...
      )
      if (is.list(out)) {
        return(out$recf)
      } else {
        return(out)
      }
    }
  } else {
    r_u <- NROW(C)
    Ut <- cbind(diag(1, r_u), -C)
    hts_mod <- function(...) {
      out <- htsrec(
        C = C, comb = hts_comb, mse = mse,
        corpcor = corpcor, W = W, type = type, keep = "recf",
        nn = nn, settings = settings, sol = sol, ...
      )
      if (is.list(out)) {
        return(out$recf)
      } else {
        return(out)
      }
    }
  }

  thf_mod <- function(...) {
    out <- thfrec(
      m = m, comb = thf_comb, mse = mse, corpcor = corpcor,
      type = type, nn = nn, settings = settings, sol = sol,
      Omega = Omega, ...
    )
    if (is.list(out)) {
      return(out$recf)
    } else {
      return(out)
    }
  }

  n <- NROW(basef)
  # Tools for cross-sectional reconciliation
  h <- ncol(basef) / thf_tools(m)$kt
  tmp <- thf_tools(m, h = h)
  kt <- tmp$kt
  ks <- tmp$ks
  kset <- tmp$kset
  Zt <- tmp$Zt

  # Incoherence vector
  d_cs_thf <- rep(NA, maxit + 1)
  d_te_thf <- rep(NA, maxit + 1)
  d_cs_thf[1] <- abs(t(rep(1, r_u)) %*% (Ut %*% basef) %*% rep(1, h * kt))
  d_te_thf[1] <- abs(t(rep(1, h * ks)) %*% (Zt %*% t(basef)) %*% rep(1, n))

  d_cs_hts <- rep(NA, maxit + 1)
  d_te_hts <- rep(NA, maxit + 1)
  d_cs_hts[1] <- abs(t(rep(1, r_u)) %*% (Ut %*% basef) %*% rep(1, h * kt))
  d_te_hts[1] <- abs(t(rep(1, h * ks)) %*% (Zt %*% t(basef)) %*% rep(1, n))

  # Check: reconciliation needed?
  if (d_cs_thf[1] < tol & d_te_thf[1] < tol) {
    warning("basef is already reconciled!")
    out <- list()
    out$recf <- basef
    out$flag <- 3
    return(out)
  }

  # Matrix with ALL incoherence at any step and any iteration
  d_thf <- t(c(0, d_cs_thf[1], d_te_thf[1]))
  d_hts <- t(c(0, d_cs_hts[1], d_te_hts[1]))

  # Choice of the starting reconciliation dimension
  if (start_rec == "thf") {
    start_thf <- TRUE
    start_hts <- FALSE

    # Settings
    Y2_thf <- basef
    flag_thf <- -1
    Y2_hts <- basef
    flag_hts <- -3
  } else if (start_rec == "hts") {
    start_thf <- FALSE
    start_hts <- TRUE

    # Settings
    Y2_thf <- basef
    flag_thf <- -3
    Y2_hts <- basef
    flag_hts <- -1
  } else {
    start_thf <- TRUE
    start_hts <- TRUE

    # Settings
    Y2_thf <- basef
    flag_thf <- -1
    Y2_hts <- basef
    flag_hts <- -1
  }

  dist_thf <- matrix(0, maxit, 5)
  Yi_thf <- Y2_thf
  colnames(dist_thf) <- c("iter", "norm2", "norm1", "rel2", "rel1")
  dist_hts <- matrix(0, maxit, 5)
  Yi_hts <- Y2_hts
  colnames(dist_hts) <- c("iter", "norm2", "norm1", "rel2", "rel1")

  # Time starting
  start_time <- Sys.time()
  if (note == TRUE) {
    message("---------------------------------------")
  }
  if (start_thf) {
    if (note == TRUE) {
      message("Iter #  (Cross-sectional incoherence) (Temporal incoherence)", sep = "")
      message("thf: ", formatC(0, width = 3, format = "d", flag = "0"), "  (", d_cs_thf[1], ")  (",
        d_te_thf[1], ")",
        sep = ""
      )
    }

    for (i in 1:maxit) {
      # Step1
      Y1_thf <- step_thf(basef = Y2_thf, res = res, thf_mod = thf_mod)
      d_cs_thf[i + 1] <- abs(t(rep(1, r_u)) %*% (Ut %*% Y1_thf) %*% rep(1, h * kt))
      check <- as.vector(abs(t(rep(1, h * ks)) %*% (Zt %*% t(Y1_thf)) %*% rep(1, n)))

      d_thf <- rbind(d_thf, c(i, d_cs_thf[i + 1], check))
      if (check > tol) {
        Y2_thf <- Y1_thf
        flag_thf <- -2
        warning(paste("thf: Program (starting from thf) stops at iteration number ", i,
          "! \nthf: Temporal reconciliation does not work. (", check, ")",
          sep = ""
        ), call. = FALSE)
        break
      } else if (i > 1) {
        if (d_cs_thf[(i + 1)] > d_cs_thf[i] & d_cs_thf[i] > d_cs_thf[(i - 1)] & flag_thf < 2) {
          flag_thf <- 2
        } else if (d_cs_thf[(i + 1)] > d_cs_thf[i] & flag_thf < 1) {
          flag_thf <- 1
        }
      }

      # Step 2
      Y2_thf <- step_hts(
        basef = Y1_thf, kset = kset, h = h,
        res = res, hts_mod = hts_mod
      )
      d_te_thf[i + 1] <- abs(t(rep(1, h * ks)) %*% (Zt %*% t(Y2_thf)) %*% rep(1, n))
      check <- as.vector(abs(t(rep(1, r_u)) %*% (Ut %*% Y2_thf) %*% rep(1, h * kt)))

      d_thf <- rbind(d_thf, c(i, check, d_te_thf[i + 1]))
      if (check > tol) {
        flag_thf <- -2
        warning(paste("thf: Program (starting from thf) stops at iteration number ", i,
          "! \nthf: Cross-sectional reconciliation does not work.(", check, ")",
          sep = ""
        ), call. = FALSE)
        break
      } else if (d_te_thf[i + 1] < tol) {
        if (flag_thf == -1) flag_thf <- 0
        if (note == TRUE) {
          message("thf: Convergence (starting from thf) achieved at iteration number ", i,
            "! \nthf: Temporal incoherence ", d_te_thf[i + 1], " < ", tol,
            " tolerance",
            sep = ""
          )
        }
        rel <- ((Y2_thf - Yi_thf) / abs(Yi_thf))
        dist_thf[i, ] <- c(
          i, sqrt(sum((Yi_thf - Y2_thf)^2)),
          sum(abs(Yi_thf - Y2_thf)),
          sum((rel * is.finite(rel))^2, na.rm = TRUE),
          sum(abs(rel * is.finite(rel)), na.rm = TRUE)
        )
        Yi_thf <- Y2_thf
        break
      } else if (i > 1) {
        if (d_te_thf[(i + 1)] > d_te_thf[i] & d_te_thf[i] > d_te_thf[(i - 1)] & flag_thf < 2) {
          flag_thf <- 2
        } else if (d_te_thf[(i + 1)] > d_te_thf[i] & flag_thf < 1) {
          flag_thf <- 1
        }
      }
      rel <- ((Y2_thf - Yi_thf) / abs(Yi_thf))
      dist_thf[i, ] <- c(
        i, sqrt(sum((Yi_thf - Y2_thf)^2)),
        sum(abs(Yi_thf - Y2_thf)),
        sum((rel * is.finite(rel))^2, na.rm = TRUE),
        sum(abs(rel * is.finite(rel)), na.rm = TRUE)
      )
      Yi_thf <- Y2_thf
      if (note == TRUE) {
        message("thf: ", formatC(i, width = 3, format = "d", flag = "0"), "  (", d_cs_thf[i + 1], ")  (",
          d_te_thf[i + 1], ")",
          sep = ""
        )
      }
    }

    if ((flag_thf == -1 | i == maxit) & start_rec != "hts") {
      flag_thf <- -1
      warning("thf: Convergence NOT achieved (starting from thf)! Maximum number of iterations reached (",
        maxit, ")",
        sep = "", call. = FALSE
      )
    }
    dist_thf <- dist_thf[1:i, , drop = FALSE]
    if (note == TRUE) {
      message("---------------------------------------")
    }
  }

  end_time_thf <- Sys.time()

  if (start_hts) {
    if (note == TRUE) {
      message("Iter #  (Cross-sectional incoherence) (Temporal incoherence)", sep = "")
      message("hts: ", formatC(0, width = 3, format = "d", flag = "0"), "  (", d_cs_hts[1], ")  (",
        d_te_hts[1], ")",
        sep = ""
      )
    }

    for (j in 1:maxit) {
      # Step 1
      Y1_hts <- step_hts(
        basef = Y2_hts, kset = kset, h = h,
        res = res, hts_mod = hts_mod
      )
      d_te_hts[j + 1] <- abs(t(rep(1, h * ks)) %*% (Zt %*% t(Y1_hts)) %*% rep(1, n))
      check <- as.vector(abs(t(rep(1, r_u)) %*% (Ut %*% Y1_hts) %*% rep(1, h * kt)))

      d_hts <- rbind(d_hts, c(j, check, d_te_hts[j + 1]))
      if (check > tol) {
        flag_hts <- -2
        warning(paste("hts: Program (starting from hts) stops at iteration number ", j,
          "! \nhts: Cross-sectional reconciliation does not work.(", check, ")",
          sep = ""
        ), call. = FALSE)
        break
      } else if (j > 1) {
        if (d_te_hts[(j + 1)] > d_te_hts[j] & d_te_hts[j] > d_te_hts[(j - 1)] & flag_hts < 2) {
          flag_hts <- 2
        } else if (d_te_hts[(j + 1)] > d_te_hts[j] & flag_hts < 1) {
          flag_hts <- 1
        }
      }

      # Step 2
      Y2_hts <- step_thf(basef = Y1_hts, res = res, thf_mod = thf_mod)
      d_cs_hts[j + 1] <- abs(t(rep(1, r_u)) %*% (Ut %*% Y2_hts) %*% rep(1, h * kt))
      check <- as.vector(abs(t(rep(1, h * ks)) %*% (Zt %*% t(Y2_hts)) %*% rep(1, n)))

      d_hts <- rbind(d_hts, c(j, d_cs_hts[j + 1], check))
      if (check > tol) {
        flag_hts <- -2
        warning(paste("hts: Program (starting from hts) stops at iteration number ", j,
          "! \nhts: Temporal reconciliation does not work.(", check, ")",
          sep = ""
        ))
        break
      } else if (d_cs_hts[j + 1] < tol) {
        if (flag_hts == -1) flag_hts <- 0
        if (note == TRUE) {
          message("hts: Convergence (starting from hts) achieved at iteration number ", j,
            "! \nhts: Cross-sectional incoherence ", d_cs_hts[j + 1], " < ", tol,
            " tolerance",
            sep = ""
          )
        }
        rel <- ((Y2_hts - Yi_hts) / abs(Yi_hts))
        dist_hts[j, ] <- c(
          j, sqrt(sum((Yi_hts - Y2_hts)^2)),
          sum(abs(Yi_hts - Y2_hts)),
          sum((rel * is.finite(rel))^2, na.rm = TRUE),
          sum(abs(rel * is.finite(rel)), na.rm = TRUE)
        )
        Yi_hts <- Y2_hts
        break
      } else if (j > 1) {
        if (d_cs_hts[(j + 1)] > d_cs_hts[j] & d_cs_hts[j] > d_cs_hts[(j - 1)] & flag_hts < 2) {
          flag_hts <- 2
        } else if (d_cs_hts[(j + 1)] > d_cs_hts[j] & flag_hts < 1) {
          flag_hts <- 1
        }
      }

      if (note == TRUE) {
        message("hts: ", formatC(j, width = 3, format = "d", flag = "0"), "  (", d_cs_hts[j + 1], ")  (",
          d_te_hts[j + 1], ")",
          sep = ""
        )
      }

      rel <- ((Y2_hts - Yi_hts) / abs(Yi_hts))
      dist_hts[j, ] <- c(
        j, sqrt(sum((Yi_hts - Y2_hts)^2)),
        sum(abs(Yi_hts - Y2_hts)),
        sum((rel * is.finite(rel))^2, na.rm = TRUE),
        sum(abs(rel * is.finite(rel)), na.rm = TRUE)
      )
      Yi_hts <- Y2_hts
    }

    if ((flag_hts == -1 | j == maxit) & start_rec != "thf") {
      flag_hts <- -1
      warning("hts: Convergence NOT achieved (starting from hts)! Maximum number of iterations reached (",
        maxit, ")",
        sep = "", call. = FALSE
      )
    }
    dist_hts <- dist_hts[1:j, , drop = FALSE]
    if (note == TRUE) {
      message("---------------------------------------")
    }
  }

  end_time_hts <- Sys.time()

  out <- list()
  if (start_rec == "thf") {
    tcs <- list()
    tcs$recf <- Y2_thf
    rownames(tcs$recf) <- if (is.null(rownames(basef))) paste("serie", 1:NROW(tcs$recf), sep = "") else rownames(basef)
    colnames(tcs$recf) <- paste("k", rep(kset, h * rev(kset)), "h",
      do.call("c", as.list(sapply(
        rev(kset) * h,
        function(x) seq(1:x)
      ))),
      sep = ""
    )
    tcs$d_cs <- d_cs_thf[!is.na(d_cs_thf)]
    tcs$d_te <- d_te_thf[!is.na(d_te_thf)]
    tcs$start_rec <- start_rec
    tcs$iter <- i
    tcs$tol <- tol
    tcs$flag <- flag_thf
    tcs$diff_iter <- dist_thf
    tcs$time <- as.numeric(end_time_thf - start_time, units = "secs")

    if (note == TRUE) {
      if (plot == "mti" | plot == "all") {
        plot(0, 0,
          type = "n", pch = 19, xlab = "iteration", ylab = "Incoherence",
          main = paste("Starting thf, flag = ", tcs$flag, ", tol = ", tol, sep = ""),
          ylim = c(0, max(tcs$d_cs, tcs$d_te)),
          xlim = c(0, i),
          xaxt = "n"
        )
        graphics::axis(1, at = 0:i)
        graphics::lines(x = 0:(length(tcs$d_cs) - 1), y = tcs$d_cs, type = "b", pch = 1, col = "#fe5f55")
        graphics::lines(x = 0:(length(tcs$d_te) - 1), y = tcs$d_te, type = "b", pch = 4, col = 1)
        if (tcs$flag == 0) {
          graphics::legend("topright",
            inset = 0.01, legend = c("Cross-sectional", "Temporal"),
            col = c("#fe5f55", 1), lty = 1, pch = c(1, 4), bty = "n"
          )
        }
        if (tcs$flag == 1) {
          graphics::legend("bottomright",
            inset = 0.01, legend = c("Cross-sectional", "Temporal"),
            col = c("#fe5f55", 1), lty = 1, pch = c(1, 4), bty = "n"
          )
        }
      }

      if (plot == "pat" | plot == "all") {
        plot(d_thf[, 2:3],
          type = "b", col = (d_thf[, 1] + 1), pch = 20, lty = "dotted",
          xlim = c(0, max(d_thf[, 2])),
          main = paste("Starting thf, flag = ", tcs$flag, ", tol = ", tol, sep = ""),
          xlab = "Cross-sectional Incoherence", ylab = "Temporal Incoherence"
        )
        graphics::text(d_thf[, 2:3], labels = d_thf[, 1], cex = 0.7, pos = 4)
      }

      if (plot == "distf" | plot == "all") {
        plot(1, 1,
          type = "n", pch = 19, xlab = "iteration", ylab = "distance forecasts iteration i and i-1",
          main = paste("Starting thf, flag = ", tcs$flag, ", tol = ", tol, sep = ""),
          ylim = c(min(dist_thf[, -c(1, 2, 3)]), max(dist_thf[, -c(1, 2, 3)])),
          xlim = c(1, i),
          xaxt = "n"
        )

        graphics::axis(1, at = 1:i)
        graphics::lines(x = dist_thf[, 1], y = dist_thf[, 2], type = "b", pch = 1, col = 1)
        graphics::lines(x = dist_thf[, 1], y = dist_thf[, 3], type = "b", pch = 2, col = 2)
        graphics::lines(x = dist_thf[, 1], y = dist_thf[, 4], type = "b", pch = 3, col = 3)
        graphics::lines(x = dist_thf[, 1], y = dist_thf[, 5], type = "b", pch = 4, col = 4)

        graphics::legend("topright",
          inset = 0.01, legend = c("norm2", "norm1", "rel1", "rel2"),
          col = c(1:4), lty = 1, pch = c(1:4), bty = "n"
        )
      }
    }
    out <- tcs
  } else if (start_rec == "hts") {
    cst <- list()
    cst$recf <- Y2_hts
    rownames(cst$recf) <- if (is.null(rownames(basef))) paste("serie", 1:NROW(cst$recf), sep = "") else rownames(basef)
    colnames(cst$recf) <- paste("k", rep(kset, h * rev(kset)), "h",
      do.call("c", as.list(sapply(
        rev(kset) * h,
        function(x) seq(1:x)
      ))),
      sep = ""
    )
    cst$d_cs <- d_cs_hts[!is.na(d_cs_hts)]
    cst$d_te <- d_te_hts[!is.na(d_te_hts)]
    cst$start_rec <- start_rec
    cst$iter <- j
    cst$tol <- tol
    cst$flag <- flag_hts
    cst$diff_iter <- dist_hts
    cst$time <- as.numeric(end_time_hts - start_time, units = "secs")

    if (note == TRUE) {
      if (plot == "mti" | plot == "all") {
        plot(0, 0,
          type = "n", pch = 19, xlab = "iteration", ylab = "Incoherence",
          main = paste("Starting hts, flag = ", cst$flag, ", tol = ", tol, sep = ""),
          ylim = c(0, max(cst$d_cs, cst$d_te)),
          xlim = c(0, j),
          xaxt = "n"
        )
        graphics::axis(1, at = 0:j)
        graphics::lines(x = 0:(length(cst$d_cs) - 1), y = cst$d_cs, type = "b", pch = 1, col = "#fe5f55")
        graphics::lines(x = 0:(length(cst$d_te) - 1), y = cst$d_te, type = "b", pch = 4, col = 1)
        if (cst$flag == 0) {
          graphics::legend("topright",
            inset = 0.01, legend = c("Cross-sectional", "Temporal"),
            col = c("#fe5f55", 1), lty = 1, pch = c(1, 4), bty = "n"
          )
        }
        if (cst$flag == 1) {
          graphics::legend("bottomright",
            inset = 0.01, legend = c("Cross-sectional", "Temporal"),
            col = c("#fe5f55", 1), lty = 1, pch = c(1, 4), bty = "n"
          )
        }
      }

      if (plot == "pat" | plot == "all") {
        plot(d_hts[, 2:3],
          type = "b", col = (d_hts[, 1] + 1), pch = 20, lty = "dotted",
          xlim = c(0, max(d_hts[, 2])),
          main = paste("Starting hts, flag = ", cst$flag, ", tol = ", tol, sep = ""),
          xlab = "Cross-sectional Incoherence", ylab = "Temporal Incoherence"
        )
        graphics::text(d_hts[, 2:3], labels = d_hts[, 1], cex = 0.7, pos = 4)
      }

      if (plot == "distf" | plot == "all") {
        plot(1, 1,
          type = "n", pch = 19, xlab = "iteration", ylab = "distance forecasts iteration i and i-1",
          main = paste("Starting hts, flag = ", cst$flag, ", tol = ", tol, sep = ""),
          ylim = c(min(dist_hts[, -c(1, 2, 3)]), max(dist_hts[, -c(1, 2, 3)])),
          xlim = c(1, j),
          xaxt = "n"
        )

        graphics::axis(1, at = 1:j)
        graphics::lines(x = dist_hts[, 1], y = dist_hts[, 2], type = "b", pch = 1, col = 1)
        graphics::lines(x = dist_hts[, 1], y = dist_hts[, 3], type = "b", pch = 2, col = 2)
        graphics::lines(x = dist_hts[, 1], y = dist_hts[, 4], type = "b", pch = 3, col = 3)
        graphics::lines(x = dist_hts[, 1], y = dist_hts[, 5], type = "b", pch = 4, col = 4)

        graphics::legend("topright",
          inset = 0.01, legend = c("norm2", "norm1", "rel1", "rel2"),
          col = c(1:4), lty = 1, pch = c(1:4), bty = "n"
        )
      }
    }
    out <- cst
  } else {
    norm2 <- c(
      sqrt(sum((Y2_hts - Y2_thf)^2)),
      sqrt(sum((Y2_hts - basef)^2)),
      sqrt(sum((Y2_thf - basef)^2))
    )
    norm1 <- c(
      sum(abs(Y2_hts - Y2_thf)),
      sum(abs(Y2_hts - basef)),
      sum(abs(Y2_thf - basef))
    )

    tcs <- list()
    tcs$recf <- Y2_thf
    rownames(tcs$recf) <- if (is.null(rownames(basef))) paste("serie", 1:NROW(tcs$recf), sep = "") else rownames(basef)
    colnames(tcs$recf) <- paste("k", rep(kset, h * rev(kset)), "h",
      do.call("c", as.list(sapply(
        rev(kset) * h,
        function(x) seq(1:x)
      ))),
      sep = ""
    )
    tcs$d_cs <- d_cs_thf[!is.na(d_cs_thf)]
    tcs$d_te <- d_te_thf[!is.na(d_te_thf)]
    tcs$start_rec <- "thf"
    tcs$iter <- i
    tcs$tol <- tol
    tcs$flag <- flag_thf
    tcs$diff_iter <- dist_thf
    tcs$time <- as.numeric(end_time_thf - start_time, units = "secs")

    cst <- list()
    cst$recf <- Y2_hts
    rownames(cst$recf) <- if (is.null(rownames(basef))) paste("serie", 1:NROW(cst$recf), sep = "") else rownames(basef)
    colnames(cst$recf) <- paste("k", rep(kset, h * rev(kset)), "h",
      do.call("c", as.list(sapply(
        rev(kset) * h,
        function(x) seq(1:x)
      ))),
      sep = ""
    )

    cst$d_cs <- d_cs_hts[!is.na(d_cs_hts)]
    cst$d_te <- d_te_hts[!is.na(d_te_hts)]
    cst$start_rec <- "hts"
    cst$iter <- j
    cst$tol <- tol
    cst$flag <- flag_hts
    cst$diff_iter <- dist_hts
    cst$time <- as.numeric(end_time_hts - end_time_thf, units = "secs")

    if (note == TRUE) {
      opar <- graphics::par(no.readonly =TRUE)
      on.exit(graphics::par(opar))

      if (plot == "mti" | plot == "all") {
        graphics::par(mfrow = c(1, 2))
        plot(0, 0,
          type = "n", pch = 19, xlab = "iteration", ylab = "Incoherence",
          main = paste("Starting hts, flag = ", cst$flag, ", tol = ", tol, sep = ""),
          ylim = c(0, max(cst$d_cs, cst$d_te)),
          xlim = c(0, j),
          xaxt = "n"
        )
        graphics::axis(1, at = 0:j)
        graphics::lines(x = 0:(length(cst$d_cs) - 1), y = cst$d_cs, type = "b", pch = 1, col = "#fe5f55")
        graphics::lines(x = 0:(length(cst$d_te) - 1), y = cst$d_te, type = "b", pch = 4, col = 1)
        if (cst$flag == 0) {
          graphics::legend("topright",
            inset = 0.01, legend = c("Cross-sectional", "Temporal"),
            col = c("#fe5f55", 1), lty = 1, pch = c(1, 4), bty = "n"
          )
        }
        if (cst$flag == 1) {
          graphics::legend("bottomright",
            inset = 0.01, legend = c("Cross-sectional", "Temporal"),
            col = c("#fe5f55", 1), lty = 1, pch = c(1, 4), bty = "n"
          )
        }

        plot(0, 0,
          type = "n", pch = 19, xlab = "iteration", ylab = "Incoherence",
          main = paste("Starting thf, flag = ", tcs$flag, ", tol = ", tol, sep = ""),
          ylim = c(0, max(tcs$d_cs, tcs$d_te)),
          xlim = c(0, i),
          xaxt = "n"
        )
        graphics::axis(1, at = 0:i)
        graphics::lines(x = 0:(length(tcs$d_cs) - 1), y = tcs$d_cs, type = "b", pch = 1, col = "#fe5f55")
        graphics::lines(x = 0:(length(tcs$d_te) - 1), y = tcs$d_te, type = "b", pch = 4, col = 1)
        if (tcs$flag == 0) {
          graphics::legend("topright",
            inset = 0.01, legend = c("Cross-sectional", "Temporal"),
            col = c("#fe5f55", 1), lty = 1, pch = c(1, 4), bty = "n"
          )
        }
        if (tcs$flag == 1) {
          graphics::legend("bottomright",
            inset = 0.01, legend = c("Cross-sectional", "Temporal"),
            col = c("#fe5f55", 1), lty = 1, pch = c(1, 4), bty = "n"
          )
        }
        graphics::par(mfrow = c(1, 1))
      }

      if (plot == "pat" | plot == "all") {
        graphics::par(mfrow = c(1, 2))
        plot(d_hts[, 2:3],
          type = "b", col = (d_hts[, 1] + 1), pch = 20, lty = "dotted",
          xlim = c(0, max(d_hts[, 2])),
          main = paste("Starting hts, flag = ", cst$flag, ", tol = ", tol, sep = ""),
          xlab = "Cross-sectional Incoherence", ylab = "Temporal Incoherence"
        )
        graphics::text(d_hts[, 2:3], labels = d_hts[, 1], cex = 0.7, pos = 4)

        plot(d_thf[, 2:3],
          type = "b", col = (d_thf[, 1] + 1), pch = 20, lty = "dotted",
          xlim = c(0, max(d_thf[, 2])),
          main = paste("Starting thf, flag = ", tcs$flag, ", tol = ", tol, sep = ""),
          xlab = "Cross-sectional Incoherence", ylab = "Temporal Incoherence"
        )
        graphics::text(d_thf[, 2:3], labels = d_thf[, 1], cex = 0.7, pos = 4)

        graphics::par(mfrow = c(1, 1))
      }

      if (plot == "distf" | plot == "all") {
        graphics::par(mfrow = c(1, 2))
        plot(1, 1,
          type = "n", pch = 19, xlab = "iteration", ylab = "distance forecasts iteration i and i-1",
          main = paste("Starting hts, flag = ", cst$flag, ", tol = ", tol, sep = ""),
          ylim = c(min(dist_hts[, -c(1, 2, 3)]), max(dist_hts[, -c(1, 2, 3)])),
          xlim = c(1, j),
          xaxt = "n"
        )

        graphics::axis(1, at = 1:j)
        graphics::lines(x = dist_hts[, 1], y = dist_hts[, 2], type = "b", pch = 1, col = 1)
        graphics::lines(x = dist_hts[, 1], y = dist_hts[, 3], type = "b", pch = 2, col = 2)
        graphics::lines(x = dist_hts[, 1], y = dist_hts[, 4], type = "b", pch = 3, col = 3)
        graphics::lines(x = dist_hts[, 1], y = dist_hts[, 5], type = "b", pch = 4, col = 4)

        graphics::legend("topright",
          inset = 0.01, legend = c("norm2", "norm1", "rel1", "rel2"),
          col = c(1:4), lty = 1, pch = c(1:4), bty = "n"
        )

        plot(1, 1,
          type = "n", pch = 19, xlab = "iteration", ylab = "distance forecasts iteration i and i-1",
          main = paste("Starting thf, flag = ", tcs$flag, ", tol = ", tol, sep = ""),
          ylim = c(min(dist_thf[, -c(1, 2, 3)]), max(dist_thf[, -c(1, 2, 3)])),
          xlim = c(1, i),
          xaxt = "n"
        )

        graphics::axis(1, at = 1:i)
        graphics::lines(x = dist_thf[, 1], y = dist_thf[, 2], type = "b", pch = 1, col = 1)
        graphics::lines(x = dist_thf[, 1], y = dist_thf[, 3], type = "b", pch = 2, col = 2)
        graphics::lines(x = dist_thf[, 1], y = dist_thf[, 4], type = "b", pch = 3, col = 3)
        graphics::lines(x = dist_thf[, 1], y = dist_thf[, 5], type = "b", pch = 4, col = 4)

        graphics::legend("topright",
          inset = 0.01, legend = c("norm2", "norm1", "rel1", "rel2"),
          col = c(1:4), lty = 1, pch = c(1:4), bty = "n"
        )
        graphics::par(mfrow = c(1, 1))
      }
    }

    if ((norm2[2] < norm2[3] | flag_thf < 0) & flag_hts >= 0) {
      out$hts <- cst
      out$thf <- tcs
      out$best <- "hts"
    } else {
      out$thf <- tcs
      out$hts <- cst
      out$best <- "thf"
    }
    out$dist <- rbind(norm2, norm1)
    colnames(out$dist) <- c("hts-thf", "hts-base", "thf-base")
    rownames(out$dist) <- c("norm2", "norm1")
  }

  return(out)
}

# Function for the cross-sectional reconciliation step
step_hts <- function(basef, kset, h, res, hts_mod) {
  Y <- lapply(kset, function(x) basef[, rep(kset, rev(kset) * h) == x, drop = FALSE])


  if (missing(res)) {
    Y1 <- lapply(Y, function(x) t(hts_mod(basef = t(x))))
  } else {
    # Create list with lenght p, with time by time temporally reconcilked residuals matrices
    r <- NCOL(res) / sum(kset)
    E <- lapply(kset, function(x) res[, rep(kset, rev(kset) * r) == x, drop = FALSE])

    ## list of time by time cross sectional M matrix
    Y1 <- mapply(function(Y, E) t(hts_mod(basef = t(Y), res = t(E))),
      Y = Y, E = E
    )
  }
  Y1 <- do.call("cbind", Y1)
  dimnames(Y1) <- NULL
  return(Y1)
}

# Function for the temporal reconciliation step
step_thf <- function(basef, res, thf_mod) {
  if (is.null(res)) {
    Y1 <- t(apply(basef, 1, function(x) thf_mod(basef = x)))
  } else {
    Y1 <- t(mapply(function(Y, X) thf_mod(basef = Y, res = X),
      Y = split(basef, row(basef)), X = split(res, row(res))
    ))
  }
}
