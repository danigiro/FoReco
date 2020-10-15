#' @title Measuring forecasting accuracy
#'
#' @description Function to calculate the accuracy indices of the reconciled point
#' forecasts of a cross-temporal (not only, see examples) system (more in
#' \href{https://danigiro.github.io/FoReco/articles/accuracy_indices.html}{Average relative accuracy indices}).
#'
#' @param recf list of q (forecast origins) reconciled forecasts' matrices (\code{[n x h(k* + m)]} in the
#' cross-temporal, \code{[h x n]} in cross-sectional and vectors of length \code{[h(k* + m)]} in temporal framework).
#' @param base list of q (forecast origins) base forecasts' matrices (\code{[n x h(k* + m)]} in the
#' cross-temporal, \code{[n x h]} in cross-sectional and \code{[h(k* + m) x 1]} in temporal framework).
#' @param test list of q (forecast origins) test observations' matrices (\code{[n x h(k* + m)]} in the
#' cross-temporal, \code{[n x h]} in cross-sectional and \code{[h(k* + m) x 1]} in temporal framework).
#' @param type type of accuracy measure ("\code{mse}" Mean Square Error, "\code{rmse}" Root Mean Square Error
#' or "\code{mae}" Mean Absolute Error)
#' @param m highest frequency of the forecasted time series.
#' @param nb number of bottom time series in the cross-sectional framework.
#' @param compact if TRUE return only the summary matrix.
#'
#' @return
#' It returns a summary table called \code{Avg_mat} (if \code{compact} option is \code{TRUE},
#' \emph{default}), otherwise it returns a list of four tables (more in
#' \href{https://danigiro.github.io/FoReco/articles/accuracy_indices.html}{Average relative accuracy indices}).
#'
#' @references
#' Di Fonzo, T., Girolimetto, D. (2020), Cross-Temporal Forecast Reconciliation:
#' Optimal Combination Method and Heuristic Alternatives, Department of Statistical
#' Sciences, University of Padua, \href{https://arxiv.org/abs/2006.08570}{arXiv:2006.08570}.
#'
#' @examples
#' \donttest{
#' data(FoReco_data)
#'
#' # Cross-temporal framework
#' oct_recf <- octrec(FoReco_data$base, m = 12, C = FoReco_data$C,
#'                    comb = "bdshr", res = FoReco_data$res)$recf
#' oct_score <- score_index(recf = oct_recf,
#'                          base = FoReco_data$base,
#'                          test = FoReco_data$test, m = 12, nb = 5)
#'
#' # Cross-sectional framework
#' # monthly base forecasts
#' id <- which(simplify2array(strsplit(colnames(FoReco_data$base), split = "_"))[1, ] == "k1")
#' mbase <- t(FoReco_data$base[, id])
#' # monthly test set
#' mtest <- t(FoReco_data$test[, id])
#' # monthly residuals
#' id <- which(simplify2array(strsplit(colnames(FoReco_data$res), split = "_"))[1, ] == "k1")
#' mres <- t(FoReco_data$res[, id])
#' # monthly reconciled forecasts
#' mrecf <- htsrec(mbase, C = FoReco_data$C, comb = "shr", res = mres)$recf
#' # score
#' hts_score <- score_index(recf = mrecf, base = mbase, test = mtest, nb = 5)
#'
#' # Temporal framework
#' data(FoReco_data)
#' # top ts base forecasts ([lowest_freq' ...  highest_freq']')
#' topbase <- FoReco_data$base[1, ]
#' # top ts residuals ([lowest_freq' ...  highest_freq']')
#' topres <- FoReco_data$res[1, ]
#' # top ts test ([lowest_freq' ...  highest_freq']')
#' toptest <- FoReco_data$test[1, ]
#' # top ts recf ([lowest_freq' ...  highest_freq']')
#' toprecf <- thfrec(topbase, m = 12, comb = "acov", res = topres)$recf
#' # score
#' thf_score <- score_index(recf = toprecf, base = topbase, test = toptest, m = 12)
#' }
#'
#' @keywords utilities
#'
#' @export
score_index <- function(recf, base, test, m, nb, type = "mse", compact = TRUE) {

  type <- match.arg(type, c("mse", "mae", "rmse"))

  if (missing(recf) | missing(base) | missing(test)) {
    stop("The arguments recf, base and/or test are not specified.")
  }

  if (!is.list(recf)) {
    recf <- list(recf)
  }

  if (!is.list(base)) {
    base <- list(base)
  }

  if (!is.list(test)) {
    test <- list(test)
  }

  if (missing(m)) {
    m <- 1
    recf <- lapply(recf, t)
    base <- lapply(base, t)
    test <- lapply(test, t)
  }

  if (missing(nb)) {
    nb <- 1
  }

  if (is.vector(recf[[1]])) {
    recf <- lapply(recf, t)
  }

  if (is.vector(base[[1]])) {
    base <- lapply(base, t)
  }

  if (is.vector(test[[1]])) {
    test <- lapply(test, t)
  }

  test <- simplify2array(test)
  base <- simplify2array(base)
  recf <- simplify2array(recf)

  Ebase <- base - test
  Erecf <- recf - test

  if (type == "mse" | type == "rmse") {
    Qbase <- Ebase * Ebase
    Qrecf <- Erecf * Erecf
  } else if (type == "mae") {
    Qbase <- abs(Ebase)
    Qrecf <- abs(Erecf)
  }

  # Some usefull tools
  kset <- rev(divisors(m))
  kt <- sum(kset)
  p <- length(kset)
  n <- dim(Ebase)[1]
  na <- n - nb
  h <- dim(Ebase)[2] / kt
  q <- dim(Ebase)[3]
  kpos <- rep(kset, rep(rev(kset * h))) # position of k in colums of E[,,q]
  kpos <- factor(kpos, kset, ordered = TRUE)
  if (m == 1) {
    kh <- paste("k", kpos, "h", 1:h, sep = "")
  } else {
    kh <- paste("k", kpos, "h",
      do.call("c", as.list(sapply(rev(kset) * h, function(x) seq(1:x)))),
      sep = ""
    )
  }

  IND_base <- apply(Qbase, c(1, 2), mean)
  IND_recf <- apply(Qrecf, c(1, 2), mean)

  if (type == "rmse") {
    IND_base <- sqrt(IND_base)
    IND_recf <- sqrt(IND_recf)
  }

  RelIND_ikh <- IND_recf / IND_base
  colnames(RelIND_ikh) <- kh
  logRelIND_ikh <- log(RelIND_ikh)
  logRelIND_ikh[abs(logRelIND_ikh) == Inf] <- NA


  # 24
  AvgRelIND <- exp(mean(logRelIND_ikh, na.rm = TRUE))
  # 25
  AvgRelIND_i <- exp(rowMeans(logRelIND_ikh, na.rm = TRUE))
  # 27
  AvgRelIND__kh <- exp(colMeans(logRelIND_ikh, na.rm = TRUE))
  AvgRelIND_akh <- exp(colMeans(logRelIND_ikh[1:na, , drop = FALSE], na.rm = TRUE))
  AvgRelIND_bkh <- exp(colMeans(logRelIND_ikh[-c(1:na), , drop = FALSE], na.rm = TRUE))
  Avg_k <- rbind(AvgRelIND__kh, AvgRelIND_akh, AvgRelIND_bkh)
  rownames(Avg_k) <- c("all", "uts", "bts")

  # solo k
  AvgRelIND__k <- exp(tapply(colMeans(logRelIND_ikh, na.rm = TRUE), kpos, mean, na.rm = TRUE))

  # bottom and aggregate
  AvgRelIND_ak <- exp(tapply(colMeans(logRelIND_ikh[1:na, , drop = FALSE], na.rm = TRUE),
                             kpos, mean, na.rm = TRUE))
  AvgRelIND_bk <- exp(tapply(colMeans(logRelIND_ikh[-c(1:na), , drop = FALSE], na.rm = TRUE),
                             kpos, mean, na.rm = TRUE))
  AvgRelIND_a <- exp(mean(logRelIND_ikh[1:na, , drop = FALSE], na.rm = TRUE))
  AvgRelIND_b <- exp(mean(logRelIND_ikh[-c(1:na), , drop = FALSE], na.rm = TRUE))

  Avg_matrix <- data.frame(
    all = c(AvgRelIND__k, AvgRelIND),
    uts = c(AvgRelIND_ak, AvgRelIND_a), bts = c(AvgRelIND_bk, AvgRelIND_b)
  )
  rownames(Avg_matrix) <- c(unique(kset), "all")

  # per ogni serie
  AvgRelIND_ik <- exp(t(apply(logRelIND_ikh, 1, function(z) tapply(z, kpos, mean, na.rm = TRUE))))
  all <- AvgRelIND_i

  if (NROW(AvgRelIND_ik) == 1) {
    AvgRelIND_ik <- cbind(t(AvgRelIND_ik), all)
  } else {
    AvgRelIND_ik <- cbind(AvgRelIND_ik, all)
  }

  if (compact == TRUE) {
    if (nb == 1) {
      return(Avg_matrix[, 1, drop = FALSE])
    } else if (m == 1) {
      return(Avg_matrix[2, , drop = FALSE])
    } else {
      return(Avg_matrix)
    }
  } else {
    if (nb == 1) {
      out <- list()
      out$Avg_mat <- Avg_matrix[, 1, drop = FALSE]
      out$Avg_k <- Avg_k[1, , drop = FALSE]
      return(out)
    } else if (m == 1) {
      out <- list()
      out$Avg_mat <- Avg_matrix[2, , drop = FALSE]
      out$Avg_ik <- AvgRelIND_ik[, 2, drop = FALSE]
      out$Rel_mat <- RelIND_ikh
      out$Avg_k <- Avg_k
      return(out)
    } else {
      out <- list()
      out$Avg_mat <- Avg_matrix
      out$Avg_ik <- AvgRelIND_ik
      out$Rel_mat <- RelIND_ikh
      out$Avg_k <- Avg_k
      return(out)
    }
  }
}
