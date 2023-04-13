#' @title Measuring accuracy in a rolling forecast experiment
#'
#' @description
#' \loadmathjax
#' Function to calculate the accuracy indices of the reconciled point
#' forecasts of a cross-temporal (not only, see examples) system (more in
#' \href{https://danigiro.github.io/FoReco/articles/accuracy_indices.html}{Average relative accuracy indices}).
#' (\emph{Experimental version})
#'
#' @param recf list of q (forecast origins) reconciled forecasts' matrices
#' (\mjseqn{[n \times h(k^\ast + m)]} in the cross-temporal case,
#' \mjseqn{[h \times n]} in the cross-sectional case, and vectors of length
#' \mjseqn{[h(k^\ast \times m)]} in the temporal framework).
#' @param base list of q (forecast origins) base forecasts' matrices
#' (\mjseqn{[n \times h(k^\ast + m)]} in the cross-temporal case,
#' \mjseqn{[h \times n]} in the cross-sectional case, and vectors of length
#' \mjseqn{[h(k^\ast \times m)]} in the temporal framework).
#' @param test list of q (forecast origins) test observations' matrices
#' (\mjseqn{[n \times h(k^\ast + m)]} in the cross-temporal case,
#' \mjseqn{[h \times n]} in the cross-sectional case, and vectors of length
#' \mjseqn{[h(k^\ast \times m)]} in the temporal framework).
#' @param type type of accuracy measure ("\code{mse}" Mean Square Error,
#' "\code{rmse}" Root Mean Square Error or "\code{mae}" Mean Absolute Error).
#' @param m Highest available sampling frequency per seasonal cycle (max. order
#' of temporal aggregation, \mjseqn{m}), or a subset of \mjseqn{p} factors
#' of \mjseqn{m}.
#' @param nb number of bottom time series in the cross-sectional framework.
#' @param compact if TRUE returns only the summary matrix.
#' @param nl (\mjseqn{L \times 1}) vector containing the number of time series
#' in each cross-sectional level of the hierarchy (\code{nl[1] = 1}).
#'
#' @return
#' It returns a summary table called \code{Avg_mat} (if \code{compact} option is \code{TRUE},
#' \emph{default}), otherwise it returns a list of six tables (more in
#' \href{https://danigiro.github.io/FoReco/articles/accuracy_indices.html}{Average relative accuracy indices}).
#'
#' @references
#' Di Fonzo, T., and Girolimetto, D. (2023), Cross-temporal forecast reconciliation:
#' Optimal combination method and heuristic alternatives, \emph{International Journal
#' of Forecasting}, 39(1), 39-57.
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
#' mbase <- FoReco2matrix(FoReco_data$base, m = 12)$k1
#' # monthly test set
#' mtest <- FoReco2matrix(FoReco_data$test, m = 12)$k1
#' # monthly residuals
#' mres <- FoReco2matrix(FoReco_data$res, m = 12)$k1
#' # monthly reconciled forecasts
#' mrecf <- htsrec(mbase, C = FoReco_data$C, comb = "shr", res = mres)$recf
#' # score
#' hts_score <- score_index(recf = mrecf, base = mbase, test = mtest, nb = 5)
#'
#' # Temporal framework
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
#' @family utilities
#'
#' @export
score_index <- function(recf, base, test, m, nb, nl, type = "mse", compact = TRUE) {

  type <- match.arg(type, c("mse", "mae", "rmse"))

  if (missing(recf) | missing(base) | missing(test)) {
    stop("The arguments recf, base and/or test are not specified.", call. = FALSE)
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
  if(length(m)==1){
    if(m == 1){
      kset <- 1
    }else{
      thf_obj <- thf_tools(m = m)
      m <- thf_obj$m
      kset <- thf_obj$kset
    }
  }else{
    thf_obj <- thf_tools(m = m)
    m <- thf_obj$m
    kset <- thf_obj$kset
  }
  kt <- sum(kset)
  p <- length(kset)
  n <- dim(Ebase)[1]
  na <- n - nb
  h <- dim(Ebase)[2] / sum(m/kset)
  q <- dim(Ebase)[3]
  kpos <- rep(kset, rep(m/kset * h)) # position of k in colums of E[,,q]
  kpos <- factor(kpos, kset, ordered = TRUE)
  if (m == 1) {
    kh <- paste("h", 1:h, sep = "")
    khc <- paste("h1:", 1:h, sep = "")
  } else {
    kh <- paste("k", kpos, "h",
                do.call("c", lapply(m/kset * h, function(x) seq(1:x))),
                sep = ""
    )

    khc <- paste("k", kpos, "h1:",
                 do.call("c", lapply(m/kset * h, function(x) seq(1:x))),
                 sep = ""
    )
  }

  if(missing(nl)){
    nl <- na
  }else if(sum(nl) != na){
    stop("Please, provide a valid nl vector s.t. sum(nl) == na", call. = FALSE)
  }

  nlnb <- c(nl, nb)
  L <- length(nlnb)
  levpos <- rep(1:L, nlnb)

  if(L>2){
    csname <- paste0("L", 1:L)
  }else{
    csname <- c("uts", "bts")
  }

  IND_base <- apply(Qbase, c(1, 2), mean, na.rm = TRUE)
  IND_recf <- apply(Qrecf, c(1, 2), mean, na.rm = TRUE)

  if (type == "rmse") {
    IND_base <- sqrt(IND_base)
    IND_recf <- sqrt(IND_recf)
  }

  RelIND_ikh <- IND_recf / IND_base
  colnames(RelIND_ikh) <- kh
  RelIND_ikh[is.nan(RelIND_ikh)] <- 1
  logRelIND_ikh <- log(RelIND_ikh)
  #logRelIND_ikh[abs(logRelIND_ikh) == Inf] <- NA

  logRelIND_ikh[logRelIND_ikh == -Inf] <- log(1e-16)
  if(any(logRelIND_ikh == Inf)){
    if(all(logRelIND_ikh == Inf)){
      logRelIND_ikh[abs(logRelIND_ikh) == Inf] <- NA
      warning("Base forecasts perfect prediction",call. = FALSE)
    }else{
      logRelIND_ikh[logRelIND_ikh == Inf] <- 10*max(logRelIND_ikh[logRelIND_ikh != Inf])
    }
  }

  # 24
  AvgRelIND <- exp(mean(logRelIND_ikh, na.rm = TRUE))
  # 25
  AvgRelIND_i <- exp(rowMeans(logRelIND_ikh, na.rm = TRUE))
  # 27
  AvgRelIND__kh <- exp(colMeans(logRelIND_ikh, na.rm = TRUE))
  AvgRelIND_lkh <- lapply(1:L, function(x) exp(colMeans(logRelIND_ikh[levpos == x, , drop = FALSE], na.rm = TRUE)))
  AvgRelIND_lkh <- do.call(rbind, AvgRelIND_lkh)
  Avg_k <- rbind(AvgRelIND__kh, AvgRelIND_lkh)
  rownames(Avg_k) <- c("all", csname)

  # only k
  AvgRelIND__k <- exp(tapply(colMeans(logRelIND_ikh, na.rm = TRUE), kpos, mean, na.rm = TRUE))

  # bottom and aggregate
  AvgRelIND_lk <- lapply(1:L, function(x) exp(tapply(colMeans(logRelIND_ikh[levpos == x, , drop = FALSE], na.rm = TRUE),
                                                     kpos, mean, na.rm = TRUE)))
  AvgRelIND_l <- lapply(1:L, function(x) exp(mean(logRelIND_ikh[levpos == x, , drop = FALSE], na.rm = TRUE)))
  Avg_matrix <- cbind(c(AvgRelIND__k, AvgRelIND),
                      as.data.frame(t(cbind(do.call(rbind, AvgRelIND_lk),
                                            do.call(rbind, AvgRelIND_l)))))
  colnames(Avg_matrix) <- c("all", csname)
  rownames(Avg_matrix) <- c(unique(kset), "all")

  # for each series
  AvgRelIND_ik <- exp(t(apply(logRelIND_ikh, 1, function(z) tapply(z, kpos, mean, na.rm = TRUE))))
  all <- AvgRelIND_i

  if (NROW(AvgRelIND_ik) == 1) {
    AvgRelIND_ik <- cbind(t(AvgRelIND_ik), all)
  } else {
    AvgRelIND_ik <- cbind(AvgRelIND_ik, all)
  }

  # CUMULATIVE forecast horizon
  AvgRelIND__khc <- exp(do.call(c, lapply(unique(kpos), function(y){
    cummean(colMeans(logRelIND_ikh[, kpos == y, drop = FALSE], na.rm = TRUE))
  })))
  AvgRelIND_lkhc <- lapply(1:L, function(x) exp(do.call(c, lapply(unique(kpos), function(y){
    cummean(colMeans(logRelIND_ikh[levpos == x, kpos == y, drop = FALSE], na.rm = TRUE))
  }))))
  AvgRelIND_lkhc <- do.call(rbind, AvgRelIND_lkhc)
  Avg_kc <- rbind(AvgRelIND__khc, AvgRelIND_lkhc)
  rownames(Avg_kc) <- c("all", csname)
  colnames(Avg_kc) <- khc

  RelIND_ikhc <- t(apply(logRelIND_ikh, 1,
                         function(x)
                           exp(do.call(c, lapply(unique(kpos),
                                                 function(y){
                                                   cummean(x[kpos == y])
                                                 })))))
  if(NCOL(RelIND_ikhc) != length(khc)){
    RelIND_ikhc <- t(RelIND_ikhc)
    colnames(RelIND_ikhc) <- khc
  }else{
    colnames(RelIND_ikhc) <- khc
  }

  if (compact == TRUE) {
    if (nb == 1) {
      return(Avg_matrix[, "all", drop = FALSE])
    } else if (m == 1) {
      return(Avg_matrix["all", , drop = FALSE])
    } else {
      return(Avg_matrix)
    }
  } else {
    if (nb == 1) {
      out <- list()
      out$Avg_mat <- Avg_matrix[, "all", drop = FALSE]
      out$Avg_k <- Avg_k["all", , drop = FALSE]
      out$Avg_k_cum <- Avg_kc["all", , drop = FALSE]
      return(out)
    } else if (m == 1) {
      out <- list()
      out$Avg_mat <- Avg_matrix["all", , drop = FALSE]
      out$Avg_ik <- AvgRelIND_ik[, "all", drop = FALSE]
      out$Rel_mat <- RelIND_ikh
      out$Avg_k <- Avg_k
      out$Rel_mat_cum <- RelIND_ikhc
      out$Avg_k_cum <- Avg_kc
      return(out)
    } else {
      out <- list()
      out$Avg_mat <- Avg_matrix
      out$Avg_ik <- AvgRelIND_ik
      out$Rel_mat <- RelIND_ikh
      out$Avg_k <- Avg_k
      out$Rel_mat_cum <- RelIND_ikhc
      out$Avg_k_cum <- Avg_kc
      return(out)
    }
  }
}

cummean <- function(x){
  #x[is.na(x)] <- 0
  cumsum(x) / seq_along(x)
}
