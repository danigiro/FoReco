#' Reconciled forecasts matrix/vector to time-series class
#'
#' @description
#' \loadmathjax
#' Function to transform the matrix/vector of FoReco forecasts input and output into a list of
#' time series/matrix/vector objects.
#'
#' @param recf (\mjseqn{h(k^\ast + m) \times 1}) forecasts vector from \link[FoReco]{thfrec},
#' (\mjseqn{h \times n}) forecasts matrix from \link[FoReco]{htsrec} or
#' (\mjseqn{n \times h(k^\ast + m)}) forecasts matrix from \link[FoReco]{octrec},
#' \link[FoReco]{tcsrec}, \link[FoReco]{cstrec}, \link[FoReco]{iterec},
#' \link[FoReco]{ctbu}.
#' @param m Highest available sampling frequency per seasonal cycle (max. order
#' of temporal aggregation, \mjseqn{m}), or a subset of \mjseqn{p} factors
#' of \mjseqn{m}.
#' @param ... optional arguments to \link[stats]{ts} (i.e. starting date);
#' frequency is required only for the cross-sectional case.
#'
#' @return
#' A list of class \code{"ts"} objects for \link[FoReco]{FoReco2ts} and a list of
#' matrix/vector for \link[FoReco]{FoReco2matrix}
#'
#' @examples
#' data(FoReco_data)
#' # Cross-temporal framework
#' oct_recf <- octrec(FoReco_data$base, m = 12, C = FoReco_data$C,
#'                    comb = "bdshr", res = FoReco_data$res)$recf
#' ts_oct <- FoReco2ts(recf = oct_recf, m = 12, start = c(15, 1))
#' mat_oct <- FoReco2matrix(recf = oct_recf, m = 12)
#'
#' # Cross-sectional framework
#' # monthly base forecasts
#' mbase <- FoReco2matrix(FoReco_data$base, m = 12)$k1
#' # monthly residuals
#' mres <- FoReco2matrix(FoReco_data$res, m = 12)$k1
#' hts_recf <- htsrec(mbase, C = FoReco_data$C, comb = "shr", res = mres)$recf
#' ts_hts <- FoReco2ts(recf = hts_recf, start = c(15, 1), frequency = 12)
#' mat_hts <- FoReco2matrix(recf = hts_recf)
#'
#' # Temporal framework
#' # top ts base forecasts ([lowest_freq' ...  highest_freq']')
#' topbase <- FoReco_data$base[1, ]
#' # top ts residuals ([lowest_freq' ...  highest_freq']')
#' topres <- FoReco_data$res[1, ]
#' thf_recf <- thfrec(topbase, m = 12, comb = "acov", res = topres)$recf
#' ts_thf <- FoReco2ts(recf = thf_recf, m = 12, start = c(15, 1))
#' mat_thf <- FoReco2matrix(recf = thf_recf, m = 12)
#'
#' @keywords utilities
#' @family utilities
#'
#' @export
FoReco2ts <- function(recf, m, ...) {
  if(missing(m)){
    return(stats::ts(recf, ...))
  }else{
    info <- thf_tools(m = m)
    m <- info$m
    kset <- info$kset
    if(is.vector(recf)) {
      h <- length(recf)/info$kt
      nm <- names(recf) <-  paste0("k", rep(kset, h*m/kset),
                                      "h", unlist(sapply(m/kset, function(x) 1:(x*h))))
      numkh <- strsplit(gsub("k", "", nm), "h")
      indkh <- data.frame(apply(t(simplify2array(numkh)), 2, as.numeric))
      colnames(indkh) <- c("k", "h")
      ts_vec <- function(x) stats::ts(x, frequency = length(x) / indkh$h[1], ...)
      names(recf) <- NULL
      out <- tapply(recf, indkh$k, ts_vec, simplify = FALSE)
      names(out) <- paste("k", names(out), sep = "")
      return(rev(out))
    }else if (is.matrix(recf)) {
      h <- NCOL(recf)/info$kt
      nm <- colnames(recf) <-  paste0("k", rep(kset, h*m/kset),
                                      "h", unlist(sapply(m/kset, function(x) 1:(x*h))))
      numkh <- strsplit(gsub("k", "", nm), "h")
      indkh <- data.frame(apply(t(simplify2array(numkh)), 2, as.numeric))
      colnames(indkh) <- c("k", "h")
      recf_t <- t(recf)
      out <- lapply(unique(indkh$k), function(x) stats::ts(recf_t[indkh$k == x, , drop = FALSE], frequency = sum(indkh$k == x) / indkh$h[1], ...))
      names(out) <- paste("k", unique(indkh$k), sep = "")
      return(out)
    }else{
      stop("recf must be a vector or a matrix", call. = FALSE)
    }
  }
}

#' @rdname FoReco2ts
#' @importFrom stats tsp<-
#' @export
FoReco2matrix <- function(recf, m){
  if(missing(m)){
    obj <- FoReco2ts(recf = recf)
  }else{
    obj <- FoReco2ts(recf = recf, m = m)
  }
  if(is.list(obj)){
    lapply(obj, function(x) `tsp<-`(x, NULL))
  }else{
    `tsp<-`(obj, NULL)
  }
}
