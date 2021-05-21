#' Reconciled forecasts matrix/vector to time-series class
#'
#' @description
#' \loadmathjax
#' Function to transform the matrix/vector of reconciled forecasts
#' (\code{recf} from \link[FoReco]{htsrec}, \link[FoReco]{thfrec}, \link[FoReco]{tdrec},
#' \link[FoReco]{octrec}, \link[FoReco]{lccrec}, \link[FoReco]{tcsrec}, \link[FoReco]{cstrec},
#' \link[FoReco]{iterec}, \link[FoReco]{ctbu}) into a list of time series objects.
#'
#' @param recf (\mjseqn{h(k^\ast + m) \times 1}) reconciled forecasts vector from \link[FoReco]{thfrec},
#' (\mjseqn{h \times n}) reconciled forecasts matrix from \link[FoReco]{htsrec} or
#' (\mjseqn{n \times h(k^\ast + m)}) reconciled forecasts matrix from \link[FoReco]{octrec},
#' \link[FoReco]{tcsrec}, \link[FoReco]{cstrec}, \link[FoReco]{iterec},
#' \link[FoReco]{ctbu}.
#' @param ... optional arguments to \link[stats]{ts} (i.e. starting date);
#' frequency is required only for the cross-sectional case.
#'
#' @return
#' A list of class \code{"ts"} objects
#'
#' @examples
#' data(FoReco_data)
#' # Cross-temporal framework
#' oct_recf <- octrec(FoReco_data$base, m = 12, C = FoReco_data$C,
#'                    comb = "bdshr", res = FoReco_data$res)$recf
#' obj_oct <- FoReco2ts(recf = oct_recf, start = c(15, 1))
#'
#' # Cross-sectional framework
#' # monthly base forecasts
#' id <- which(simplify2array(strsplit(colnames(FoReco_data$base),
#'                                     split = "_"))[1, ] == "k1")
#' mbase <- t(FoReco_data$base[, id])
#' # monthly residuals
#' id <- which(simplify2array(strsplit(colnames(FoReco_data$res),
#'                                     split = "_"))[1, ] == "k1")
#' mres <- t(FoReco_data$res[, id])
#' hts_recf <- htsrec(mbase, C = FoReco_data$C, comb = "shr", res = mres)$recf
#' obj_hts <- FoReco2ts(recf = hts_recf, start = c(15, 1), frequency = 12)
#'
#' # Temporal framework
#' # top ts base forecasts ([lowest_freq' ...  highest_freq']')
#' topbase <- FoReco_data$base[1, ]
#' # top ts residuals ([lowest_freq' ...  highest_freq']')
#' topres <- FoReco_data$res[1, ]
#' thf_recf <- thfrec(topbase, m = 12, comb = "acov", res = topres)$recf
#' obj_thf <- FoReco2ts(recf = thf_recf, start = c(15, 1))
#'
#' @keywords utilities
#' @family utilities
#'
#' @export
FoReco2ts <- function(recf, ...) {
  if (is.vector(recf)) {
    nm <- names(recf)
    numkh <- strsplit(gsub("k", "", nm), "h")
    indkh <- data.frame(apply(t(simplify2array(numkh)), 2, as.numeric))
    colnames(indkh) <- c("k", "h")
    ts_vec <- function(x) stats::ts(x, frequency = length(x) / indkh$h[1], ...)
    names(recf) <- NULL
    out <- tapply(recf, indkh$k, ts_vec, simplify = FALSE)
    names(out) <- paste("k", names(out), sep = "")
    out <- rev(out)
  } else if (is.matrix(recf)) {
    nm <- colnames(recf)
    if (length(grep("^(k)[0-9]{1,}(h)[0-9]{1,}$", nm, perl = TRUE)) < length(nm)) {
      out <- stats::ts(recf, ...)
    } else {
      numkh <- strsplit(gsub("k", "", nm), "h")
      indkh <- data.frame(apply(t(simplify2array(numkh)), 2, as.numeric))
      colnames(indkh) <- c("k", "h")
      recf_t <- t(recf)
      out <- lapply(unique(indkh$k), function(x) stats::ts(recf_t[indkh$k == x, , drop = FALSE], frequency = sum(indkh$k == x) / indkh$h[1], ...))
      names(out) <- paste("k", unique(indkh$k), sep = "")
    }
  } else {
    warning("recf must be a vector or a matrix", call. = FALSE)
  }

  return(out)
}
