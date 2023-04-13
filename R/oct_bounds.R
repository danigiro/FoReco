#' Optimal cross-temporal bounds
#'
#' @description
#' \loadmathjax
#' Function to export the constraints designed for the cross-sectional and/or
#' temporal reconciled forecasts
#'
#' @param hts_bounds (\mjseqn{n \times 2}) matrix with cross-sectional bounds:
#' the first column is the lower bound, and the second column is the upper bound.
#' @param thf_bounds (\mjseqn{(k^\ast + m) \times 2}) matrix with temporal bounds:
#' the first column is the lower bound, and the second column is the upper bound.
#' @param m Highest available sampling frequency per seasonal cycle (max. order
#' of temporal aggregation, \mjseqn{m}), or a subset of \mjseqn{p} factors
#' of \mjseqn{m}.
#' @param C (\mjseqn{n_a \times n_b}) cross-sectional (contemporaneous) matrix
#' mapping the bottom level series into the higher level ones.
#' @param Ut Zero constraints cross-sectional (contemporaneous) kernel matrix
#' \mjseqn{(\mathbf{U}'\mathbf{y} = \mathbf{0})} spanning the null space valid
#' for the reconciled forecasts. It can be used instead of parameter
#' \code{C}, but \code{nb} (\mjseqn{n = n_a + n_b}) is needed if
#' \mjseqn{\mathbf{U}' \neq [\mathbf{I} \ -\mathbf{C}]}{}. If the hierarchy
#' admits a structural representation, \mjseqn{\mathbf{U}'} has dimension
#' (\mjseqn{n_a \times n}).
#'
#' @return A matrix with the cross-temporal bounds.
#'
#' @keywords utilities
#' @family utilities
#'
#' @examples
#' data(FoReco_data)
#' # monthly base forecasts
#' mbase <- FoReco2matrix(FoReco_data$base, m = 12)$k1
#' # monthly residuals
#' mres <- FoReco2matrix(FoReco_data$res, m = 12)$k1
#'
#' # For example, in FoReco_data we want that BA > 78, and C > 50
#' cs_bound <- matrix(c(rep(-Inf, 5), 78, -Inf, 50, rep(+Inf, 8)), ncol = 2)
#' ## Cross-sectional reconciliation
#' csobj <- htsrec(mbase, C = FoReco_data$C, comb = "shr", res = mres, bounds = cs_bound)
#'
#' # Extension of the constraints to the cross-temporal case
#' ct_bound <- oct_bounds(hts_bounds = cs_bound, m = 12)
#' ## Cross-temporal reconciliation
#' obj <- octrec(FoReco_data$base, m = 12, C = FoReco_data$C, comb = "bdshr",
#'               res = FoReco_data$res, bounds = ct_bound)
#'
#' @export
oct_bounds <- function(hts_bounds, thf_bounds, m, C, Ut){
  if(missing(thf_bounds) & missing(hts_bounds)){
    stop("ciao")
  }else if(missing(thf_bounds) | missing(hts_bounds)){
    print("c")
    if(missing(thf_bounds)){
      hts_nrow <- NROW(hts_bounds)
      thf_nrow <- thf_tools(m = m)$kt

      if(!missing(C)){
        n <- NROW(C)+NCOL(C)
        if(hts_nrow != n){
          stop("hts_bounds must be a matrix (", n, " x 2)", call. = FALSE)
        }
      }else if(!missing(Ut)){
        n <- NCOL(Ut)
        if(hts_nrow != n){
          stop("hts_bounds must be a matrix (", n, " x 2)", call. = FALSE)
        }
      }

      bhts <- t(simplify2array(rep(split(hts_bounds, 1:hts_nrow), each = thf_nrow)))
      dimnames(bhts) <- NULL
      colnames(bhts) <- c("lower", "upper")
      return(bhts)
    }else if(missing(hts_bounds)){
      if(!missing(C)){
        hts_nrow <- NROW(C)+NCOL(C)
      }else if(!missing(Ut)){
        hts_nrow <- NCOL(Ut)
      }else{
        stop("the argument C or Ut is not specified", call. = FALSE)
      }

      thf_nrow <- NROW(thf_bounds)

      if(!missing(m)){
        kt <- thf_tools(m = m)$kt
        if(thf_nrow != kt){
          stop("thf_bounds must be a matrix (", kt, " x 2)", call. = FALSE)
        }
      }

      bthf <- t(simplify2array(rep(split(thf_bounds, 1:thf_nrow), hts_nrow)))
      dimnames(bthf) <- NULL
      colnames(bthf) <- c("lower", "upper")
      return(bthf)
    }
  }else{
    kt <- thf_tools(m = m)$kt

    if(!missing(C)){
      n <- NROW(C)+NCOL(C)
    }else if(!missing(Ut)){
      n <- NCOL(Ut)
    }else{
      stop("the argument C or Ut is not specified", call. = FALSE)
    }
    thf_nrow <- NROW(thf_bounds)
    hts_nrow <- NROW(hts_bounds)

    if(hts_nrow != n){
      stop("hts_bounds must be a matrix (", n, " x 2)", call. = FALSE)
    }

    if(thf_nrow != kt){
      stop("thf_bounds must be a matrix (", kt, " x 2)", call. = FALSE)
    }

    bthf <- t(simplify2array(rep(split(thf_bounds, 1:thf_nrow), hts_nrow)))
    bhts <- t(simplify2array(rep(split(hts_bounds, 1:hts_nrow), each = thf_nrow)))
    out <- cbind(apply(cbind(bhts[,1], bthf[,1]), 1, max),
                 apply(cbind(bhts[,2], bthf[,2]), 1, min))
    dimnames(out) <- NULL
    colnames(out) <- c("lower", "upper")
    return(out)
  }
}

