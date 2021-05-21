#' Top-down forecast reconciliation for genuine hierarchical/grouped time series
#'
#' @description
#' \loadmathjax
#' Top-down forecast reconciliation for genuine hierarchical/grouped time series,
#' where the forecast of a `Total' (top-level series, expected to be positive)
#' is disaggregated according to a proportional scheme given by a vector
#' of proportions (weights).
#' Besides the fulfillment of any aggregation constraint,
#' the top-down reconciled forecasts should respect two main properties:
#' - the top-level value remains unchanged;
#' - all the bottom time series reconciled forecasts are non-negative.
#' The top-down procedure is extended to deal with both temporal and cross-temporal cases.
#' Since this is a post forecasting function, the weight vector must be given
#' in input by the user, and is not calculated automatically (see Examples).
#'
#' @param topf (\mjseqn{h \times 1}) vector of the top-level base forecast to be
#' disaggregated; \mjseqn{h} is the forecast horizon (for the lowest temporal
#' aggregation order in temporal and cross-temporal cases).
#' @param C (\mjseqn{n_a \times n_b}) cross-sectional (contemporaneous) matrix
#' mapping the \mjseqn{n_b} bottom level series into the \mjseqn{n_a} higher level ones.
#' @param m Highest available sampling frequency per seasonal cycle (max. order
#' of temporal aggregation, \mjseqn{m}), or a subset of the \mjseqn{p} factors
#' of \mjseqn{m}.
#' @param weights vector of weights to be used to disaggregate topf:
#' (\mjseqn{n_b \times h}) matrix in the cross-sectional framework;
#' (\mjseqn{m \times h}) matrix in the temporal framework;
#' (\mjseqn{n_b m \times h}) matrix in the cross-temporal framework.
#'
#' @details Fix \mjseqn{h = 1}, then
#' \mjsdeqn{\widetilde{\mathbf{y}} = \mathbf{S}\mathbf{w}\widehat{a}_1}
#' where \mjseqn{\widetilde{\mathbf{y}}} is the vector of reconciled forecasts,
#' \mjseqn{\mathbf{S}} is the summing matrix (whose pattern depends on which type
#' of reconciliation is being performed), \mjseqn{\mathbf{w}} is the vector of weights,
#' and \mjseqn{\widehat{a}_1} is the top-level value to be disaggregated.
#'
#' @return The function returns an (\mjseqn{h \times n}) matrix of
#' cross-sectionally reconciled forecasts, or an (\mjseqn{h(k^\ast + m) \times 1})
#' vector of top-down temporally reconciled forecasts, or an
#' (\mjseqn{n \times h (k^\ast + m)}) matrix of top-down
#' cross-temporally reconciled forecasts.
#'
#' @references
#' Athanasopoulos, G., Ahmed, R.A., Hyndman, R.J. (2009), Hierarchical
#' forecasts for Australian domestic tourism, \emph{International Journal of
#' Forecasting}, 25, 1, 146–166.
#'
#' @keywords top-down
#' @family reconciliation procedures
#'
#' @examples
#' data(FoReco_data)
#' ### CROSS-SECTIONAL TOP-DOWN RECONCILIATION
#' # Cross sectional aggregation matrix
#' C <- FoReco_data$C
#' # monthly base forecasts
#' id <- which(simplify2array(strsplit(colnames(FoReco_data$base), split = "_"))[1, ] == "k1")
#' mbase <- t(FoReco_data$base[, id])
#' obs_1 <- FoReco_data$obs$k1
#' # average historical proportions
#' props <- colMeans(obs_1[1:168,-c(1:3)]/obs_1[1:168,1])
#' cs_td <- tdrec(topf = mbase[,1], C = C, weights = props)
#'
#' ### TEMPORAL TOP-DOWN RECONCILIATION
#' # top ts base forecasts ([lowest_freq' ...  highest_freq']')
#' top_obs12 <- FoReco_data$obs$k12[1:14,1]
#' bts_obs1 <- FoReco_data$obs$k1[1:168,1]
#' # average historical proportions
#' props <- colMeans(matrix(bts_obs1, ncol = 12, byrow = TRUE)/top_obs12)
#' topbase <- FoReco_data$base[1, 1]
#' t_td <- tdrec(topf = topbase, m = 12, weights = props)
#'
#' ### CROSS-TEMPORAL TOP-DOWN RECONCILIATION
#' top_obs <- FoReco_data$obs$k12[1:14,1]
#' bts_obs <- FoReco_data$obs$k1[1:168,-c(1:3)]
#' bts_obs <- lapply(1:5, function(x) matrix(bts_obs[,x], nrow=14, byrow = TRUE))
#' bts_obs <- do.call(cbind, bts_obs)
#' # average historical proportions
#' props <- colMeans(bts_obs/top_obs)
#' ct_td <- tdrec(topf = topbase, m = 12, C = C, weights = props)
#'
#' @export
#'
#' @import Matrix
tdrec <- function(topf, C, m, weights){
  if(missing(C) & missing(m)){
    stop("Missing C or/and m. Remember:\n - put only C for a cross-sectional top-down reconciliation\n - put only m for a temporal top-down reconciliation\n - put C AND m for a cross-temporal top-down reconciliation.")
  }else if(missing(C)){
    obj_thf <- thf_tools(m = m)
    m <- obj_thf$m
    S <- obj_thf$R
    h <- length(topf)
    kset <- obj_thf$kset

    Dh <- Dmat(h = h, m = kset, n = 1)
    vnames <- paste("k", rep(kset, h * rev(kset)), "h",
                    do.call("c", as.list(sapply(
                      rev(kset) * h,
                      function(x) seq(1:x)))),
                    sep = "")
  }else if(missing(m)){
    S <- hts_tools(C = C)$S

    cnames <- if (is.null(rownames(C)) | is.null(colnames(C))) {
      paste("serie", 1:(NROW(C)+NCOL(C)), sep = "")
    } else {
      c(rownames(C), colnames(C))
    }
    rnames <- paste("h", 1:length(topf), sep="")
  }else{
    obj_ctf <- ctf_tools(C = C, m = m, sparse = TRUE)
    nb <- obj_ctf$hts$nb
    na <- obj_ctf$hts$na
    kt <- obj_ctf$thf$kt
    ks <- obj_ctf$thf$ks
    m <- obj_ctf$thf$m
    h <- length(topf)
    kset <- obj_ctf$thf$kset
    S <- obj_ctf$ctf$Stilde

    Dh <- Dmat(h = h, m = kset, n = na+nb)
    rnames <- if (is.null(rownames(C)) | is.null(colnames(C))) {
      paste("serie", 1:(na+nb), sep = "")
    } else {
      c(rownames(C), colnames(C))
    }
    cnames <- paste("k", rep(kset, h * rev(kset)), "h",
                    do.call("c", as.list(sapply(
                      rev(kset) * h,
                      function(x) seq(1:x)))),
                    sep = "")
  }

  if(is.vector(weights) | NCOL(weights) == 1 | NROW(weights) == 1){
    if(NCOL(weights) == 1 | NROW(weights) == 1){
      weights <- as.vector(weights)
    }

    if(NROW(weights) != NCOL(S)){
      stop(paste0("The weights vector must have ", NCOL(S)," elements"))
    }

    p <- weights/sum(weights)
    recf <- t(S %*% p %*% t(topf))
  }else{
    if(NROW(weights) != NCOL(S)){
      stop(paste0("The weights matrix must have ", NCOL(S)," rows"))
    }

    if(any(abs(colSums(weights)-1)<sqrt(.Machine$double.eps))){
      weights <- weights/colSums(weights)
    }
    p <- weights
    recf <- t(S %*% p %*% diag(topf))
  }

  if(missing(C)){
    recf <- as.vector(t(Dh) %*% as.vector(t(recf)))
    recf <- stats::setNames(recf, vnames)
  }else if(missing(m)){
    recf <- as.matrix(recf)
    colnames(recf) <- cnames
    rownames(recf) <- rnames
  }else{
    recf <- matrix(t(Dh) %*% as.vector(t(recf)), nrow = na + nb, byrow = TRUE)
    colnames(recf) <- cnames
    rownames(recf) <- rnames
  }
  return(recf)
}
