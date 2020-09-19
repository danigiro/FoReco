#' Simple examples to compare \code{FoReco} and \code{hts} packages
#'
#' Two datasets of the \pkg{hts} package are used to show how to get the same results
#' using \pkg{FoReco}. First, we consider the \code{htseg1} datatset
#' (a simulated three level hierarchy, with a total of 8 series, each of length 10).
#' Then, we take the \code{htseg2} datatset (a simulated four level hierarchy with
#' a total of 17 series, each of length 16). \code{htseg1} and \code{htseg2}
#' are objects of class \code{hts} in \pkg{hts}.
#'
#' @name FoReco-hts
#' @docType data
#' @references
#' Hyndman, R. J., Lee, A., Wang, E., and Wickramasuriya, S. (2020).
#' hts: Hierarchical and Grouped Time Series, \emph{R package version 6.0.1},
#' \href{https://CRAN.R-project.org/package=hts}{https://CRAN.R-project.org/package=hts}.
#'
#' @keywords datasets
#'
#' @examples
#' \dontrun{
#' library(hts)
#' require(FoReco)
#'
#' ####### htseg1 #######
#' data <- allts(htseg1)
#' n <- NCOL(data)
#' nb <- NCOL(htseg1$bts)
#' na <- n-nb
#' C <- smatrix(htseg1)[1:na, ]
#'
#' # List containing the base forecasts
#' # Forecast horizon: 10
#' base <- list()
#' for (i in 1:n) {
#'   base[[i]] <- forecast(auto.arima(data[, i]))
#' }
#'
#' # Create the matrix of base forecasts
#' BASE <- NULL
#' for (i in 1:n) {
#'   BASE <- cbind(BASE, base[[i]]$mean)
#' }
#' colnames(BASE) <- colnames(data)
#'
#' # Create the matrix of residuals
#' res <- NULL
#' for (i in 1:n) {
#'   res <- cbind(res, base[[i]]$residuals)
#' }
#' colnames(res) <- colnames(data)
#'
#' # ols
#' Y_hts_forecast <- forecast(htseg1, method = "comb", fmethod = "arima", weights = "ols")
#' Y_hts_ols <- combinef(BASE, nodes = get_nodes(htseg1), keep = "all")
#' sum((allts(Y_hts_forecast) - Y_hts_ols) > 1e-10)
#' Y_FoReco_ols <- htsrec(BASE, C = C, comb = "ols")$recf
#' sum((Y_hts_ols - Y_FoReco_ols) > 1e-10)
#'
#' # struc
#' w <- 1 / apply(smatrix(htseg1), 1, sum)
#' Y_hts_struc <- combinef(BASE, nodes = get_nodes(htseg1), weights = w, keep = "all")
#' Y_FoReco_struc <- htsrec(BASE, C = C, comb = "struc")$recf
#' sum((Y_hts_struc - Y_FoReco_struc) > 1e-10)
#'
#' # shr
#' Y_hts_shr <- MinT(BASE, nodes = get_nodes(htseg1), keep = "all", covariance = "shr", residual = res)
#' Y_FoReco_shr <- htsrec(BASE, C = C, comb = "shr", res = res)$recf
#' sum((Y_hts_shr - Y_FoReco_shr) > 1e-10)
#'
#' # sam - hts error "MinT needs covariance matrix to be positive definite"
#' #       The covariance matrix is non-singular, but its condition number is very low,
#' #       and hts considers it as non-invertible
#' Y_hts_sam <- MinT(BASE, nodes = get_nodes(htseg1), keep = "all", covariance = "sam", residual = res)
#' Y_FoReco_sam <- htsrec(BASE, C = C, comb = "sam", res = res)$recf
#' # sum((Y_hts_sam-Y_FoReco_sam)>1e-10)
#'
#' ####### htseg2 #######
#' data <- allts(htseg2)
#' n <- NCOL(data)
#' nb <- NCOL(htseg2$bts)
#' na <- n-nb
#' C <- smatrix(htseg2)[1:na, ]
#'
#' ## In FoReco, forecasts must be obtained externally
#' # List containing the base forecasts
#' # Forecast horizon: 10
#' base <- list()
#' for (i in 1:n) {
#'   base[[i]] <- forecast(auto.arima(data[, i]))
#' }
#'
#' # Create the matrix of base forecasts
#' BASE <- NULL
#' for (i in 1:n) {
#'   BASE <- cbind(BASE, base[[i]]$mean)
#' }
#' colnames(BASE) <- colnames(data)
#'
#' # Create the matrix of residuals
#' res <- NULL
#' for (i in 1:n) {
#'   res <- cbind(res, base[[i]]$residuals)
#' }
#' colnames(res) <- colnames(data)
#'
#' ## Comparison
#' # ols
#' Y_hts_forecast <- forecast(htseg2, method = "comb", fmethod = "arima", weights = "ols")
#' Y_hts_ols <- combinef(BASE, nodes = get_nodes(htseg2), keep = "all")
#' sum((allts(Y_hts_forecast) - Y_hts_ols) > 1e-10)
#' Y_FoReco_ols <- htsrec(BASE, C = C, comb = "ols")$recf
#' sum(abs(Y_hts_ols - Y_FoReco_ols) > 1e-10)
#'
#' # struc
#' w <- 1 / apply(smatrix(htseg2), 1, sum)
#' Y_hts_struc <- combinef(BASE, nodes = get_nodes(htseg2), weights = w, keep = "all")
#' Y_FoReco_struc <- htsrec(BASE, C = C, comb = "struc")$recf
#' sum(abs(Y_hts_struc - Y_FoReco_struc) > 1e-10)
#'
#' # shr
#' Y_hts_shr <- MinT(BASE, nodes = get_nodes(htseg2), keep = "all", covariance = "shr", residual = res)
#' Y_FoReco_shr <- htsrec(BASE, C = C, comb = "shr", res = res)$recf
#' sum(abs(Y_hts_shr- Y_FoReco_shr) > 1e-10)
#' }
NULL
