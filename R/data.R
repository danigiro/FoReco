#' Italian Quarterly National Accounts
#'
#' @description
#' A subset of the data used by Girolimetto et al. (2023) from the Italian Quarterly
#' National Accounts (output, income and expenditure sides) spanning the period
#' 2000:Q1-2019:Q4.
#'
#' @usage
#' # 21 time series of the Italian Quarterly National Accounts
#' itagdp
#'
#' @format
#' \code{itagdp} is a \eqn{(80 \times 21)} \code{ts} object, corresponding to
#' 21 time series of the Italian Quarterly National Accounts (2000:Q1-2019:Q4).
#'
#' @references
#' Girolimetto, D. and Di Fonzo, T. (2023), Point and probabilistic forecast reconciliation
#' for general linearly constrained multiple time series,
#' \emph{Statistical Methods & Applications}, 33, 581-607. \doi{10.1007/s10260-023-00738-6}.
#'
#' @source <https://ec.europa.eu/eurostat/web/national-accounts/>
"itagdp"

#' @rdname itagdp
#'
#' @usage
#' # 'agg_mat' and 'cons_mat' for the output side
#' outside
#'
#' @format \code{outside}, \code{income} and \code{expenditure} are lists with two elements:
#' \itemize{
#' \item{\code{agg_mat} contains the \eqn{(1 \times 2)}, \eqn{(2 \times 4)}, or \eqn{(6 \times 8)}
#' aggregation matrix according to output, income or expenditure side, respectively.}
#' \item{\code{cons_mat} contains the \eqn{(1 \times 3)}, \eqn{(2 \times 6)}, or \eqn{(6 \times 14)}
#' zero constraints matrix according to output, income or expenditure side, respectively.}
#' }
#'
"outside"

#' @rdname itagdp
#'
#' @usage
#' # 'agg_mat' and 'cons_mat' for the expenditure side
#' expside
#'
#' @format NULL
"expside"

#' @rdname itagdp
#'
#' @usage
#' # 'agg_mat' and 'cons_mat' for the income side
#' incside
#'
#' @format NULL
"incside"

#' @rdname itagdp
#'
#' @usage
#' # zero constraints matrix encompassing output, expenditure and income sides
#' gdpconsmat
#'
#' @format \code{gdpconsmat} is the complete \eqn{(9 \times 21)} zero constraints matrix
#' encompassing output, expenditure and income sides.
"gdpconsmat"

#' Australian Tourism Demand dataset
#'
#' @description
#' The Australian Tourism Demand dataset (Wickramasuriya et al. 2019) measures the number of
#' nights Australians spent away from home. It includes 228 monthly observations of Visitor
#' Nights (VNs) from January 1998 to December 2016, and has a cross-sectional grouped
#' structure based on a geographic hierarchy crossed by purpose of travel. The geographic
#' hierarchy comprises 7 states, 27 zones, and 76 regions, for a total of 111 nested geographic
#' divisions. Six of these zones are each formed by a single region, resulting in 105 unique
#' nodes in the hierarchy. The purpose of travel comprises four categories: holiday, visiting
#' friends and relatives, business, and other. To avoid redundancies (Girolimetto et al. 2023),
#' 24 nodes (6 zones are formed by a single region) are not considered, resulting in an unbalanced
#' hierarchy of 525 (304 bottom and 221 upper time series) unique nodes instead of the theoretical
#' 555 with duplicated nodes.
#'
#' @usage
#' # 525 time series of the Australian Tourism Demand dataset
#' vndata
#'
#' @format
#' \code{vndata} is a \eqn{(228 \times 525)} \code{ts} object, corresponding to
#' 525 time series of the Australian Tourism Demand dataset (1998:01-2016:12).
#'
#' @references
#' Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J. (2024),
#' Cross-temporal probabilistic forecast reconciliation: Methodological and
#' practical issues. \emph{International Journal of Forecasting},  40, 3, 1134-1151.
#' \doi{10.1016/j.ijforecast.2023.10.003}
#'
#' Wickramasuriya, S.L., Athanasopoulos, G. and Hyndman, R.J. (2019), Optimal forecast
#' reconciliation for hierarchical and grouped time series through trace minimization,
#' \emph{Journal of the American Statistical Association}, 114, 526, 804-819.
#' \doi{10.1080/01621459.2018.1448825}
#'
#' @source <https://robjhyndman.com/publications/mint/>
"vndata"

#' @rdname vndata
#'
#' @usage
#' # aggregation matrix
#' vnaggmat
#'
#' @format \code{vnaggmat} is the \eqn{(221 \times 304)} aggregation matrix.
"vnaggmat"
