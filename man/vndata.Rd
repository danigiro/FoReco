% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{vndata}
\alias{vndata}
\alias{vnaggmat}
\title{Australian Tourism Demand dataset}
\format{
\code{vndata} is a \eqn{(228 \times 525)} \code{ts} object, corresponding to
525 time series of the Australian Tourism Demand dataset (1998:01-2016:12).

\code{vnaggmat} is the \eqn{(221 \times 304)} aggregation matrix.
}
\source{
\url{https://robjhyndman.com/publications/mint/}
}
\usage{
# 525 time series of the Australian Tourism Demand dataset
vndata

# aggregation matrix
vnaggmat
}
\description{
The Australian Tourism Demand dataset (Wickramasuriya et al. 2019) measures the number of
nights Australians spent away from home. It includes 228 monthly observations of Visitor
Nights (VNs) from January 1998 to December 2016, and has a cross-sectional grouped
structure based on a geographic hierarchy crossed by purpose of travel. The geographic
hierarchy comprises 7 states, 27 zones, and 76 regions, for a total of 111 nested geographic
divisions. Six of these zones are each formed by a single region, resulting in 105 unique
nodes in the hierarchy. The purpose of travel comprises four categories: holiday, visiting
friends and relatives, business, and other. To avoid redundancies (Girolimetto et al. 2023),
24 nodes (6 zones are formed by a single region) are not considered, resulting in an unbalanced
hierarchy of 525 (304 bottom and 221 upper time series) unique nodes instead of the theoretical
555 with duplicated nodes.
}
\references{
Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J. (2024),
Cross-temporal probabilistic forecast reconciliation: Methodological and
practical issues. \emph{International Journal of Forecasting},  40, 3, 1134-1151.
\doi{10.1016/j.ijforecast.2023.10.003}

Wickramasuriya, S.L., Athanasopoulos, G. and Hyndman, R.J. (2019), Optimal forecast
reconciliation for hierarchical and grouped time series through trace minimization,
\emph{Journal of the American Statistical Association}, 114, 526, 804-819.
\doi{10.1080/01621459.2018.1448825}
}
\keyword{datasets}
