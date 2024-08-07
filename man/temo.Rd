% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/middleout.R
\name{temo}
\alias{temo}
\title{Temporal middle-out reconciliation}
\usage{
temo(base, agg_order, order = max(agg_order), weights, tew = "sum",
     normalize = TRUE)
}
\arguments{
\item{base}{A (\eqn{hk \times 1}) numeric vector containing the temporal
aggregated base forecasts of order \eqn{k}; \eqn{k} is an aggregation
order (a factor of \eqn{m}, and \eqn{1<k<m}), \eqn{m} is the max aggregation
order, and \eqn{h} is the forecast horizon for the lowest frequency time series.}

\item{agg_order}{Highest available sampling frequency per seasonal cycle (max. order
of temporal aggregation, \eqn{m}), or a vector representing a subset of \eqn{p} factors
of \eqn{m}.}

\item{order}{The intermediate fixed aggregation order \eqn{k}.}

\item{weights}{A (\eqn{hm \times 1}) numeric vector containing the proportions for the
high-frequency time series; \eqn{m} is the max aggregation order, and \eqn{h} is the
forecast horizon for the lowest frequency time series.}

\item{tew}{A string specifying the type of temporal aggregation. Options include:
"\code{sum}" (simple summation, \emph{default}), "\code{avg}" (average),
"\code{first}" (first value of the period), and "\code{last}"
(last value of the period).}

\item{normalize}{If \code{TRUE} (\emph{default}), the \code{weights} will sum to 1.}
}
\value{
A (\eqn{h(k^\ast+m) \times 1}) numeric vector of temporal reconciled forecasts.
}
\description{
The middle-out forecast reconciliation for temporal hierarchies
combines top-down (\link{tetd}) and bottom-up (\link{tebu}) methods. Given
the base forecasts of an intermediate temporal aggregation order \eqn{k}, it performs
\itemize{
\item a top-down approach for the aggregation orders \eqn{<k};
\item a bottom-up approach for the aggregation orders \eqn{>k}.
}
}
\examples{
set.seed(123)
# (6 x 1) base forecasts vector (simulated), forecast horizon = 3
# and intermediate aggregation order k = 2 (max agg order = 4)
basek2 <- rnorm(3*2, 5)
# Same weights for different forecast horizons
fix_weights <- runif(4)
reco <- temo(base = basek2, order = 2, agg_order = 4, weights = fix_weights)

# Different weights for different forecast horizons
h_weights <- runif(4*3)
recoh <- temo(base = basek2, order = 2, agg_order = 4, weights = h_weights)

}
\seealso{
Middle-out reconciliation: 
\code{\link{csmo}()},
\code{\link{ctmo}()}

Temporal framework: 
\code{\link{teboot}()},
\code{\link{tebu}()},
\code{\link{tecov}()},
\code{\link{telcc}()},
\code{\link{terec}()},
\code{\link{tetd}()},
\code{\link{tetools}()}
}
\concept{Framework: temporal}
\concept{Reco: middle-out}
