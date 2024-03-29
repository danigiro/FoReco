% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tdrec.R
\name{tdrec}
\alias{tdrec}
\title{Top-down forecast reconciliation for genuine hierarchical/grouped time series}
\usage{
tdrec(topf, C, m, weights)
}
\arguments{
\item{topf}{(\mjseqn{h \times 1}) vector of the top-level base forecast to be
disaggregated; \mjseqn{h} is the forecast horizon (for the lowest temporal
aggregation order in temporal and cross-temporal cases).}

\item{C}{(\mjseqn{n_a \times n_b}) cross-sectional (contemporaneous) matrix
mapping the \mjseqn{n_b} bottom level series into the \mjseqn{n_a} higher level ones.}

\item{m}{Highest available sampling frequency per seasonal cycle (max. order
of temporal aggregation, \mjseqn{m}), or a subset of the \mjseqn{p} factors
of \mjseqn{m}.}

\item{weights}{vector of weights to be used to disaggregate topf:
(\mjseqn{n_b \times h}) matrix in the cross-sectional framework;
(\mjseqn{m \times h}) matrix in the temporal framework;
(\mjseqn{n_b m \times h}) matrix in the cross-temporal framework.}
}
\value{
The function returns an (\mjseqn{h \times n}) matrix of
cross-sectionally reconciled forecasts, or an (\mjseqn{h(k^\ast + m) \times 1})
vector of top-down temporally reconciled forecasts, or an
(\mjseqn{n \times h (k^\ast + m)}) matrix of top-down
cross-temporally reconciled forecasts.
}
\description{
\loadmathjax
Top-down forecast reconciliation for genuine hierarchical/grouped time series,
where the forecast of a `Total' (top-level series, expected to be positive)
is disaggregated according to a proportional scheme given by a vector
of proportions (weights).
Besides the fulfillment of any aggregation constraint,
the top-down reconciled forecasts should respect two main properties:
\itemize{
  \item the top-level value remains unchanged;
  \item all the bottom time series reconciled forecasts are non-negative.
}
The top-down procedure is extended to deal with both temporal and cross-temporal cases.
Since this is a post-forecasting function, the vector of weights must be given
in input by the user, and is not calculated automatically (see Examples).
}
\details{
Fix \mjseqn{h = 1}, then
\mjsdeqn{\widetilde{\mathbf{y}} = \mathbf{S}\mathbf{w}\widehat{a}_1}
where \mjseqn{\widetilde{\mathbf{y}}} is the vector of reconciled forecasts,
\mjseqn{\mathbf{S}} is the summing matrix (whose pattern depends on which type
of reconciliation is being performed), \mjseqn{\mathbf{w}} is the vector of weights,
and \mjseqn{\widehat{a}_1} is the top-level value to be disaggregated.
}
\examples{
data(FoReco_data)
### CROSS-SECTIONAL TOP-DOWN RECONCILIATION
# Cross sectional aggregation matrix
C <- FoReco_data$C
# monthly base forecasts
mbase <- FoReco2matrix(FoReco_data$base, m = 12)$k1
obs_1 <- FoReco_data$obs$k1
# average historical proportions
props <- colMeans(obs_1[1:168,-c(1:3)]/obs_1[1:168,1])
cs_td <- tdrec(topf = mbase[,1], C = C, weights = props)

### TEMPORAL TOP-DOWN RECONCILIATION
# top ts base forecasts ([lowest_freq' ...  highest_freq']')
top_obs12 <- FoReco_data$obs$k12[1:14,1]
bts_obs1 <- FoReco_data$obs$k1[1:168,1]
# average historical proportions
props <- colMeans(matrix(bts_obs1, ncol = 12, byrow = TRUE)/top_obs12)
topbase <- FoReco_data$base[1, 1]
t_td <- tdrec(topf = topbase, m = 12, weights = props)

### CROSS-TEMPORAL TOP-DOWN RECONCILIATION
top_obs <- FoReco_data$obs$k12[1:14,1]
bts_obs <- FoReco_data$obs$k1[1:168,-c(1:3)]
bts_obs <- lapply(1:5, function(x) matrix(bts_obs[,x], nrow=14, byrow = TRUE))
bts_obs <- do.call(cbind, bts_obs)
# average historical proportions
props <- colMeans(bts_obs/top_obs)
ct_td <- tdrec(topf = topbase, m = 12, C = C, weights = props)

}
\references{
Athanasopoulos, G., Ahmed, R.A., Hyndman, R.J. (2009), Hierarchical
forecasts for Australian domestic tourism, \emph{International Journal of
Forecasting}, 25, 1, 146–166.
}
\seealso{
Other reconciliation procedures: 
\code{\link{cstrec}()},
\code{\link{ctbu}()},
\code{\link{htsrec}()},
\code{\link{iterec}()},
\code{\link{lccrec}()},
\code{\link{octrec}()},
\code{\link{tcsrec}()},
\code{\link{thfrec}()}
}
\concept{reconciliation procedures}
\keyword{top-down}
