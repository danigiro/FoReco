% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reco_heu.R
\name{tcsrec}
\alias{tcsrec}
\alias{cstrec}
\title{Heuristic cross-temporal reconciliation}
\usage{
# First-temporal-then-cross-sectional forecast reconciliation
tcsrec(base, cslist, telist, res = NULL, avg = "KA")

# First-cross-sectional-then-temporal forecast reconciliation
cstrec(base, cslist, telist, res = NULL)
}
\arguments{
\item{base}{A (\eqn{n \times h(k^\ast+m)}) numeric matrix containing the base forecasts to
be reconciled; \eqn{n} is the total number of variables, \eqn{m} is the max. order of temporal
aggregation, \eqn{k^\ast} is the sum of (a subset of) (\eqn{p-1}) factors of \eqn{m},
excluding \eqn{m}, and \eqn{h} is the forecast horizon for the lowest frequency time series.
The row identifies a time series, and the forecasts in each row are ordered from the
lowest frequency (most temporally aggregated) to the highest frequency.}

\item{cslist}{A list of elements for the cross-sectional reconciliation.
See \link{csrec} for a complete list (excluded \code{base} and \code{res}).}

\item{telist}{A list of elements for the temporal reconciliation.
See \link{terec} for a complete list (excluded \code{base} and \code{res}).}

\item{res}{A (\eqn{n \times N(k^\ast+m)}) optional numeric matrix containing the
in-sample residuals at all the temporal frequencies ordered from the lowest frequency
to the highest frequency (columns) for each variable (rows). This matrix is used
to compute some covariance matrices.}

\item{avg}{If \code{avg = "KA"} (\emph{default}), the final projection
matrix \eqn{\mathbf{M}} is the one proposed by Kourentzes and
Athanasopoulos (2019), otherwise it is calculated as simple average of
all the involved projection matrices at step 2 of the procedure (see
Di Fonzo and Girolimetto, 2023).}
}
\value{
A (\eqn{n \times h(k^\ast+m)}) numeric matrix of cross-temporal reconciled forecasts.
}
\description{
\link{tcsrec} replicates the procedure by Kourentzes and Athanasopoulos (2019):
(i) for each time series the forecasts at any temporal aggregation order are
reconciled using temporal hierarchies; (ii) time-by-time cross-sectional
reconciliation is performed; and (iii) the projection matrices obtained at
step (ii) are then averaged and used to cross-sectionally reconcile the
forecasts obtained at step (i). In \link{cstrec}, the order of application of the
two reconciliation steps (temporal first, then cross-sectional), is inverted
compared to \link{tcsrec} (Di Fonzo and Girolimetto, 2023).
}
\section{Warning}{

The two-step heuristic reconciliation allows considering
non negativity constraints only in the first step. This means that non-negativity
is not guaranteed in the final reconciled values.
}

\examples{
set.seed(123)
# (3 x 7) base forecasts matrix (simulated), Z = X + Y and m = 4
base <- rbind(rnorm(7, rep(c(20, 10, 5), c(1, 2, 4))),
              rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))),
              rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))))
# (3 x 70) in-sample residuals matrix (simulated)
res <- rbind(rnorm(70), rnorm(70), rnorm(70))

A <- t(c(1,1)) # Aggregation matrix for Z = X + Y
m <- 4 # from quarterly to annual temporal aggregation

rtcs <- tcsrec(base = base,
               cslist = list(agg_mat = A, comb = "shr"),
               telist = list(agg_order = m, comb = "wlsv"),
               res = res)

rcst <- tcsrec(base = base,
               cslist = list(agg_mat = A, comb = "shr"),
               telist = list(agg_order = m, comb = "wlsv"),
               res = res)

}
\references{
Di Fonzo, T. and Girolimetto, D. (2023), Cross-temporal forecast reconciliation:
Optimal combination method and heuristic alternatives, \emph{International Journal
of Forecasting}, 39, 1, 39-57. \doi{10.1016/j.ijforecast.2021.08.004}

Kourentzes, N. and Athanasopoulos, G. (2019), Cross-temporal coherent forecasts
for Australian tourism, \emph{Annals of Tourism Research}, 75, 393-409.
\doi{10.1016/j.annals.2019.02.001}
}
\seealso{
Cross-temporal framework: 
\code{\link{ctboot}()},
\code{\link{ctbu}()},
\code{\link{ctcov}()},
\code{\link{ctlcc}()},
\code{\link{ctmo}()},
\code{\link{ctrec}()},
\code{\link{cttd}()},
\code{\link{cttools}()},
\code{\link{iterec}()}
}
\concept{Framework: cross-temporal}
