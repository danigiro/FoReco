% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot.R
\name{teboot}
\alias{teboot}
\title{Temporal joint block bootstrap}
\usage{
teboot(model_list, boot_size, agg_order, block_size = 1, seed = NULL)
}
\arguments{
\item{model_list}{A list of all the \eqn{(k^\ast+m)} base forecasts models ordered
from the lowest frequency (most temporally aggregated) to the highest frequency.
A \code{simulate()} function for each model has to be available and implemented
according to the package \href{https://CRAN.R-project.org/package=forecast}{\pkg{forecast}},
with the following mandatory parameters: \emph{object}, \emph{innov},
\emph{future}, and \emph{nsim}.}

\item{boot_size}{The number of bootstrap replicates.}

\item{agg_order}{Highest available sampling frequency per seasonal cycle (max. order
of temporal aggregation, \eqn{m}), or a vector representing a subset of \eqn{p} factors
of \eqn{m}.}

\item{block_size}{Block size of the bootstrap, which is typically equivalent
to the forecast horizon for the most temporally aggregated series.}

\item{seed}{An integer seed.}
}
\value{
A list with two elements: the seed used to sample the errors and
a (\eqn{\text{boot\_size}\times (k^\ast+m)\text{block\_size}}) matrix.
}
\description{
Joint block bootstrap for generating probabilistic base forecasts that take into account
the correlation between different temporal aggregation orders (Girolimetto et al. 2023).
}
\references{
Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J. (2023), Cross-temporal
probabilistic forecast reconciliation: Methodological and practical issues.
\emph{International Journal of Forecasting}, 40(3), 1134-1151. \doi{10.1016/j.ijforecast.2023.10.003}
}
\seealso{
Bootstrap samples: 
\code{\link{csboot}()},
\code{\link{ctboot}()}

Temporal framework: 
\code{\link{tebu}()},
\code{\link{tecov}()},
\code{\link{telcc}()},
\code{\link{temo}()},
\code{\link{terec}()},
\code{\link{tetd}()},
\code{\link{tetools}()}
}
\concept{Bootstrap samples}
\concept{Framework: temporal}
