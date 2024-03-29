% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commat.R
\name{commat}
\alias{commat}
\title{Commutation matrix}
\usage{
commat(r, c)
}
\arguments{
\item{r}{Number of rows of \mjseqn{\mathbf{Y}}.}

\item{c}{Number of columns of \mjseqn{\mathbf{Y}}.}
}
\value{
A sparse (\mjseqn{r c \times r c}) matrix, \mjseqn{\mathbf{P}}.
}
\description{
\loadmathjax
This function returns the (\mjseqn{r c \times r c})
commutation matrix \mjseqn{\mathbf{P}} such that
\mjsdeqn{\mathbf{P} \mbox{vec}(\mathbf{Y}) = \mbox{vec}(\mathbf{Y}'),}
where \mjseqn{\mathbf{Y}} is a (\mjseqn{r \times c}) matrix.
}
\examples{
Y <- matrix(rnorm(30), 5, 6)
P <- commat(5, 6)
P \%*\% as.vector(Y) == as.vector(t(Y)) # check
}
\references{
Magnus, J.R., Neudecker, H. (2019), Matrix Differential Calculus
with Applications in Statistics and Econometrics, third edition, New York,
Wiley, pp. 54-55.
}
\seealso{
Other utilities: 
\code{\link{Cmatrix}()},
\code{\link{FoReco2ts}()},
\code{\link{agg_ts}()},
\code{\link{arrange_hres}()},
\code{\link{ctf_tools}()},
\code{\link{hts_tools}()},
\code{\link{lcmat}()},
\code{\link{oct_bounds}()},
\code{\link{residuals_matrix}()},
\code{\link{score_index}()},
\code{\link{shrink_estim}()},
\code{\link{thf_tools}()}
}
\concept{utilities}
\keyword{utilities}
