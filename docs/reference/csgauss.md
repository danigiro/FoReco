# Cross-sectional Gaussian probabilistic reconciliation

Cross-sectional Gaussian probabilistic reconciliation

## Usage

``` r
csgauss(
  base,
  agg_mat,
  cons_mat,
  comb = "ols",
  comb_base = comb,
  res = NULL,
  approach = "proj",
  reduce_form = FALSE,
  ...
)
```

## Arguments

- base:

  A (\\h \times n\\) numeric matrix or multivariate time series (`mts`
  class) containing the base forecasts to be reconciled; \\h\\ is the
  forecast horizon, and \\n\\ is the total number of time series (\\n =
  n_a + n_b\\).

- agg_mat:

  A (\\n_a \times n_b\\) numeric matrix representing the cross-sectional
  aggregation matrix. It maps the \\n_b\\ bottom-level (free) variables
  into the \\n_a\\ upper (constrained) variables.

- cons_mat:

  A (\\n_a \times n\\) numeric matrix representing the cross-sectional
  zero constraints: each row represents a constraint equation, and each
  column represents a variable. The matrix can be of full rank, meaning
  the rows are linearly independent, but this is not a strict
  requirement, as the function allows for redundancy in the constraints.

- comb:

  A string specifying the reconciliation method. For a complete list,
  see [cscov](https://danigiro.github.io/FoReco/reference/cscov.md).

- comb_base:

  A string specifying the reconciliation method. For a complete list,
  see [cscov](https://danigiro.github.io/FoReco/reference/cscov.md).

- res:

  An (\\N \times n\\) optional numeric matrix containing the in-sample
  residuals. This matrix is used to compute some covariance matrices.

- approach:

  A string specifying the approach used to compute the reconciled mean
  and covariance matrix. Options include:

  - "`proj`" (*default*): Projection approach according to Byron (1978,
    1979).

  - "`strc`": Structural approach as proposed by Hyndman et al. (2011).

- reduce_form:

  A logical parameter indicating whether the function should return the
  full distribution (`FALSE`, *default*) or only the distribution
  corresponding to the bottom-level time series (`TRUE`).

- ...:

  Arguments passed on to
  [`cscov`](https://danigiro.github.io/FoReco/reference/cscov.md)

  `mse`

  :   If `TRUE` (*default*) the residuals used to compute the covariance
      matrix are not mean-corrected.

  `shrink_fun`

  :   Shrinkage function of the covariance matrix,
      [shrink_estim](https://danigiro.github.io/FoReco/reference/shrink_estim.md)
      (*default*).

## Value

A
[distributional::dist_multivariate_normal](https://pkg.mitchelloharawild.com/distributional/reference/dist_multivariate_normal.html)
object.

## References

Byron, R.P. (1978), The estimation of large social account matrices,
*Journal of the Royal Statistical Society, Series A*, 141, 3, 359-367.
[doi:10.2307/2344807](https://doi.org/10.2307/2344807)

Byron, R.P. (1979), Corrigenda: The estimation of large social account
matrices, *Journal of the Royal Statistical Society, Series A*, 142(3),
405. [doi:10.2307/2982515](https://doi.org/10.2307/2982515)

Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J.
(2024), Cross-temporal probabilistic forecast reconciliation:
Methodological and practical issues. *International Journal of
Forecasting*, 40, 3, 1134-1151.
[doi:10.1016/j.ijforecast.2023.10.003](https://doi.org/10.1016/j.ijforecast.2023.10.003)

Hyndman, R.J., Ahmed, R.A., Athanasopoulos, G. and Shang, H.L. (2011),
Optimal combination forecasts for hierarchical time series,
*Computational Statistics & Data Analysis*, 55, 9, 2579-2589.
[doi:10.1016/j.csda.2011.03.006](https://doi.org/10.1016/j.csda.2011.03.006)

Panagiotelis, A., Gamakumara, P., Athanasopoulos, G. and Hyndman, R.J.
(2023), Probabilistic forecast reconciliation: Properties, evaluation
and score optimisation, *European Journal of Operational Research*
306(2), 693–706.
[doi:10.1016/j.ejor.2022.07.040](https://doi.org/10.1016/j.ejor.2022.07.040)

## See also

Probabilistic reconciliation:
[`cssample`](https://danigiro.github.io/FoReco/reference/cssample.md)`()`,
[`ctgauss`](https://danigiro.github.io/FoReco/reference/ctgauss.md)`()`,
[`ctsample`](https://danigiro.github.io/FoReco/reference/ctsample.md)`()`,
[`tegauss`](https://danigiro.github.io/FoReco/reference/tegauss.md)`()`,
[`tesample`](https://danigiro.github.io/FoReco/reference/tesample.md)`()`

Cross-sectional framework:
[`csboot`](https://danigiro.github.io/FoReco/reference/csboot.md)`()`,
[`csbu`](https://danigiro.github.io/FoReco/reference/csbu.md)`()`,
[`cscov`](https://danigiro.github.io/FoReco/reference/cscov.md)`()`,
[`cslcc`](https://danigiro.github.io/FoReco/reference/cslcc.md)`()`,
[`csmo`](https://danigiro.github.io/FoReco/reference/csmo.md)`()`,
[`csrec`](https://danigiro.github.io/FoReco/reference/csrec.md)`()`,
[`cssample`](https://danigiro.github.io/FoReco/reference/cssample.md)`()`,
[`cstd`](https://danigiro.github.io/FoReco/reference/cstd.md)`()`,
[`cstools`](https://danigiro.github.io/FoReco/reference/cstools.md)`()`

## Examples

``` r
set.seed(123)
# (2 x 3) base forecasts matrix (simulated), Z = X + Y
base <- matrix(rnorm(6, mean = c(20, 10, 10)), 2, byrow = TRUE)
# (10 x 3) in-sample residuals matrix (simulated)
res <- t(matrix(rnorm(n = 30), nrow = 3))

# Aggregation matrix for Z = X + Y
A <- t(c(1,1))
reco_dist <- csgauss(base = base, agg_mat = A, comb = "shr", res = res)
```
