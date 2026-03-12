# Cross-sectional covariance matrix approximation

This function provides an approximation of the cross-sectional base
forecasts errors covariance matrix using different reconciliation
methods (see Wickramasuriya et al., 2019 and Di Fonzo and Girolimetto,
2023).

## Usage

``` r
cscov(comb = "ols", agg_mat = NULL, res = NULL, n = NULL, mse = TRUE,
      shrink_fun = shrink_estim, ...)
```

## Arguments

- comb:

  A string specifying the covariance approximation method.

  - For ordinary least squares reconciliation:

    - "`ols`" (*default*) - identity error covariance matrix.

  - For weighted least squares reconciliation:

    - "`str`" - structural variances.

    - "`wls`" - series variances (uses `res`).

  - For generalized least squares (uses `res`) reconciliation:

    - "`shr`" - shrunk covariance (Wickramasuriya et al., 2019).

    - "`oasd`" - oracle shrunk covariance (Ando and Xiao, 2023).

    - "`sam`" - sample covariance.

  - Others (no for reconciliation):

    - "`bu`" - bottom-up covariance.

- agg_mat:

  A (\\n_a \times n_b\\) numeric matrix representing the cross-sectional
  aggregation matrix. It maps the \\n_b\\ bottom-level (free) variables
  into the \\n_a\\ upper (constrained) variables.

- res:

  An (\\N \times n\\) optional numeric matrix containing the in-sample
  residuals or validation errors. This matrix is used to compute some
  covariance matrices.

- n:

  Number of variables (\\n = n_a + n_b\\).

- mse:

  If `TRUE` (*default*) the errors used to compute the covariance matrix
  are not mean-corrected.

- shrink_fun:

  Shrinkage function of the covariance matrix,
  [shrink_estim](https://danigiro.github.io/FoReco/reference/shrink_estim.md)
  (*default*).

- ...:

  Not used.

## Value

A (\\n \times n\\) symmetric positive (semi-)definite matrix.

## References

Ando, S., and Xiao, M. (2023), High-dimensional covariance matrix
estimation: shrinkage toward a diagonal target. *IMF Working Papers*,
2023(257), A001.

Di Fonzo, T. and Girolimetto, D. (2023a), Cross-temporal forecast
reconciliation: Optimal combination method and heuristic alternatives,
*International Journal of Forecasting*, 39, 1, 39-57.
[doi:10.1016/j.ijforecast.2021.08.004](https://doi.org/10.1016/j.ijforecast.2021.08.004)

Wickramasuriya, S.L., Athanasopoulos, G. and Hyndman, R.J. (2019),
Optimal forecast reconciliation for hierarchical and grouped time series
through trace minimization, *Journal of the American Statistical
Association*, 114, 526, 804-819.
[doi:10.1080/01621459.2018.1448825](https://doi.org/10.1080/01621459.2018.1448825)

## See also

Cross-sectional framework:
[`csboot()`](https://danigiro.github.io/FoReco/reference/csboot.md),
[`csbu()`](https://danigiro.github.io/FoReco/reference/csbu.md),
[`cslcc()`](https://danigiro.github.io/FoReco/reference/cslcc.md),
[`csmo()`](https://danigiro.github.io/FoReco/reference/csmo.md),
[`csmvn()`](https://danigiro.github.io/FoReco/reference/csmvn.md),
[`csrec()`](https://danigiro.github.io/FoReco/reference/csrec.md),
[`cssmp()`](https://danigiro.github.io/FoReco/reference/cssmp.md),
[`cstd()`](https://danigiro.github.io/FoReco/reference/cstd.md),
[`cstools()`](https://danigiro.github.io/FoReco/reference/cstools.md)

## Examples

``` r
# Aggregation matrix for Z = X + Y
A <- t(c(1,1))
# (10 x 3) in-sample residuals matrix (simulated)
res <- t(matrix(rnorm(n = 30), nrow = 3))

cov1 <- cscov("ols", n = 3)          # OLS
cov2 <- cscov("str", agg_mat = A)    # STR
cov3 <- cscov("wls", res = res)      # WLS
cov4 <- cscov("shr", res = res)      # SHR
cov5 <- cscov("sam", res = res)      # SAM

# Custom covariance matrix
cscov.ols2 <- function(comb, x) diag(x)
cscov(comb = "ols2", x = 3) # == cscov("ols", n = 3)
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
```
