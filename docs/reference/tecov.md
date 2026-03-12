# Temporal covariance matrix approximation

This function provides an approximation of the temporal base forecasts
errors covariance matrix using different reconciliation methods (see Di
Fonzo and Girolimetto, 2023).

## Usage

``` r
tecov(comb, agg_order = NULL, tew = "sum", res = NULL, mse = TRUE,
      shrink_fun = shrink_estim, ...)
```

## Arguments

- comb:

  A string specifying the covariance approximation method.

  - For ordinary least squares reconciliation:

    - "`ols`" (*default*) - identity error covariance.

  - For weighted least squares reconciliation:

    - "`str`" - structural variances.

    - "`wlsh`" - hierarchy variances (uses `res`).

    - "`wlsv`" - series variances (uses `res`).

  - For generalized least squares (uses `res`) reconciliation:

    - "`acov`" - series auto-covariance.

    - "`strar1`" - structural Markov covariance.

    - "`sar1`" - series Markov covariance.

    - "`har1`" - hierarchy Markov covariance.

    - "`shr`"/"`sam`" - shrunk/sample covariance.

  - Others (no for reconciliation):

    - "`bu`" - bottom-up covariance.

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

- tew:

  A string specifying the type of temporal aggregation. Options include:
  "`sum`" (simple summation, *default*), "`avg`" (average), "`first`"
  (first value of the period), and "`last`" (last value of the period).

- res:

  A (\\N(k^\ast+m) \times 1\\) optional numeric vector containing the
  in-sample residuals or validation errors ordered from the lowest
  frequency to the highest frequency. This vector is used to compute
  come covariance matrices.

- mse:

  If `TRUE` (*default*) the errors used to compute the covariance matrix
  are not mean-corrected.

- shrink_fun:

  Shrinkage function of the covariance matrix,
  [shrink_estim](https://danigiro.github.io/FoReco/reference/shrink_estim.md)
  (*default*)

- ...:

  Not used.

## Value

A (\\(k^\ast+m) \times (k^\ast+m)\\) symmetric matrix.

## References

Di Fonzo, T. and Girolimetto, D. (2023a), Cross-temporal forecast
reconciliation: Optimal combination method and heuristic alternatives,
*International Journal of Forecasting*, 39, 1, 39-57.
[doi:10.1016/j.ijforecast.2021.08.004](https://doi.org/10.1016/j.ijforecast.2021.08.004)

## See also

Temporal framework:
[`teboot()`](https://danigiro.github.io/FoReco/reference/teboot.md),
[`tebu()`](https://danigiro.github.io/FoReco/reference/tebu.md),
[`telcc()`](https://danigiro.github.io/FoReco/reference/telcc.md),
[`temo()`](https://danigiro.github.io/FoReco/reference/temo.md),
[`temvn()`](https://danigiro.github.io/FoReco/reference/temvn.md),
[`terec()`](https://danigiro.github.io/FoReco/reference/terec.md),
[`tesmp()`](https://danigiro.github.io/FoReco/reference/tesmp.md),
[`tetd()`](https://danigiro.github.io/FoReco/reference/tetd.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md)

## Examples

``` r
# (7 x 70) in-sample residuals matrix (simulated), agg_order = 4
res <- rnorm(70)

cov1 <- tecov("ols", agg_order = 4)                 # OLS
cov2 <- tecov("str", agg_order = 4)                 # STRC
cov3 <- tecov("wlsv", agg_order = 4, res = res)     # WLSv
cov4 <- tecov("wlsh", agg_order = 4, res = res)     # WLSh
cov5 <- tecov("acov", agg_order = 4, res = res)     # ACOV
cov6 <- tecov("strar1", agg_order = 4, res = res)   # STRAR1
cov7 <- tecov("har1", agg_order = 4, res = res)     # HAR1
cov8 <- tecov("sar1", agg_order = 4, res = res)     # SAR1
cov9 <- tecov("shr", agg_order = 4, res = res)      # SHR
cov10 <- tecov("sam", agg_order = 4, res = res)     # SAM

# Custom covariance matrix
tecov.ols2 <- function(comb, x) diag(x)
tecov(comb = "ols2", x = 7) # == tecov("ols", agg_order = 4)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> [1,]    1    0    0    0    0    0    0
#> [2,]    0    1    0    0    0    0    0
#> [3,]    0    0    1    0    0    0    0
#> [4,]    0    0    0    1    0    0    0
#> [5,]    0    0    0    0    1    0    0
#> [6,]    0    0    0    0    0    1    0
#> [7,]    0    0    0    0    0    0    1
```
