# Cross-temporal covariance matrix approximation

This function provides an approximation of the cross-temporal base
forecasts errors covariance matrix using different reconciliation
methods (Di Fonzo and Girolimetto, 2023, and Girolimetto et al., 2023).

## Usage

``` r
ctcov(comb = "ols", agg_mat = NULL, agg_order = NULL, tew = "sum",
      res = NULL, n = NULL, mse = TRUE, shrink_fun = shrink_estim, ...)
```

## Arguments

- comb:

  A string specifying the reconciliation method.

  - For ordinary least squares reconciliation:

    - "`ols`" (*default*) - identity error covariance.

  - For weighted least squares reconciliation:

    - "`str`" - structural variances.

    - "`csstr`" - cross-sectional structural variances.

    - "`testr`" - temporal structural variances.

    - "`wlsh`" - hierarchy variances (uses `res`).

    - "`wlsv`" - series variances (uses `res`).

  - For generalized least squares (uses `res`) reconciliation:

    - "`acov`" - series auto-covariance.

    - "`bdshr`"/"`bdsam`" - shrunk/sample block diagonal cross-sectional
      covariance.

    - "`Sshr`"/"`Ssam`" - series shrunk/sample covariance.

    - "`shr`"/"`sam`" - shrunk/sample covariance.

    - "`hbshr`"/"`hbsam`" - shrunk/sample high frequency bottom time
      series covariance.

    - "`bshr`"/"`bsam`" - shrunk/sample bottom time series covariance.

    - "`hshr`"/"`hsam`" - shrunk/sample high frequency covariance.

  - Others (no for reconciliation):

    - "`bu`" - bottom-up covariance.

- agg_mat:

  A (\\n_a \times n_b\\) numeric matrix representing the cross-sectional
  aggregation matrix. It maps the \\n_b\\ bottom-level (free) variables
  into the \\n_a\\ upper (constrained) variables.

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

- tew:

  A string specifying the type of temporal aggregation. Options include:
  "`sum`" (simple summation, *default*), "`avg`" (average), "`first`"
  (first value of the period), and "`last`" (last value of the period).

- res:

  A (\\n \times N(k^\ast+m)\\) optional numeric matrix containing the
  in-sample residuals or validation errors ordered from the lowest
  frequency to the highest frequency (columns) for each variable (rows).
  This matrix is used to compute some covariance matrices.

- n:

  Cross-sectional number of variables.

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

A (\\n(k^\ast+m) \times n(k^\ast+m)\\) symmetric matrix.

## References

Di Fonzo, T. and Girolimetto, D. (2023a), Cross-temporal forecast
reconciliation: Optimal combination method and heuristic alternatives,
*International Journal of Forecasting*, 39, 1, 39-57.
[doi:10.1016/j.ijforecast.2021.08.004](https://doi.org/10.1016/j.ijforecast.2021.08.004)

Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J.
(2024), Cross-temporal probabilistic forecast reconciliation:
Methodological and practical issues. *International Journal of
Forecasting*, 40, 3, 1134-1151.
[doi:10.1016/j.ijforecast.2023.10.003](https://doi.org/10.1016/j.ijforecast.2023.10.003)

## See also

Cross-temporal framework:
[`ctboot()`](https://danigiro.github.io/FoReco/reference/ctboot.md),
[`ctbu()`](https://danigiro.github.io/FoReco/reference/ctbu.md),
[`ctlcc()`](https://danigiro.github.io/FoReco/reference/ctlcc.md),
[`ctmo()`](https://danigiro.github.io/FoReco/reference/ctmo.md),
[`ctmvn()`](https://danigiro.github.io/FoReco/reference/ctmvn.md),
[`ctrec()`](https://danigiro.github.io/FoReco/reference/ctrec.md),
[`ctsmp()`](https://danigiro.github.io/FoReco/reference/ctsmp.md),
[`cttd()`](https://danigiro.github.io/FoReco/reference/cttd.md),
[`cttools()`](https://danigiro.github.io/FoReco/reference/cttools.md),
[`iterec()`](https://danigiro.github.io/FoReco/reference/iterec.md),
[`tcsrec()`](https://danigiro.github.io/FoReco/reference/heuristic-reco.md)

## Examples

``` r
set.seed(123)
# Aggregation matrix for Z = X + Y
A <- t(c(1,1))
# (3 x 70) in-sample residuals matrix (simulated),
# agg_order = 4 (annual-quarterly)
res <- rbind(rnorm(70), rnorm(70), rnorm(70))

cov1 <- ctcov("ols", n = 3, agg_order = 4)                     # OLS
cov2 <- ctcov("str", agg_mat = A, agg_order = 4)               # STR
cov3 <- ctcov("csstr", agg_mat = A, agg_order = 4)             # CSSTR
cov4 <- ctcov("testr", n = 3, agg_order = 4)                   # TESTR
cov5 <- ctcov("wlsv", agg_order = 4, res = res)                # WLSv
cov6 <- ctcov("wlsh", agg_order = 4, res = res)                # WLSh
cov7 <- ctcov("shr", agg_order = 4, res = res)                 # SHR
cov8 <- ctcov("sam", agg_order = 4, res = res)                 # SAM
cov9 <- ctcov("acov", agg_order = 4, res = res)                # ACOV
cov10 <- ctcov("Sshr", agg_order = 4, res = res)               # Sshr
cov11 <- ctcov("Ssam", agg_order = 4, res = res)               # Ssam
cov12 <- ctcov("hshr", agg_order = 4, res = res)               # Hshr
cov13 <- ctcov("hsam", agg_order = 4, res = res)               # Hsam
cov14 <- ctcov("hbshr", agg_mat = A, agg_order = 4, res = res) # HBshr
cov15 <- ctcov("hbsam", agg_mat = A, agg_order = 4, res = res) # HBsam
cov16 <- ctcov("bshr", agg_mat = A, agg_order = 4, res = res)  # Bshr
cov17 <- ctcov("bsam", agg_mat = A, agg_order = 4, res = res)  # Bsam
cov18 <- ctcov("bdshr", agg_order = 4, res = res)              # BDshr
cov19 <- ctcov("bdsam", agg_order = 4, res = res)              # BDsam

# Custom covariance matrix
ctcov.ols2 <- function(comb, x) diag(x)
cov20 <- ctcov(comb = "ols2", x = 21) # == ctcov("ols", n = 3, agg_order = 4)
```
