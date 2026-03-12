# Iterative cross-temporal reconciliation

This function performs the iterative procedure described in Di Fonzo and
Girolimetto (2023), which produces cross-temporally reconciled forecasts
by alternating forecast reconciliation along one single dimension
(either cross-sectional or temporal) at each iteration step.

## Usage

``` r
iterec(base, cslist, telist, res = NULL, itmax = 100, tol = 1e-5,
       type = "tcs", norm = "inf", verbose = TRUE)
```

## Arguments

- base:

  A (\\n \times h(k^\ast+m)\\) numeric matrix containing the base
  forecasts to be reconciled; \\n\\ is the total number of variables,
  \\m\\ is the maximum aggregation order, and \\k^\ast\\ is the sum of a
  chosen subset of the \\p - 1\\ factors of \\m\\ (excluding \\m\\
  itself), and \\h\\ is the forecast horizon for the lowest frequency
  time series. The row identifies a time series, and the forecasts in
  each row are ordered from the lowest frequency (most temporally
  aggregated) to the highest frequency.

- cslist:

  A list of elements for the cross-sectional reconciliation. See
  [csrec](https://danigiro.github.io/FoReco/reference/csrec.md) for a
  complete list (excluded `base` and `res`).

- telist:

  A list of elements for the temporal reconciliation. See
  [terec](https://danigiro.github.io/FoReco/reference/terec.md) for a
  complete list (excluded `base` and `res`).

- res:

  A (\\n \times N(k^\ast+m)\\) optional numeric matrix containing the
  in-sample residuals or validation errors ordered from the lowest
  frequency to the highest frequency (columns) for each variable (rows).
  This matrix is used to compute some covariance matrices.

- itmax:

  Max number of iteration (`100`, *default*).

- tol:

  Convergence tolerance (`1e-5`, *default*).

- type:

  A string specifying the uni-dimensional reconciliation order: temporal
  and then cross-sectional ("`tcs`") or cross-sectional and then
  temporal ("`cst`").

- norm:

  Norm used to calculate the temporal and the cross-sectional
  incoherence: infinity norm ("`inf`", *default*), one norm ("`one`"),
  and 2-norm ("`two`").

- verbose:

  If `TRUE`, reconciliation information are printed.

## Value

A (\\n \times h(k^\ast+m)\\) numeric matrix of cross-temporal reconciled
forecasts.

## References

Di Fonzo, T. and Girolimetto, D. (2023), Cross-temporal forecast
reconciliation: Optimal combination method and heuristic alternatives,
*International Journal of Forecasting*, 39, 1, 39-57.
[doi:10.1016/j.ijforecast.2021.08.004](https://doi.org/10.1016/j.ijforecast.2021.08.004)

## See also

Cross-temporal framework:
[`ctboot()`](https://danigiro.github.io/FoReco/reference/ctboot.md),
[`ctbu()`](https://danigiro.github.io/FoReco/reference/ctbu.md),
[`ctcov()`](https://danigiro.github.io/FoReco/reference/ctcov.md),
[`ctlcc()`](https://danigiro.github.io/FoReco/reference/ctlcc.md),
[`ctmo()`](https://danigiro.github.io/FoReco/reference/ctmo.md),
[`ctmvn()`](https://danigiro.github.io/FoReco/reference/ctmvn.md),
[`ctrec()`](https://danigiro.github.io/FoReco/reference/ctrec.md),
[`ctsmp()`](https://danigiro.github.io/FoReco/reference/ctsmp.md),
[`cttd()`](https://danigiro.github.io/FoReco/reference/cttd.md),
[`cttools()`](https://danigiro.github.io/FoReco/reference/cttools.md),
[`tcsrec()`](https://danigiro.github.io/FoReco/reference/heuristic-reco.md)

## Examples

``` r
set.seed(123)
# (3 x 7) base forecasts matrix (simulated), Z = X + Y and m = 4
base <- rbind(rnorm(7, rep(c(20, 10, 5), c(1, 2, 4))),
              rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))),
              rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))))
# (3 x 70) in-sample residuals matrix (simulated)
res <- rbind(rnorm(70), rnorm(70), rnorm(70))

A <- t(c(1,1)) # Aggregation matrix for Z = X + Y
m <- 4 # from quarterly to annual temporal aggregation

rite <- iterec(base = base,
               cslist = list(agg_mat = A, comb = "shr"),
               telist = list(agg_order = m, comb = "wlsv"),
               res = res)
#> ── Iterative heuristic cross-temporal forecast reconciliation ──────────────────
#> Legend: i = iteration; s = step. Norm = "inf".
#> 
#>   i.s |        Temporal | Cross-sectional |
#>     0 |            3.36 |            1.79 |
#>   1.1 |            0.00 |            1.91 |
#>   1.2 |        2.02e-01 |            0.00 |
#>   2.1 |            0.00 |        2.42e-02 |
#>   2.2 |        2.56e-03 |            0.00 |
#>   3.1 |            0.00 |        3.18e-04 |
#>   3.2 |        3.37e-05 |            0.00 |
#>   4.1 |            0.00 |        4.50e-06 |
#>   4.2 |        4.62e-07 |            0.00 |
#> 
#> ✔ Convergence achieved at iteration 4.
#> ────────────────────────────────────────────────────────────────────────────────
```
