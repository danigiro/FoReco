# Temporal Gaussian probabilistic reconciliation

This function performs temporal probabilistic forecast reconciliation
assuming a multivariate normal base forecast distribution (Girolimetto
et al., 2024) for a single time series using temporal hierarchies
(Athanasopoulos et al., 2017).

## Usage

``` r
temvn(base, agg_order, tew = "sum", comb = "ols", res = NULL,
      approach = "proj", comb_base = comb, reduce_form = FALSE, ...)
```

## Arguments

- base:

  A (\\h(k^\ast + m) \times 1\\) numeric vector containing the base
  forecasts to be reconciled, ordered from lowest to highest frequency;
  \\m\\ is the maximum aggregation order, \\k^\ast\\ is the sum of a
  chosen subset of the \\p - 1\\ factors of \\m\\ (excluding \\m\\
  itself) and \\h\\ is the forecast horizon for the lowest frequency
  time series.

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

- tew:

  A string specifying the type of temporal aggregation. Options include:
  "`sum`" (simple summation, *default*), "`avg`" (average), "`first`"
  (first value of the period), and "`last`" (last value of the period).

- comb:

  A string specifying the reconciliation method. For a complete list,
  see [tecov](https://danigiro.github.io/FoReco/reference/tecov.md).

- res:

  A (\\N(k^\ast+m) \times 1\\) optional numeric vector containing the
  in-sample residuals or validation errors ordered from the lowest
  frequency to the highest frequency. This vector is used to compute
  come covariance matrices.

- approach:

  A string specifying the approach used to compute the reconciled
  forecasts. Options include:

  - "`proj`" (*default*): Projection approach according to Byron (1978,
    1979).

  - "`strc`": Structural approach as proposed by Hyndman et al. (2011).

  - "`proj_osqp`": Numerical solution using [osqp](https://osqp.org/)
    for projection approach.

  - "`strc_osqp`": Numerical solution using [osqp](https://osqp.org/)
    for structural approach.

- comb_base:

  A string specifying the base covariance matrix approach. For a
  complete list, see
  [tecov](https://danigiro.github.io/FoReco/reference/tecov.md). Default
  is the equal to `comb`.

- reduce_form:

  A logical parameter indicating whether the function should return the
  full distribution (`FALSE`, *default*) or only the distribution
  corresponding to the high-frequency time series (`TRUE`).

- ...:

  Arguments passed on to
  [`tecov`](https://danigiro.github.io/FoReco/reference/tecov.md)

  `mse`

  :   If `TRUE` (*default*) the errors used to compute the covariance
      matrix are not mean-corrected.

  `shrink_fun`

  :   Shrinkage function of the covariance matrix,
      [shrink_estim](https://danigiro.github.io/FoReco/reference/shrink_estim.md)
      (*default*)

## Value

A
[distributional::dist_multivariate_normal](https://pkg.mitchelloharawild.com/distributional/reference/dist_multivariate_normal.html)
object.

## References

Athanasopoulos, G., Hyndman, R.J., Kourentzes, N. and Petropoulos, F.
(2017), Forecasting with Temporal Hierarchies, *European Journal of
Operational Research*, 262, 1, 60-74.
[doi:10.1016/j.ejor.2017.02.046](https://doi.org/10.1016/j.ejor.2017.02.046)

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

## See also

Probabilistic reconciliation:
[`csmvn()`](https://danigiro.github.io/FoReco/reference/csmvn.md),
[`cssmp()`](https://danigiro.github.io/FoReco/reference/cssmp.md),
[`ctmvn()`](https://danigiro.github.io/FoReco/reference/ctmvn.md),
[`ctsmp()`](https://danigiro.github.io/FoReco/reference/ctsmp.md),
[`tesmp()`](https://danigiro.github.io/FoReco/reference/tesmp.md)

Temporal framework:
[`teboot()`](https://danigiro.github.io/FoReco/reference/teboot.md),
[`tebu()`](https://danigiro.github.io/FoReco/reference/tebu.md),
[`tecov()`](https://danigiro.github.io/FoReco/reference/tecov.md),
[`telcc()`](https://danigiro.github.io/FoReco/reference/telcc.md),
[`temo()`](https://danigiro.github.io/FoReco/reference/temo.md),
[`terec()`](https://danigiro.github.io/FoReco/reference/terec.md),
[`tesmp()`](https://danigiro.github.io/FoReco/reference/tesmp.md),
[`tetd()`](https://danigiro.github.io/FoReco/reference/tetd.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md)

## Examples

``` r
set.seed(123)
# (7 x 1) base forecasts vector (simulated), m = 4
base <- rnorm(7*2, rep(c(20, 10, 5), 2*c(1, 2, 4)))
# (70 x 1) in-sample residuals vector (simulated)
res <- rnorm(70)

m <- 4 # from quarterly to annual temporal aggregation
reco_dist <- terec(base = base, agg_order = m, comb = "wlsv", res = res)
```
