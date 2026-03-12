# Cross-sectional probabilistic reconciliation (sample approach)

Cross-sectional probabilistic reconciliation (sample approach)

## Usage

``` r
cssample(sample, fun = csrec, ...)
```

## Arguments

- sample:

  A (\\h \times n \times L\\) numeric array containing the base
  forecasts samples to be reconciled; \\h\\ is the forecast horizon,
  \\n\\ is the total number of time series (\\n = n_a + n_b\\), and
  \\L\\ is the sample size.

- fun:

  A string specifying the reconciliation function to be used, as
  implemented in FoReco.

- ...:

  Arguments passed on to `fun`

## Value

A
[distributional::dist_sample](https://pkg.mitchelloharawild.com/distributional/reference/dist_sample.html)
object.

## References

Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J.
(2024), Cross-temporal probabilistic forecast reconciliation:
Methodological and practical issues. *International Journal of
Forecasting*, 40, 3, 1134-1151.
[doi:10.1016/j.ijforecast.2023.10.003](https://doi.org/10.1016/j.ijforecast.2023.10.003)

Panagiotelis, A., Gamakumara, P., Athanasopoulos, G. and Hyndman, R.J.
(2023), Probabilistic forecast reconciliation: Properties, evaluation
and score optimisation, *European Journal of Operational Research*
306(2), 693–706.
[doi:10.1016/j.ejor.2022.07.040](https://doi.org/10.1016/j.ejor.2022.07.040)

## See also

Probabilistic reconciliation:
[`csgauss`](https://danigiro.github.io/FoReco/reference/csgauss.md)`()`,
[`ctgauss`](https://danigiro.github.io/FoReco/reference/ctgauss.md)`()`,
[`ctsample`](https://danigiro.github.io/FoReco/reference/ctsample.md)`()`,
[`tegauss`](https://danigiro.github.io/FoReco/reference/tegauss.md)`()`,
[`tesample`](https://danigiro.github.io/FoReco/reference/tesample.md)`()`

Cross-sectional framework:
[`csboot`](https://danigiro.github.io/FoReco/reference/csboot.md)`()`,
[`csbu`](https://danigiro.github.io/FoReco/reference/csbu.md)`()`,
[`cscov`](https://danigiro.github.io/FoReco/reference/cscov.md)`()`,
[`csgauss`](https://danigiro.github.io/FoReco/reference/csgauss.md)`()`,
[`cslcc`](https://danigiro.github.io/FoReco/reference/cslcc.md)`()`,
[`csmo`](https://danigiro.github.io/FoReco/reference/csmo.md)`()`,
[`csrec`](https://danigiro.github.io/FoReco/reference/csrec.md)`()`,
[`cstd`](https://danigiro.github.io/FoReco/reference/cstd.md)`()`,
[`cstools`](https://danigiro.github.io/FoReco/reference/cstools.md)`()`

## Examples

``` r
set.seed(123)
A <- t(c(1,1)) # Aggregation matrix for Z = X + Y

# (100 x 3) base forecasts sample (simulated) for h = 1
base_h1 <- matrix(rnorm(100*3, mean = c(20, 10, 10)), 100, byrow = TRUE)
# (100 x 3) base forecasts sample (simulated) for h = 2
base_h2 <- matrix(rnorm(100*3, mean = c(20, 10, 10)), 100, byrow = TRUE)
# (2 x 3 x 100) base forecasts sample array with
# 2 forecast horizons, 3 time series and 100 sample
base_sample <- aperm(simplify2array(list(base_h1, base_h2)), c(3,2,1))

# Optimal cross-sectional probabilistic reconciliation
reco_dist_opt <- cssample(base_sample, agg_mat = A)

# Bottom-up probabilistic reconciliation
reco_dist_bu <- cssample(base_sample[,-1,], agg_mat = A, fun = csbu)

# Level conditional coherent probabilistic reconciliation
reco_dist_lcc <- cssample(base_sample, agg_mat = A, fun = cslcc)
```
