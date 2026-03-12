# Temporal probabilistic reconciliation (sample approach)

Temporal probabilistic reconciliation (sample approach)

## Usage

``` r
tesample(sample, agg_order, fun = terec, ...)
```

## Arguments

- sample:

  A (\\L \times h(k^\ast + m)\\) numeric matrix containing the base
  forecasts samples to be reconciled; \\m\\ is the max aggregation
  order, \\k^\ast\\ is the sum of (a subset of) (\\p-1\\) factors of
  \\m\\, excluding \\m\\, \\h\\ is the forecast horizon for the lowest
  frequency time series, and \\L\\ is the sample size.

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

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
[`cssample`](https://danigiro.github.io/FoReco/reference/cssample.md)`()`,
[`ctgauss`](https://danigiro.github.io/FoReco/reference/ctgauss.md)`()`,
[`ctsample`](https://danigiro.github.io/FoReco/reference/ctsample.md)`()`,
[`tegauss`](https://danigiro.github.io/FoReco/reference/tegauss.md)`()`

Temporal framework:
[`teboot`](https://danigiro.github.io/FoReco/reference/teboot.md)`()`,
[`tebu`](https://danigiro.github.io/FoReco/reference/tebu.md)`()`,
[`tecov`](https://danigiro.github.io/FoReco/reference/tecov.md)`()`,
[`tegauss`](https://danigiro.github.io/FoReco/reference/tegauss.md)`()`,
[`telcc`](https://danigiro.github.io/FoReco/reference/telcc.md)`()`,
[`temo`](https://danigiro.github.io/FoReco/reference/temo.md)`()`,
[`terec`](https://danigiro.github.io/FoReco/reference/terec.md)`()`,
[`tetd`](https://danigiro.github.io/FoReco/reference/tetd.md)`()`,
[`tetools`](https://danigiro.github.io/FoReco/reference/tetools.md)`()`

## Examples

``` r
set.seed(123)
m <- 4 # from quarterly to annual temporal aggregation

# (100 x 14) base forecasts sample matrix (simulated), m = 4, h = 2
sample <- t(sapply(1:100, function(x) rnorm(14, rep(c(20, 10, 5), 2*c(1, 2, 4)))))
# (70 x 1) in-sample residuals vector (simulated)
res <- rnorm(70)

# Optimal cross-sectional probabilistic reconciliation
reco_dist_opt <- tesample(sample, agg_order = m, res = res, comb = "shr")

# Bottom-up probabilistic reconciliation
reco_dist_bu <- tesample(sample[,-c(1:6)], agg_order = m, fun = tebu)

# Level conditional coherent probabilistic reconciliation
reco_dist_lcc <- tesample(sample, agg_order = m, fun = telcc)
```
