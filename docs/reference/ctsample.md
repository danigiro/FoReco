# Cross-temporal probabilistic reconciliation (sample approach)

Cross-temporal probabilistic reconciliation (sample approach)

## Usage

``` r
ctsample(sample, agg_order, fun = ctrec, ...)
```

## Arguments

- sample:

  A (\\n \times h(k^\ast + m) \times L\\) numeric array containing the
  base forecasts samples to be reconciled; \\n\\ is the total number of
  variables, \\m\\ is the max. order of temporal aggregation, \\k^\ast\\
  is the sum of (a subset of) (\\p-1\\) factors of \\m\\, excluding
  \\m\\, \\h\\ is the forecast horizon for the lowest frequency time
  series, and \\L\\ is the sample size. The row identifies a time
  series, and the forecasts in each row are ordered from the lowest
  frequency (most temporally aggregated) to the highest frequency.

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
[`tegauss`](https://danigiro.github.io/FoReco/reference/tegauss.md)`()`,
[`tesample`](https://danigiro.github.io/FoReco/reference/tesample.md)`()`

Cross-temporal framework:
[`ctboot`](https://danigiro.github.io/FoReco/reference/ctboot.md)`()`,
[`ctbu`](https://danigiro.github.io/FoReco/reference/ctbu.md)`()`,
[`ctcov`](https://danigiro.github.io/FoReco/reference/ctcov.md)`()`,
[`ctgauss`](https://danigiro.github.io/FoReco/reference/ctgauss.md)`()`,
[`ctlcc`](https://danigiro.github.io/FoReco/reference/ctlcc.md)`()`,
[`ctmo`](https://danigiro.github.io/FoReco/reference/ctmo.md)`()`,
[`ctrec`](https://danigiro.github.io/FoReco/reference/ctrec.md)`()`,
[`cttd`](https://danigiro.github.io/FoReco/reference/cttd.md)`()`,
[`cttools`](https://danigiro.github.io/FoReco/reference/cttools.md)`()`,
[`iterec`](https://danigiro.github.io/FoReco/reference/iterec.md)`()`,
[`tcsrec`](https://danigiro.github.io/FoReco/reference/heuristic-reco.md)`()`

## Examples

``` r
set.seed(123)
A <- t(c(1,1)) # Aggregation matrix for Z = X + Y
m <- 4 # from quarterly to annual temporal aggregation

# (100 x 14) base forecasts sample matrix (simulated), m = 4, h = 2, n = 3
sample <- simplify2array(lapply(1:100, function(x){
  rbind(rnorm(14, rep(c(20, 10, 5), 2*c(1, 2, 4))),
        rnorm(14, rep(c(10, 5, 2.5), 2*c(1, 2, 4))),
        rnorm(14, rep(c(10, 5, 2.5), 2*c(1, 2, 4))))
}))
# (3 x 70) in-sample residuals matrix (simulated)
res <- rbind(rnorm(70), rnorm(70), rnorm(70))

# Optimal cross-sectional probabilistic reconciliation
reco_dist_opt <- ctsample(sample, agg_order = m, agg_mat = A, res = res, comb = "bdshr")

# Bottom-up probabilistic reconciliation
reco_dist_bu <- ctsample(sample[-1,-c(1:6), ], agg_order = m, agg_mat = A, fun = ctbu)

# Level conditional coherent probabilistic reconciliation
reco_dist_lcc <- ctsample(sample, agg_order = m, agg_mat = A, fun = ctlcc)
```
