# Cross-temporal Gaussian probabilistic reconciliation

This function performs cross-temporal probabilistic forecast
reconciliation assuming a multivariate normal base forecast distribution
(Girolimetto et al., 2024) for linearly constrained multiple time series
observed across both cross-sectional and temporal dimensions (Di Fonzo
and Girolimetto, 2023).

## Usage

``` r
ctmvn(base, agg_mat, cons_mat, agg_order, tew = "sum", comb = "ols",
      res = NULL, approach = "proj", comb_base = comb,
      reduce_form = FALSE, ...)
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
  see [ctcov](https://danigiro.github.io/FoReco/reference/ctcov.md).

- res:

  A (\\n \times N(k^\ast+m)\\) optional numeric matrix containing the
  in-sample residuals or validation errors ordered from the lowest
  frequency to the highest frequency (columns) for each variable (rows).
  This matrix is used to compute some covariance matrices.

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
  [ctcov](https://danigiro.github.io/FoReco/reference/ctcov.md). Default
  is the equal to `comb`.

- reduce_form:

  A logical parameter indicating whether the function should return the
  full distribution (`FALSE`, *default*) or only the distribution
  corresponding to the high-frequency bottom time series (`TRUE`).

- ...:

  Arguments passed on to
  [`ctcov`](https://danigiro.github.io/FoReco/reference/ctcov.md)

  `mse`

  :   If `TRUE` (*default*) the errors used to compute the covariance
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
[`csmvn()`](https://danigiro.github.io/FoReco/reference/csmvn.md),
[`cssmp()`](https://danigiro.github.io/FoReco/reference/cssmp.md),
[`ctsmp()`](https://danigiro.github.io/FoReco/reference/ctsmp.md),
[`temvn()`](https://danigiro.github.io/FoReco/reference/temvn.md),
[`tesmp()`](https://danigiro.github.io/FoReco/reference/tesmp.md)

Cross-temporal framework:
[`ctboot()`](https://danigiro.github.io/FoReco/reference/ctboot.md),
[`ctbu()`](https://danigiro.github.io/FoReco/reference/ctbu.md),
[`ctcov()`](https://danigiro.github.io/FoReco/reference/ctcov.md),
[`ctlcc()`](https://danigiro.github.io/FoReco/reference/ctlcc.md),
[`ctmo()`](https://danigiro.github.io/FoReco/reference/ctmo.md),
[`ctrec()`](https://danigiro.github.io/FoReco/reference/ctrec.md),
[`ctsmp()`](https://danigiro.github.io/FoReco/reference/ctsmp.md),
[`cttd()`](https://danigiro.github.io/FoReco/reference/cttd.md),
[`cttools()`](https://danigiro.github.io/FoReco/reference/cttools.md),
[`iterec()`](https://danigiro.github.io/FoReco/reference/iterec.md),
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
A <- t(c(1,1))
reco_dist <- ctmvn(base = base, res = res, agg_mat = A, agg_order = 4)
#> Warning: Argument `res` is ignored.
#> ℹ When `res` is provided, `comb = 'bdshr'` is suggested.
```
