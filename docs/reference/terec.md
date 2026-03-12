# Optimal combination temporal reconciliation

This function performs forecast reconciliation for a single time series
using temporal hierarchies (Athanasopoulos et al., 2017, Nystrup et al.,
2020). The reconciled forecasts can be computed using either a
projection approach (Byron, 1978, 1979) or the equivalent structural
approach by Hyndman et al. (2011). Non-negative (Di Fonzo and
Girolimetto, 2023) and immutable reconciled forecasts can be considered.

## Usage

``` r
terec(base, agg_order, tew = "sum", comb = "ols", res = NULL,
      approach = "proj", nn = NULL, settings = NULL, bounds = NULL,
      immutable = NULL, ...)
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

- nn:

  A string specifying the algorithm to compute non-negative forecasts:

  - "`osqp`": quadratic programming optimization
    ([osqp](https://osqp.org/) solver, Girolimetto 2025).

  - "`bpv`": block principal pivoting algorithm (Wickramasuriya et al.,
    2020).

  - "`nfca`": negative forecasts correction algorithm (Kourentzes and
    Athanasopoulos, 2021; Girolimetto 2025).

  - "`nnic`": iterative non-negative reconciliation with immutable
    constraints (Girolimetto 2025).

  - "`sntz`": heuristic "set-negative-to-zero" (Di Fonzo and
    Girolimetto, 2023; Girolimetto 2025).

- settings:

  A list of control parameters.

  - `nn = "osqp"` An object of class `osqpSettings` specifying settings
    for the [osqp](https://osqp.org/) solver. For details, refer to the
    [osqp documentation](https://osqp.org/) (Stellato et al., 2020)

  - `nn = "bpv"`

    - `ptype = "fixed"`: permutation method: "`random`" or "`fixed`"

    - `par = 10`: the number of full exchange rules that may be
      attempted

    - `tol = sqrt(.Machine$double.eps)`: the tolerance criteria

    - `gtol = sqrt(.Machine$double.eps)`: the gradient tolerance
      criteria

    - `itmax = 100`: the maximum number of algorithm iterations

  - `nn = "nfca"` and `nn = "nnic"`

    - `tol = sqrt(.Machine$double.eps)`: the tolerance criteria

    - `itmax = 100`: the maximum number of algorithm iterations

  - `nn = "sntz"`

    - `type = "bu"`: the type of set-negative-to-zero heuristic: "`bu`"
      for bottom-up, "`tdp`" for top-down proportional, "`tdsp`" for
      top-down square proportional, "`tdvw`" for top-down variance
      weighted (the `res` param is used). See Girolimetto (2025) for
      details.

    - `tol = sqrt(.Machine$double.eps)`: the tolerance identification of
      negative values

- bounds:

  A matrix (see
  [set_bounds](https://danigiro.github.io/FoReco/reference/set_bounds.md))
  with 4 columns (\\k,j,lower,upper\\), such that

  - Column 1 represents the temporal aggregation order (\\k =
    m,\dots,1\\).

  - Column 2 represents the temporal forecast horizon (\\j =
    1,\dots,m/k\\).

  - Columns 3 and 4 indicates the *lower* and *upper* bounds,
    respectively.

- immutable:

  A matrix with 2 columns (\\k,j\\), such that

  - Column 1 represents the temporal aggregation order (\\k =
    m,\dots,1\\).

  - Column 2 represents the temporal forecast horizon (\\j =
    1,\dots,m/k\\).

  For example, when working with a quarterly time series:

  - `t(c(4, 1))` - Fix the one step ahead annual forecast.

  - `t(c(1, 2))` - Fix the two step ahead quarterly forecast.

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

A (\\h(k^\ast+m) \times 1\\) numeric vector of temporal reconciled
forecasts.

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

Di Fonzo, T. and Girolimetto, D. (2023), Spatio-temporal reconciliation
of solar forecasts, *Solar Energy*, 251, 13–29.
[doi:10.1016/j.solener.2023.01.003](https://doi.org/10.1016/j.solener.2023.01.003)

Girolimetto, D. (2025), Non-negative forecast reconciliation: Optimal
methods and operational solutions. *Forecasting*, 7(4), 64;
[doi:10.3390/forecast7040064](https://doi.org/10.3390/forecast7040064)

Hyndman, R.J., Ahmed, R.A., Athanasopoulos, G. and Shang, H.L. (2011),
Optimal combination forecasts for hierarchical time series,
*Computational Statistics & Data Analysis*, 55, 9, 2579-2589.
[doi:10.1016/j.csda.2011.03.006](https://doi.org/10.1016/j.csda.2011.03.006)

Kourentzes, N. and Athanasopoulos, G. (2021) Elucidate structure in
intermittent demand series. *European Journal of Operational Research*,
288, 141-152.
[doi:10.1016/j.ejor.2020.05.046](https://doi.org/10.1016/j.ejor.2020.05.046)

Nystrup, P., Lindström, E., Pinson, P. and Madsen, H. (2020), Temporal
hierarchies with autocorrelation for load forecasting, *European Journal
of Operational Research*, 280, 1, 876-888.
[doi:10.1016/j.ejor.2019.07.061](https://doi.org/10.1016/j.ejor.2019.07.061)

Stellato, B., Banjac, G., Goulart, P., Bemporad, A. and Boyd, S. (2020),
OSQP: An Operator Splitting solver for Quadratic Programs, *Mathematical
Programming Computation*, 12, 4, 637-672.
[doi:10.1007/s12532-020-00179-2](https://doi.org/10.1007/s12532-020-00179-2)

Wickramasuriya, S. L., Turlach, B. A., and Hyndman, R. J. (2020).
Optimal non-negative forecast reconciliation. *Statistics and
Computing*, 30(5), 1167–1182.
[doi:10.1007/s11222-020-09930-0](https://doi.org/10.1007/s11222-020-09930-0)

## See also

Regression-based reconciliation:
[`csrec()`](https://danigiro.github.io/FoReco/reference/csrec.md),
[`ctrec()`](https://danigiro.github.io/FoReco/reference/ctrec.md)

Temporal framework:
[`teboot()`](https://danigiro.github.io/FoReco/reference/teboot.md),
[`tebu()`](https://danigiro.github.io/FoReco/reference/tebu.md),
[`tecov()`](https://danigiro.github.io/FoReco/reference/tecov.md),
[`telcc()`](https://danigiro.github.io/FoReco/reference/telcc.md),
[`temo()`](https://danigiro.github.io/FoReco/reference/temo.md),
[`temvn()`](https://danigiro.github.io/FoReco/reference/temvn.md),
[`tesmp()`](https://danigiro.github.io/FoReco/reference/tesmp.md),
[`tetd()`](https://danigiro.github.io/FoReco/reference/tetd.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md)

## Examples

``` r
set.seed(123)
# (7 x 1) base forecasts vector (simulated), m = 4
base <- rnorm(7, rep(c(20, 10, 5), c(1, 2, 4)))
# (70 x 1) in-sample residuals vector (simulated)
res <- rnorm(70)

m <- 4 # from quarterly to annual temporal aggregation
reco <- terec(base = base, agg_order = m, comb = "wlsv", res = res)

# Immutable reconciled forecast
# E.g. fix all the quarterly forecasts
imm_q <- expand.grid(k = 1, j = 1:4)
immreco <- terec(base = base, agg_order = m, comb = "wlsv",
                 res = res, immutable = imm_q)

# Non negative reconciliation
base[7] <- -base[7] # Making negative one of the quarterly base forecasts
nnreco <- terec(base = base, agg_order = m, comb = "wlsv",
                res = res, nn = "osqp")
recoinfo(nnreco, verbose = FALSE)$info
#>     obj_val   run_time iter     prim_res status status_polish
#> 1 -421.8914 3.9707e-05   25 9.816134e-16      1             1
```
