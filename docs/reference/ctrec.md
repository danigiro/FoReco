# Optimal combination cross-temporal reconciliation

This function performs optimal (in least squares sense) combination
cross-temporal forecast reconciliation (Di Fonzo and Girolimetto 2023a,
Girolimetto et al. 2023). The reconciled forecasts are calculated using
either a projection approach (Byron, 1978, 1979) or the equivalent
structural approach by Hyndman et al. (2011). Non-negative (Di Fonzo and
Girolimetto, 2023) and immutable reconciled forecasts can be considered.

## Usage

``` r
ctrec(base, agg_mat, cons_mat, agg_order, tew = "sum", comb = "ols",
      res = NULL, approach = "proj", nn = NULL, settings = NULL,
      bounds = NULL, immutable = NULL, ...)
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
  with 5 columns (\\i,k,j,lower,upper\\), such that

  - Column 1 represents the cross-sectional series (\\i = 1, \dots,
    n\\).

  - Column 2 represents the temporal aggregation order (\\k =
    m,\dots,1\\).

  - Column 3 represents the temporal forecast horizon (\\j =
    1,\dots,m/k\\).

  - Columns 4 and 5 indicates the *lower* and *upper* bounds,
    respectively.

- immutable:

  A matrix with three columns (\\i,k,j\\), such that

  - Column 1 represents the cross-sectional series (\\i = 1, \dots,
    n\\).

  - Column 2 represents the temporal aggregation order (\\k =
    m,\dots,1\\).

  - Column 3 represents the temporal forecast horizon (\\j =
    1,\dots,m/k\\).

  For example, when working with a quarterly multivariate time series
  (\\n = 3\\):

  - `t(c(1, 4, 1))` - Fix the one step ahead annual forecast of the
    first time series.

  - `t(c(2, 1, 2))` - Fix the two step ahead quarterly forecast of the
    second time series.

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

A (\\n \times h(k^\ast+m)\\) numeric matrix of cross-temporal reconciled
forecasts.

## References

Byron, R.P. (1978), The estimation of large social account matrices,
*Journal of the Royal Statistical Society, Series A*, 141, 3, 359-367.
[doi:10.2307/2344807](https://doi.org/10.2307/2344807)

Byron, R.P. (1979), Corrigenda: The estimation of large social account
matrices, *Journal of the Royal Statistical Society, Series A*, 142(3),
405. [doi:10.2307/2982515](https://doi.org/10.2307/2982515)

Di Fonzo, T. and Girolimetto, D. (2023a), Cross-temporal forecast
reconciliation: Optimal combination method and heuristic alternatives,
*International Journal of Forecasting*, 39, 1, 39-57.
[doi:10.1016/j.ijforecast.2021.08.004](https://doi.org/10.1016/j.ijforecast.2021.08.004)

Di Fonzo, T. and Girolimetto, D. (2023), Spatio-temporal reconciliation
of solar forecasts, *Solar Energy*, 251, 13–29.
[doi:10.1016/j.solener.2023.01.003](https://doi.org/10.1016/j.solener.2023.01.003)

Girolimetto, D. (2025), Non-negative forecast reconciliation: Optimal
methods and operational solutions. *Forecasting*, 7(4), 64;
[doi:10.3390/forecast7040064](https://doi.org/10.3390/forecast7040064)

Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J.
(2024), Cross-temporal probabilistic forecast reconciliation:
Methodological and practical issues. *International Journal of
Forecasting*, 40, 3, 1134-1151.
[doi:10.1016/j.ijforecast.2023.10.003](https://doi.org/10.1016/j.ijforecast.2023.10.003)

Hyndman, R.J., Ahmed, R.A., Athanasopoulos, G. and Shang, H.L. (2011),
Optimal combination forecasts for hierarchical time series,
*Computational Statistics & Data Analysis*, 55, 9, 2579-2589.
[doi:10.1016/j.csda.2011.03.006](https://doi.org/10.1016/j.csda.2011.03.006)

Kourentzes, N. and Athanasopoulos, G. (2021) Elucidate structure in
intermittent demand series. *European Journal of Operational Research*,
288, 141-152.
[doi:10.1016/j.ejor.2020.05.046](https://doi.org/10.1016/j.ejor.2020.05.046)

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
[`terec()`](https://danigiro.github.io/FoReco/reference/terec.md)

Cross-temporal framework:
[`ctboot()`](https://danigiro.github.io/FoReco/reference/ctboot.md),
[`ctbu()`](https://danigiro.github.io/FoReco/reference/ctbu.md),
[`ctcov()`](https://danigiro.github.io/FoReco/reference/ctcov.md),
[`ctlcc()`](https://danigiro.github.io/FoReco/reference/ctlcc.md),
[`ctmo()`](https://danigiro.github.io/FoReco/reference/ctmo.md),
[`ctmvn()`](https://danigiro.github.io/FoReco/reference/ctmvn.md),
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

A <- t(c(1,1)) # Aggregation matrix for Z = X + Y
m <- 4 # from quarterly to annual temporal aggregation
reco <- ctrec(base = base, agg_mat = A, agg_order = m,
              comb = "wlsv", res = res)

C <- t(c(1, -1, -1)) # Zero constraints matrix for Z - X - Y = 0
reco <- ctrec(base = base, cons_mat = C, agg_order = m,
              comb = "wlsv", res = res)

# Immutable reconciled forecasts
# Fix all the quarterly forecasts of the second variable.
imm_mat <- expand.grid(i = 2, k = 1, j = 1:4)
immreco <- ctrec(base = base, cons_mat = C, agg_order = m, comb = "wlsv",
                 res = res, immutable = imm_mat)

# Non negative reconciliation
# Making negative one of the quarterly base forecasts for variable X
base[2,7] <- -2*base[2,7]
nnreco <- ctrec(base = base, cons_mat = C, agg_order = m, comb = "wlsv",
                res = res, nn = "osqp")
recoinfo(nnreco, verbose = FALSE)$info
#>     obj_val   run_time iter     prim_res status status_polish
#> 1 -635.7645 8.8041e-05   25 8.881784e-16      1             1
```
