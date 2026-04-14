# Optimal combination cross-sectional reconciliation

This function performs optimal (in least squares sense) combination
cross-sectional forecast reconciliation for a linearly constrained
(e.g., hierarchical/grouped) multiple time series (Wickramasuriya et
al., 2019, Panagiotelis et al., 2022, Girolimetto and Di Fonzo, 2023).
The reconciled forecasts are calculated using either a projection
approach (Byron, 1978, 1979) or the equivalent structural approach by
Hyndman et al. (2011). Non-negative (Di Fonzo and Girolimetto, 2023) and
immutable including Zhang et al., 2023) reconciled forecasts can be
considered.

## Usage

``` r
csrec(base, agg_mat, cons_mat, comb = "ols", res = NULL, approach = "proj",
      nn = NULL, settings = NULL, bounds = NULL, immutable = NULL, ...)
```

## Arguments

- base:

  A (\\h \times n\\) numeric matrix or multivariate time series (`mts`
  class) containing the base forecasts to be reconciled; \\h\\ is the
  forecast horizon, and \\n\\ is the total number of time series (\\n =
  n_a + n_b\\).

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

- comb:

  A string specifying the reconciliation method. For a complete list,
  see [cscov](https://danigiro.github.io/FoReco/reference/cscov.md).

- res:

  An (\\N \times n\\) optional numeric matrix containing the in-sample
  residuals or validation errors. This matrix is used to compute some
  covariance matrices.

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
  with 3 columns (\\i,lower,upper\\), such that

  - Column 1 represents the cross-sectional series (\\i = 1, \dots,
    n\\).

  - Columns 2 and 3 indicates the *lower* and *upper* bounds,
    respectively.

- immutable:

  A numeric vector containing the column indices of the base forecasts
  (`base` parameter) that should be fixed.

- ...:

  Arguments passed on to
  [`cscov`](https://danigiro.github.io/FoReco/reference/cscov.md)

  `mse`

  :   If `TRUE` (*default*) the errors used to compute the covariance
      matrix are not mean-corrected.

  `shrink_fun`

  :   Shrinkage function of the covariance matrix,
      [shrink_estim](https://danigiro.github.io/FoReco/reference/shrink_estim.md)
      (*default*).

## Value

A (\\h \times n\\) numeric matrix of cross-sectional reconciled
forecasts.

## References

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

Girolimetto, D. and Di Fonzo, T. (2023), Point and probabilistic
forecast reconciliation for general linearly constrained multiple time
series, *Statistical Methods & Applications*, 33, 581-607.
[doi:10.1007/s10260-023-00738-6](https://doi.org/10.1007/s10260-023-00738-6)
.

Hyndman, R.J., Ahmed, R.A., Athanasopoulos, G. and Shang, H.L. (2011),
Optimal combination forecasts for hierarchical time series,
*Computational Statistics & Data Analysis*, 55, 9, 2579-2589.
[doi:10.1016/j.csda.2011.03.006](https://doi.org/10.1016/j.csda.2011.03.006)

Kourentzes, N. and Athanasopoulos, G. (2021) Elucidate structure in
intermittent demand series. *European Journal of Operational Research*,
288, 141-152.
[doi:10.1016/j.ejor.2020.05.046](https://doi.org/10.1016/j.ejor.2020.05.046)

Panagiotelis, A., Athanasopoulos, G., Gamakumara, P. and Hyndman, R.J.
(2021), Forecast reconciliation: A geometric view with new insights on
bias correction, *International Journal of Forecasting*, 37, 1, 343–359.
[doi:10.1016/j.ijforecast.2020.06.004](https://doi.org/10.1016/j.ijforecast.2020.06.004)

Stellato, B., Banjac, G., Goulart, P., Bemporad, A. and Boyd, S. (2020),
OSQP: An Operator Splitting solver for Quadratic Programs, *Mathematical
Programming Computation*, 12, 4, 637-672.
[doi:10.1007/s12532-020-00179-2](https://doi.org/10.1007/s12532-020-00179-2)

Wickramasuriya, S.L., Athanasopoulos, G. and Hyndman, R.J. (2019),
Optimal forecast reconciliation for hierarchical and grouped time series
through trace minimization, *Journal of the American Statistical
Association*, 114, 526, 804-819.
[doi:10.1080/01621459.2018.1448825](https://doi.org/10.1080/01621459.2018.1448825)

Wickramasuriya, S. L., Turlach, B. A., and Hyndman, R. J. (2020).
Optimal non-negative forecast reconciliation. *Statistics and
Computing*, 30(5), 1167–1182.
[doi:10.1007/s11222-020-09930-0](https://doi.org/10.1007/s11222-020-09930-0)

Zhang, B., Kang, Y., Panagiotelis, A. and Li, F. (2023), Optimal
reconciliation with immutable forecasts, *European Journal of
Operational Research*, 308(2), 650–660.
[doi:10.1016/j.ejor.2022.11.035](https://doi.org/10.1016/j.ejor.2022.11.035)

## See also

Regression-based reconciliation:
[`ctrec()`](https://danigiro.github.io/FoReco/reference/ctrec.md),
[`terec()`](https://danigiro.github.io/FoReco/reference/terec.md)

Cross-sectional framework:
[`csboot()`](https://danigiro.github.io/FoReco/reference/csboot.md),
[`csbu()`](https://danigiro.github.io/FoReco/reference/csbu.md),
[`cscov()`](https://danigiro.github.io/FoReco/reference/cscov.md),
[`cslcc()`](https://danigiro.github.io/FoReco/reference/cslcc.md),
[`csmo()`](https://danigiro.github.io/FoReco/reference/csmo.md),
[`csmvn()`](https://danigiro.github.io/FoReco/reference/csmvn.md),
[`cssmp()`](https://danigiro.github.io/FoReco/reference/cssmp.md),
[`cstd()`](https://danigiro.github.io/FoReco/reference/cstd.md),
[`cstools()`](https://danigiro.github.io/FoReco/reference/cstools.md)

## Examples

``` r
set.seed(123)
# (2 x 3) base forecasts matrix (simulated), Z = X + Y
base <- matrix(rnorm(6, mean = c(20, 10, 10)), 2, byrow = TRUE)
# (10 x 3) in-sample residuals matrix (simulated)
res <- t(matrix(rnorm(n = 30), nrow = 3))

# Aggregation matrix for Z = X + Y
A <- t(c(1,1))
reco <- csrec(base = base, agg_mat = A, comb = "wls", res = res)

# Zero constraints matrix for Z - X - Y = 0
C <- t(c(1, -1, -1))
reco <- csrec(base = base, cons_mat = C, comb = "wls", res = res)

# Non negative reconciliation
# Making negative one of the base forecasts for variable Y
base[1,3] <- -base[1,3]
nnreco <- csrec(base = base, agg_mat = A, comb = "wls", res = res,
                nn = "osqp")
recoinfo(nnreco, verbose = FALSE)$info
#>     obj_val   run_time iter     prim_res status status_polish
#> 1 -352.4619 3.7874e-05   25 1.242965e-15      1             1
```
