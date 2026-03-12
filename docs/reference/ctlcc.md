# Level conditional coherent reconciliation for cross-temporal hierarchies

This function implements a forecast reconciliation procedure inspired by
the original proposal by Hollyman et al. (2021) for cross-temporal
hierarchies. Level conditional coherent reconciled forecasts are
conditional on (i.e., constrained by) the base forecasts of a specific
upper level in the hierarchy (exogenous constraints). It also allows
handling the linear constraints linking the variables endogenously (Di
Fonzo and Girolimetto, 2022). The function can calculate Combined
Conditional Coherent (CCC) forecasts as simple averages of
Level-Conditional Coherent (LCC) and bottom-up reconciled forecasts,
with either endogenous or exogenous constraints.

## Usage

``` r
ctlcc(base, agg_mat, agg_order, tew = "sum", comb = "ols", res = NULL,
      approach = "proj", nn = NULL, settings = NULL, CCC = TRUE,
      const = "exogenous", hfbts = NULL, ...)
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

- CCC:

  A logical value indicating whether the Combined Conditional Coherent
  reconciled forecasts reconciliation should include bottom-up forecasts
  (`TRUE`, *default*), or not.

- const:

  A string specifying the reconciliation constraints:

  - "`exogenous`" (*default*): Fixes the top level of each
    sub-hierarchy.

  - "`endogenous`": Coherently revises both the top and bottom levels.

- hfbts:

  A (\\n \times mh\\) numeric matrix containing high frequency bottom
  base forecasts defined by the user. This parameter can be omitted if
  only base forecasts are used (see Di Fonzo and Girolimetto, 2024).

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

Di Fonzo, T. and Girolimetto, D. (2024), Forecast combination-based
forecast reconciliation: Insights and extensions, *International Journal
of Forecasting*, 40(2), 490–514.
[doi:10.1016/j.ijforecast.2022.07.001](https://doi.org/10.1016/j.ijforecast.2022.07.001)

Di Fonzo, T. and Girolimetto, D. (2023b) Spatio-temporal reconciliation
of solar forecasts. *Solar Energy* 251, 13–29.
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

Hollyman, R., Petropoulos, F. and Tipping, M.E. (2021), Understanding
forecast reconciliation. *European Journal of Operational Research*,
294, 149–160.
[doi:10.1016/j.ejor.2021.01.017](https://doi.org/10.1016/j.ejor.2021.01.017)

Stellato, B., Banjac, G., Goulart, P., Bemporad, A. and Boyd, S. (2020),
OSQP: An Operator Splitting solver for Quadratic Programs, *Mathematical
Programming Computation*, 12, 4, 637-672.
[doi:10.1007/s12532-020-00179-2](https://doi.org/10.1007/s12532-020-00179-2)

## See also

Level conditional coherent reconciliation:
[`cslcc()`](https://danigiro.github.io/FoReco/reference/cslcc.md),
[`telcc()`](https://danigiro.github.io/FoReco/reference/telcc.md)

Cross-temporal framework:
[`ctboot()`](https://danigiro.github.io/FoReco/reference/ctboot.md),
[`ctbu()`](https://danigiro.github.io/FoReco/reference/ctbu.md),
[`ctcov()`](https://danigiro.github.io/FoReco/reference/ctcov.md),
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
# Aggregation matrix for Z = X + Y, X = XX + XY and Y = YX + YY
A <- matrix(c(1,1,1,1,1,1,0,0,0,0,1,1), 3, byrow = TRUE)
# (7 x 7) base forecasts matrix (simulated), agg_order = 4
base <- rbind(rnorm(7, rep(c(40, 20, 10), c(1, 2, 4))),
              rnorm(7, rep(c(20, 10, 5), c(1, 2, 4))),
              rnorm(7, rep(c(20, 10, 5), c(1, 2, 4))),
              rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))),
              rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))),
              rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))),
              rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))))
# (7 x 70) in-sample residuals matrix (simulated)
res <- matrix(rnorm(70*7), nrow = 7)
# (4 x 4) Naive high frequency bottom base forecasts vector:
# all forecasts are set equal to 2.5
naive <- matrix(2.5, 4, 4)

## EXOGENOUS CONSTRAINTS (Hollyman et al., 2021)
# Level Conditional Coherent (LCC) reconciled forecasts
exo_LC <- ctlcc(base = base, agg_mat = A, agg_order = 4, comb = "wlsh", nn = "osqp",
                hfbts = naive, res = res, CCC = FALSE)

# Combined Conditional Coherent (CCC) reconciled forecasts
exo_CCC <- ctlcc(base = base, agg_mat = A, agg_order = 4, comb = "wlsh",
                hfbts = naive, res = res, CCC = TRUE)

# Results detailed by level:
info_exo <- recoinfo(exo_CCC, verbose = FALSE)
# info_exo$lcc

## ENDOGENOUS CONSTRAINTS (Di Fonzo and Girolimetto, 2024)
# Level Conditional Coherent (LCC) reconciled forecasts
endo_LC <- ctlcc(base = base, agg_mat = A, agg_order = 4, comb = "wlsh",
                 res = res, CCC = FALSE, const = "endogenous")

# Combined Conditional Coherent (CCC) reconciled forecasts
endo_CCC <- ctlcc(base = base, agg_mat = A, agg_order = 4, comb = "wlsh",
                  res = res, CCC = TRUE, const = "endogenous")

# Results detailed by level:
info_endo <- recoinfo(endo_CCC, verbose = FALSE)
# info_endo$lcc
```
