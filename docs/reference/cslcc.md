# Level conditional coherent reconciliation for genuine hierarchical/grouped time series

This function implements the cross-sectional forecast reconciliation
procedure that extends the original proposal by Hollyman et al. (2021).
Level conditional coherent reconciled forecasts are conditional on
(i.e., constrained by) the base forecasts of a specific upper level in
the hierarchy (exogenous constraints). It also allows handling the
linear constraints linking the variables endogenously (Di Fonzo and
Girolimetto, 2022). The function can calculate Combined Conditional
Coherent (CCC) forecasts as simple averages of Level-Conditional
Coherent (LCC) and bottom-up reconciled forecasts, with either
endogenous or exogenous constraints.

## Usage

``` r
cslcc(base, agg_mat, comb = "ols", res = NULL, approach = "proj", nn = NULL,
      settings = NULL, CCC = TRUE, const = "exogenous", bts = NULL, ...)
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

- CCC:

  A logical value indicating whether the Combined Conditional Coherent
  reconciled forecasts reconciliation should include bottom-up forecasts
  (`TRUE`, *default*), or not.

- const:

  A string specifying the reconciliation constraints:

  - "`exogenous`" (*default*): Fixes the top level of each
    sub-hierarchy.

  - "`endogenous`": Coherently revises both the top and bottom levels.

- bts:

  A (\\h \times n_b\\) numeric matrix or multivariate time series (`mts`
  class) containing bottom base forecasts defined by the user (e.g.,
  seasonal averages, as in Hollyman et al., 2021). This parameter can be
  omitted if only base forecasts are used (see Di Fonzo and Girolimetto,
  2024).

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
[`ctlcc()`](https://danigiro.github.io/FoReco/reference/ctlcc.md),
[`telcc()`](https://danigiro.github.io/FoReco/reference/telcc.md)

Cross-sectional framework:
[`csboot()`](https://danigiro.github.io/FoReco/reference/csboot.md),
[`csbu()`](https://danigiro.github.io/FoReco/reference/csbu.md),
[`cscov()`](https://danigiro.github.io/FoReco/reference/cscov.md),
[`csmo()`](https://danigiro.github.io/FoReco/reference/csmo.md),
[`csmvn()`](https://danigiro.github.io/FoReco/reference/csmvn.md),
[`csrec()`](https://danigiro.github.io/FoReco/reference/csrec.md),
[`cssmp()`](https://danigiro.github.io/FoReco/reference/cssmp.md),
[`cstd()`](https://danigiro.github.io/FoReco/reference/cstd.md),
[`cstools()`](https://danigiro.github.io/FoReco/reference/cstools.md)

## Examples

``` r
set.seed(123)
# Aggregation matrix for Z = X + Y, X = XX + XY and Y = YX + YY
A <- matrix(c(1,1,1,1,1,1,0,0,0,0,1,1), 3, byrow = TRUE)
# (2 x 7) base forecasts matrix (simulated)
base <- matrix(rnorm(7*2, mean = c(40, 20, 20, 10, 10, 10, 10)), 2, byrow = TRUE)
# (10 x 7) in-sample residuals matrix (simulated)
res <- matrix(rnorm(n = 7*10), ncol = 7)
# (2 x 7) Naive bottom base forecasts matrix: all forecasts are set equal to 10
naive <- matrix(10, 2, 4)

## EXOGENOUS CONSTRAINTS (Hollyman et al., 2021)
# Level Conditional Coherent (LCC) reconciled forecasts
exo_LC <- cslcc(base = base, agg_mat = A, comb = "wls", bts = naive,
                res = res, CCC = FALSE)

# Combined Conditional Coherent (CCC) reconciled forecasts
exo_CCC <- cslcc(base = base, agg_mat = A, comb = "wls", bts = naive,
                 res = res, CCC = TRUE)

# Results detailed by level:
# L-1: Level 1 immutable reconciled forecasts for the whole hierarchy
# L-2: Middle-Out reconciled forecasts
# L-3: Bottom-Up reconciled forecasts
info_exo <- recoinfo(exo_CCC, verbose = FALSE)
info_exo$lcc
#> $`L-1`
#>          s-1      s-2      s-3      s-4      s-5      s-6      s-7
#> h-1 39.43952 19.76757 19.67196 9.893588 9.873979 9.744387 9.927570
#> h-2 38.73494 19.47537 19.25957 9.759815 9.715556 9.423051 9.836517
#> attr(,"FoReco")
#> 
#> $`L-2`
#>          s-1      s-2      s-3      s-4      s-5       s-6       s-7
#> h-1 41.32853 19.76982 21.55871 9.894620 9.875202 11.214555 10.344153
#> h-2 38.86749 19.31315 19.55434 9.685546 9.627601  9.652737  9.901601
#> attr(,"FoReco")
#> 
#> $`L-3`
#>          s-1     s-2      s-3      s-4      s-5      s-6      s-7
#> h-1 42.37578 20.1998 22.17598 10.07051 10.12929 11.71506 10.46092
#> h-2 42.09535 21.5839 20.51145 11.22408 10.35981 10.40077 10.11068
#> attr(,"FoReco")
#> 

## ENDOGENOUS CONSTRAINTS (Di Fonzo and Girolimetto, 2024)
# Level Conditional Coherent (LCC) reconciled forecasts
endo_LC <- cslcc(base = base, agg_mat = A, comb = "wls",
                 res = res, CCC = FALSE,
                 const = "endogenous")

# Combined Conditional Coherent (CCC) reconciled forecasts
endo_CCC <- cslcc(base = base, agg_mat = A, comb = "wls",
                  res = res, CCC = TRUE,
                  const = "endogenous")

# Results detailed by level:
# L-1: Level 1 reconciled forecasts for L1 + L3 (bottom level)
# L-2: Level 2 reconciled forecasts for L2 + L3 (bottom level)
# L-3: Bottom-Up reconciled forecasts
info_endo <- recoinfo(endo_CCC, verbose = FALSE)
info_endo$lcc
#> $`L-1`
#>          s-1      s-2      s-3       s-4      s-5       s-6       s-7
#> h-1 40.23685 19.31277 20.92408  9.664411 9.648359 10.739579 10.184505
#> h-2 39.64745 20.56873 19.07871 10.759321 9.809413  9.284371  9.794343
#> attr(,"FoReco")
#> 
#> $`L-2`
#>          s-1      s-2      s-3       s-4      s-5       s-6       s-7
#> h-1 41.70862 19.94714 21.76149  9.954836 9.992300 11.392087 10.369398
#> h-2 40.11832 20.24956 19.86876 10.613199 9.636364  9.899977  9.968779
#> attr(,"FoReco")
#> 
#> $`L-3`
#>          s-1     s-2      s-3      s-4      s-5      s-6      s-7
#> h-1 42.37578 20.1998 22.17598 10.07051 10.12929 11.71506 10.46092
#> h-2 42.09535 21.5839 20.51145 11.22408 10.35981 10.40077 10.11068
#> attr(,"FoReco")
#> 
```
