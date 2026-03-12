# One-step and multi-step residuals

These functions can be used to arrange residuals to reconcile temporal
or cross-temporal forecasts.

res2matrix **\[deprecated\]** takes as input a set of temporal and
cross-temporal residuals and re-organizes them into a matrix where the
rows correspond to different forecast horizons, capturing the temporal
dimension. Meanwhile, the columns are ordered based on the specific
arrangement as described in Di Fonzo and Girolimetto (2023). Please see
[as_hstack_ctlayout](https://danigiro.github.io/FoReco/reference/ctmatrix_layouts.md)
and
[as_hstack_telayout](https://danigiro.github.io/FoReco/reference/tematrix_layouts.md).

arrange_hres takes as input a list of multi-step residuals and is
designed to organize them in accordance with their time order
(Girolimetto et al. 2023). When applied, this function ensures that the
sequence of multi-step residuals aligns with the chronological order in
which they occurred.

## Usage

``` r
res2matrix(res, agg_order)

arrange_hres(list_res)
```

## Arguments

- res:

  A (\\n \times N(k^\ast+m)\\) numeric matrix (cross-temporal framework)
  or an (\\N(k^\ast+m) \times 1\\) numeric vector (temporal framework)
  representing the in-sample residuals or validation errors at all the
  temporal frequencies ordered from the lowest frequency to the highest
  frequency (columns) for each variable (rows).

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

- list_res:

  A list of \\H\\ multi-step residuals. Each element in the list can be
  either a (\\T \times 1\\) vector (temporal framework) or a (\\T \times
  n\\) matrix (cross-temporal framework).

## Value

res2matrix returns a (\\N \times n(k^\ast + m)\\) matrix, where \\n =
1\\ for the temporal framework.

arrange_hres returns a (\\N(k^\ast+m) \times 1\\) vector (temporal
framework) or a (\\n \times N(k^\ast+m)\\) matrix (cross-temporal
framework) of multi-step residuals.

## Details

Let \\Z_t\\, \\t=1,\dots,T\\, be a univariate time series. We can define
the multi-step residuals such us \$\$\widehat{\varepsilon}\_{h,t} =
Z\_{t+h} - \widehat{Z}\_{t+h\|t} \qquad h \le t \le T-h\$\$ where
\\\widehat{Z}\_{t+h\|t}\\ is the \\h\\-step fitted value, calculated as
the \\h\\-step ahead forecast condition to the information up to time
\\t\\. Given the list of errors at different steps
\$\$\left(\[\widehat{\varepsilon}\_{1,1}, \\ \dots, \\
\widehat{\varepsilon}\_{1,T}\], \dots, \[\widehat{\varepsilon}\_{H,1},
\\ \dots, \\ \widehat{\varepsilon}\_{H,T}\]\right),\$\$ arrange_hres
returns a \\T\\-vector with the residuals, organized in the following
way: \$\$\[\varepsilon\_{1,1} \\ \varepsilon\_{2,2} \\ \dots \\
\varepsilon\_{H,H} \\ \varepsilon\_{1,H+1} \\ \dots \\
\varepsilon\_{H,T-H}\]'\$\$ A similar organisation can be apply to a
multivariate time series.

## References

Di Fonzo, T. and Girolimetto, D. (2023a), Cross-temporal forecast
reconciliation: Optimal combination method and heuristic alternatives,
*International Journal of Forecasting*, 39, 1, 39-57.
[doi:10.1016/j.ijforecast.2021.08.004](https://doi.org/10.1016/j.ijforecast.2021.08.004)

Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J.
(2024), Cross-temporal probabilistic forecast reconciliation:
Methodological and practical issues. *International Journal of
Forecasting*, 40, 3, 1134-1151.
[doi:10.1016/j.ijforecast.2023.10.003](https://doi.org/10.1016/j.ijforecast.2023.10.003)

## See also

Utilities:
[`FoReco2matrix()`](https://danigiro.github.io/FoReco/reference/FoReco2matrix.md),
[`aggts()`](https://danigiro.github.io/FoReco/reference/aggts.md),
[`as_ctmatrix()`](https://danigiro.github.io/FoReco/reference/ctmatrix_layouts.md),
[`as_tevector()`](https://danigiro.github.io/FoReco/reference/tematrix_layouts.md),
[`balance_hierarchy()`](https://danigiro.github.io/FoReco/reference/balance_hierarchy.md),
[`commat()`](https://danigiro.github.io/FoReco/reference/commat.md),
[`csprojmat()`](https://danigiro.github.io/FoReco/reference/csprojmat.md),
[`cstools()`](https://danigiro.github.io/FoReco/reference/cstools.md),
[`ctprojmat()`](https://danigiro.github.io/FoReco/reference/ctprojmat.md),
[`cttools()`](https://danigiro.github.io/FoReco/reference/cttools.md),
[`df2aggmat()`](https://danigiro.github.io/FoReco/reference/df2aggmat.md),
[`lcmat()`](https://danigiro.github.io/FoReco/reference/lcmat.md),
[`recoinfo()`](https://danigiro.github.io/FoReco/reference/recoinfo.md),
[`set_bounds()`](https://danigiro.github.io/FoReco/reference/set_bounds.md),
[`shrink_estim()`](https://danigiro.github.io/FoReco/reference/shrink_estim.md),
[`shrink_oasd()`](https://danigiro.github.io/FoReco/reference/shrink_oasd.md),
[`teprojmat()`](https://danigiro.github.io/FoReco/reference/teprojmat.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md),
[`unbalance_hierarchy()`](https://danigiro.github.io/FoReco/reference/unbalance_hierarchy.md)

## Examples

``` r
h <- 10
agg_order <- 4
tmp <- tetools(agg_order)
kt <- tmp$dim["kt"]

# Simulate vector (temporal case)
vec <- rnorm(kt*h)
out <- res2matrix(vec, agg_order) # matrix h x kt
#> Warning: `res2matrix()` was deprecated in FoReco 1.0.
#> ℹ Please use `FoReco::as_hstack_telayout()` or `FoReco::as_hstack_ctlayout()`.

# Simulate (n x kt) matrix (cross-temporal case) with n = 3
mat <- rbind(rnorm(kt*h), rnorm(kt*h), rnorm(kt*h))
out <- res2matrix(mat, agg_order) # matrix h x (3*kt)

# Input: 4 (forecast horizons) vectors with 4*10 elements
input <-  list(rnorm(4*10), rnorm(4*10), rnorm(4*10), rnorm(4*10))
# Output: 1 vector with 4*10 elements
out <- arrange_hres(input)

# Matrix version
input <-  list(matrix(rnorm(4*10*3), 4*10), matrix(rnorm(4*10*3), 4*10),
               matrix(rnorm(4*10*3), 4*10), matrix(rnorm(4*10*3), 4*10))
out <- arrange_hres(input)
```
