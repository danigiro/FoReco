# Set bounds for bounded forecast reconciliation

This function defines the bounds matrix considering cross-sectional,
temporal, or cross-temporal frameworks. The output matrix can be used as
input for the `bounds` parameter in functions such as
[csrec](https://danigiro.github.io/FoReco/reference/csrec.md),
[terec](https://danigiro.github.io/FoReco/reference/terec.md), or
[ctrec](https://danigiro.github.io/FoReco/reference/ctrec.md), to
perform bounded reconciliations.

## Usage

``` r
set_bounds(n, k, h, lb = -Inf, ub = Inf, approach = "osqp", bounds = NULL)
```

## Arguments

- n:

  A (\\b \times 1\\) vector representing the \\i\\th cross-sectional
  series (\\i = 1, \dots, n\\), where \\b\\ is the number of bounds to
  be set.

- k:

  A (\\b \times 1\\) vector specifying the temporal aggregation orders
  (\\k = m, \dots, 1\\).

- h:

  A (\\b \times 1\\) vector representing the forecast horizons (\\j = 1,
  \dots, m/k\\).

- lb, ub:

  A (\\b \times 1\\) vector of lower and upper bounds.

- approach:

  A string specifying the algorithm to compute bounded reconciled
  forecasts:

  - "`osqp`": quadratic programming optimization
    ([osqp](https://osqp.org/) solver).

  - "`sftb`": heuristic "set-forecasts-to-bounds", which adjusts the
    reconciled forecasts to be within specified bounds without further
    optimization.

- bounds:

  A matrix of previous bounds to be added. If not specified, new bounds
  will be computed.

## Value

A numeric matrix representing the computed bounds, which can be:

- Cross-sectional (\\b \times 3\\) matrix for cross-sectional
  reconciliation
  ([csrec](https://danigiro.github.io/FoReco/reference/csrec.md)).

- Temporal (\\b \times 4\\) matrix for temporal reconciliation
  ([terec](https://danigiro.github.io/FoReco/reference/terec.md)).

- Cross-temporal (\\b \times 5\\) matrix for cross-temporal
  reconciliation
  ([ctrec](https://danigiro.github.io/FoReco/reference/ctrec.md)).

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
[`res2matrix()`](https://danigiro.github.io/FoReco/reference/residuals.md),
[`shrink_estim()`](https://danigiro.github.io/FoReco/reference/shrink_estim.md),
[`shrink_oasd()`](https://danigiro.github.io/FoReco/reference/shrink_oasd.md),
[`teprojmat()`](https://danigiro.github.io/FoReco/reference/teprojmat.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md),
[`unbalance_hierarchy()`](https://danigiro.github.io/FoReco/reference/unbalance_hierarchy.md)

## Examples

``` r
# Example 1
# Two cross-sectional series (i = 2,3),
# with each series required to be between 0 and 1.
n <- c(2, 3)
lb <- c(0, 0)
ub <- c(1,1)
bounds_mat <- set_bounds(n = c(2, 3),
                         lb = rep(0, 2), # or lb = 0
                         ub = rep(1, 2)) # or ub = 1

# Example 2
# All the monthly values are between 0 and 1.
bounds_mat <- set_bounds(k = rep(1, 12),  # or k = 1
                         h = 1:12,
                         lb = rep(0, 12), # or lb = 0
                         ub = rep(1, 12)) # or ub = 1

# Example 3
# For two cross-sectional series (i = 2,3),
# all the monthly values are between 0 and 1.
bounds_mat <- set_bounds(n = rep(c(2, 3), each = 12),
                         k = 1,
                         h = rep(1:12, 2),
                         lb = 0, # or lb = 0
                         ub = 1) # or ub = 1
```
