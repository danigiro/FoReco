# Cross-sectional top-down reconciliation

Top-down forecast reconciliation for genuine hierarchical/grouped time
series (Gross and Sohl, 1990), where the forecast of a \`Total'
(top-level series, expected to be positive) is disaggregated according
to a proportional scheme (weights). Besides fulfilling any aggregation
constraint, the top-down reconciled forecasts should respect two main
properties:

- the top-level value remains unchanged;

- all the bottom time series reconciled forecasts are non-negative.

## Usage

``` r
cstd(base, agg_mat, weights, normalize = TRUE)
```

## Arguments

- base:

  A (\\h \times 1\\) numeric vector containing the top-level base
  forecast; \\h\\ is the forecast horizon.

- agg_mat:

  A (\\n_a \times n_b\\) numeric matrix representing the cross-sectional
  aggregation matrix. It maps the \\n_b\\ bottom-level (free) variables
  into the \\n_a\\ upper (constrained) variables.

- weights:

  A (\\h \times n_b\\) numeric matrix containing the proportions for the
  bottom time series; \\h\\ is the forecast horizon, and \\n_b\\ is the
  total number of bottom variables.

- normalize:

  If `TRUE` (*default*), the `weights` will sum to 1.

## Value

A (\\h \times n\\) numeric matrix of cross-sectional reconciled
forecasts.

## References

Gross, C.W. and Sohl, J.E. (1990), Disaggregation methods to expedite
product line forecasting. *Journal of Forecasting* 9(3), 233–254.
[doi:10.1002/for.3980090304](https://doi.org/10.1002/for.3980090304)

## See also

Top-down reconciliation:
[`cttd()`](https://danigiro.github.io/FoReco/reference/cttd.md),
[`tetd()`](https://danigiro.github.io/FoReco/reference/tetd.md)

Cross-sectional framework:
[`csboot()`](https://danigiro.github.io/FoReco/reference/csboot.md),
[`csbu()`](https://danigiro.github.io/FoReco/reference/csbu.md),
[`cscov()`](https://danigiro.github.io/FoReco/reference/cscov.md),
[`cslcc()`](https://danigiro.github.io/FoReco/reference/cslcc.md),
[`csmo()`](https://danigiro.github.io/FoReco/reference/csmo.md),
[`csmvn()`](https://danigiro.github.io/FoReco/reference/csmvn.md),
[`csrec()`](https://danigiro.github.io/FoReco/reference/csrec.md),
[`cssmp()`](https://danigiro.github.io/FoReco/reference/cssmp.md),
[`cstools()`](https://danigiro.github.io/FoReco/reference/cstools.md)

## Examples

``` r
set.seed(123)
# Aggregation matrix for Z = X + Y, X = XX + XY and Y = YX + YY
A <- matrix(c(1,1,1,1,1,1,0,0,0,0,1,1), 3, byrow = TRUE)
# (3 x 1) top base forecasts vector (simulated), forecast horizon = 3
topf <- rnorm(3, 10)
# Same weights for different forecast horizons
fix_weights <- runif(4)
reco <- cstd(base = topf, agg_mat = A, weights = fix_weights)

# Different weights for different forecast horizons
h_weights <- matrix(runif(4*3), 3, 4)
recoh <- cstd(base = topf, agg_mat = A, weights = h_weights)
```
