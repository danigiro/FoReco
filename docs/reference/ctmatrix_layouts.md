# Convert between horizon-stacked and cross-temporal layouts

These functions convert matrix between the two canonical layouts used in
cross-temporal reconciliation. Let \\m\\ be the maximum temporal
aggregation order and \\k^\ast\\ the sum of a subset of the \\(p-1)\\
proper factors of \\m\\ (excluding \\m\\); let \\h\\ be the forecast
horizon for the lowest frequency series (e.g., most aggregated temporal
forecast horizon) and \\n\\ the number of variables:

- *Horizon-stacked layout (cross-temporal version)*: a \\h \times
  n(k^\ast + m)\\ matrix where rows are the most aggregated temporal
  forecast horizons, and the values in each row are ordered from the
  lowest frequency (most temporally aggregated) to the highest frequency
  grouped by variable.

- *Cross-temporal layout*: a \\n \times h(k^\ast + m)\\ matrix where
  rows are variables, and horizons for each temporal block appear
  consecutively. rows are variables, and the values in each row are
  ordered from the lowest frequency (most temporally aggregated) to the
  highest frequency.

Then, as_ctmatrix converts a \\(h \times n(k^\ast+m))\\ horizon-stacked
to a \\(n \times h(k^\ast+m))\\ cross-temporal matrix;
as_hstack_ctlayout performs the inverse transform.

## Usage

``` r
as_ctmatrix(hmat, agg_order, n, row_names = NULL)

as_hstack_ctlayout(ctmat, agg_order)
```

## Arguments

- hmat:

  A \\h \times n(k^\ast+m)\\ numeric matrix in *horizon-stacked* layout
  (cross-temporal version).

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

- n:

  Cross-sectional number of variables.

- row_names:

  Optional character vector of length `n` with row names for the
  *cross-temporal* output of `as_ctmatrix()`. If `NULL` (*default*) no
  custom names are assigned.

- ctmat:

  A \\n \times h(k^\ast+m)\\ numeric matrix in *cross-temporal* layout.

## Value

as_ctmatrix returns a \\n \times h(k^\ast+m)\\ numeric matrix in
*cross-temporal* layout.

as_hstack_ctlayout returns a \\h \times n(k^\ast+m)\\ numeric matrix in
*horizon-stacked* layout (cross-temporal version).

## See also

Utilities:
[`FoReco2matrix()`](https://danigiro.github.io/FoReco/reference/FoReco2matrix.md),
[`aggts()`](https://danigiro.github.io/FoReco/reference/aggts.md),
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
[`set_bounds()`](https://danigiro.github.io/FoReco/reference/set_bounds.md),
[`shrink_estim()`](https://danigiro.github.io/FoReco/reference/shrink_estim.md),
[`shrink_oasd()`](https://danigiro.github.io/FoReco/reference/shrink_oasd.md),
[`teprojmat()`](https://danigiro.github.io/FoReco/reference/teprojmat.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md),
[`unbalance_hierarchy()`](https://danigiro.github.io/FoReco/reference/unbalance_hierarchy.md)

## Examples

``` r
h <- 2   # horizons
n <- 3   # variables
m <- 4   # temporal aggregation order
kt <- tetools(m)$dim["kt"]

# Build a horizon-stacked matrix: h rows, n * k_t columns
input_ct <- matrix(seq_len(h * n * kt), nrow = n, byrow = TRUE)

hmat <- as_hstack_ctlayout(input_ct, agg_order = m)
ctmat <- as_ctmatrix(hmat, agg_order = m, n = n)
# all.equal(ctmat, input_ct, check.attributes = FALSE)
```
