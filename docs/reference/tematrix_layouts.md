# Convert between horizon-stacked and temporal layouts

These functions convert matrix between the two canonical layouts used in
temporal reconciliation. Let \\m\\ be the maximum temporal aggregation
order and \\k^\ast\\ the sum of a subset of the \\(p-1)\\ proper factors
of \\m\\ (excluding \\m\\); let \\h\\ be the forecast horizon for the
lowest frequency series (e.g., most aggregated temporal forecast
horizon):

- *Horizon-stacked layout (temporal version)*: a \\h \times (k^\ast +
  m)\\ matrix where rows are the most aggregated temporal forecast
  horizons, and the values in each row are ordered from the lowest
  frequency (most temporally aggregated) to the highest frequency.

- *Temporal layout*: a (\\h(k^\ast + m) \times 1\\) numeric vector where
  values are ordered from the lowest frequency (most temporally
  aggregated) to the highest frequency.

Then, as_tevector converts a \\(h \times (k^\ast+m))\\ horizon-stacked
matrix to a (\\h(k^\ast + m) \times 1\\) temporal vector;
as_hstack_telayout performs the inverse transform.

## Usage

``` r
as_tevector(hmat, agg_order)

as_hstack_telayout(tevec, agg_order)
```

## Arguments

- hmat:

  A \\h \times (k^\ast+m)\\ numeric matrix in *horizon-stacked* layout
  (temporal version).

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

- tevec:

  A (\\h(k^\ast + m) \times 1\\) numeric vector in *temporal* layout.

## Value

as_tevector returns a (\\h(k^\ast + m) \times 1\\) numeric vector in
*temporal* layout.

as_hstack_telayout returns a \\h \times (k^\ast+m)\\ numeric matrix in
*horizon-stacked* layout (temporal version).

## See also

Utilities:
[`FoReco2matrix()`](https://danigiro.github.io/FoReco/reference/FoReco2matrix.md),
[`aggts()`](https://danigiro.github.io/FoReco/reference/aggts.md),
[`as_ctmatrix()`](https://danigiro.github.io/FoReco/reference/ctmatrix_layouts.md),
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
m <- 4   # temporal aggregation order
kt <- tetools(m)$dim["kt"]

# Build a horizon-stacked matrix: h rows, n * k_t columns
input_te <- seq_len(h * kt)

hmat <- as_hstack_telayout(input_te, agg_order = m)
tevec <- as_tevector(hmat, agg_order = m)
# all.equal(tevec, input_te, check.attributes = FALSE)
```
