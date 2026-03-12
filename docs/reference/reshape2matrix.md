# Convert a horizon-stacked matrix to a cross-temporal form matrix

This function arranges a (\\h \times n(k^\ast+m)\\) "horizon-stacked"
matrix into a (\\n \times h(k^\ast+m)\\) cross-temporal form matrix by
reshaping horizons so that they appear consecutively for each temporal
block.

## Usage

``` r
reshape2matrix(hmat, agg_order, n, row_names = NULL)
```

## Arguments

- hmat:

  A (\\h \times n(k^\ast+m)\\) numeric matrix in horizon-stacked form;
  \\n\\ is the total number of variables, \\m\\ is the max. order of
  temporal aggregation, \\k^\ast\\ is the sum of (a subset of) (\\p-1\\)
  factors of \\m\\, excluding \\m\\, and \\h\\ is the forecast horizon
  for the lowest frequency time series. The row identifies the forecast
  horizon, and the forecasts in each row are ordered from the lowest
  frequency (most temporally aggregated) to the highest frequency
  grouped by variable.

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

- n:

  Cross-sectional number of variables.

- row_names:

  Optional character vector of length `n` giving names for the rows of
  the output matrix. If `NULL` (*default*), no row names are assigned.

## Value

A (\\n \times h(k^\ast+m)\\) numeric matrix in cross-temporal form: the
row identifies a variable, and the forecasts in each row are ordered
from the lowest frequency (most temporally aggregated) to the highest
frequency.

## Examples

``` r
# Example with small dimensions
h <- 2   # horizons
n <- 3   # series
m <- 4   # temporal aggregation order
kt <- tetools(m)$dim["kt"]

# Build a horizon-stacked matrix: h rows, n * k_t columns
hmat <- matrix(seq_len(h * n * kt), nrow = h)

# Convert to cross-temporal layout: n rows, h * k_t columns
out <- reshape2matrix(hmat, agg_order = m, n = n)
out
#>     k-4 h-1 k-4 h-2 k-2 h-1 k-2 h-2 k-2 h-3 k-2 h-4 k-1 h-1 k-1 h-2 k-1 h-3
#> s-1       1       2       3       5       4       6       7       9      11
#> s-2      15      16      17      19      18      20      21      23      25
#> s-3      29      30      31      33      32      34      35      37      39
#>     k-1 h-4 k-1 h-5 k-1 h-6 k-1 h-7 k-1 h-8
#> s-1      13       8      10      12      14
#> s-2      27      22      24      26      28
#> s-3      41      36      38      40      42
```
