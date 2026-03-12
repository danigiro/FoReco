# Convert a horizon-stacked matrix to a temporal form vector

This function arranges a (\\h \times (k^\ast+m)\\) "horizon-stacked"
matrix into a (\\h(k^\ast + m) \times 1\\) temporal form vector by
reshaping horizons so that they appear consecutively for each temporal
block.

## Usage

``` r
reshape2vector(hmat, agg_order)
```

## Arguments

- hmat:

  A (\\h \times (k^\ast+m)\\) numeric matrix in horizon-stacked form;
  \\m\\ is the max. order of temporal aggregation, \\k^\ast\\ is the sum
  of (a subset of) (\\p-1\\) factors of \\m\\, excluding \\m\\, and
  \\h\\ is the forecast horizon for the lowest frequency time series.
  The row identifies the forecast horizon, and the forecasts in each row
  are ordered from the lowest frequency (most temporally aggregated) to
  the highest frequency.

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

## Value

A (\\h(k^\ast + m) \times 1\\) numeric vector in temporal form: values
are ordered from the lowest frequency to the highest frequency.

## Examples

``` r
# Example with small dimensions
h <- 2   # horizons
m <- 4   # temporal aggregation order
kt <- tetools(m)$dim["kt"]

# Build a horizon-stacked matrix: h rows, k_t columns
hmat <- matrix(seq_len(h * kt), nrow = h)

# Convert to temporal form: h * k_t vector
out <- reshape2vector(hmat, agg_order = m)
out
#> k-4 h-1 k-4 h-2 k-2 h-1 k-2 h-2 k-2 h-3 k-2 h-4 k-1 h-1 k-1 h-2 k-1 h-3 k-1 h-4 
#>       1       2       3       5       4       6       7       9      11      13 
#> k-1 h-5 k-1 h-6 k-1 h-7 k-1 h-8 
#>       8      10      12      14 
```
