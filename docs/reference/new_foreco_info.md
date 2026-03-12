# Low-level construction for reconcilied forecasts attribute foreco_info class

`new_foreca_info()` is the contructor for the `foreca_info` class, which
acompany the output from the reconciliation functions in the attribute
`FoReco`. This is exported for extension purposes and for expert use
only.

## Usage

``` r
new_foreco_info(x = list())
```

## Arguments

- x:

  A list of information related to the reconcilied forecasts.

## Examples

``` r
new_foreco_info(list(
  framework = "Cross-sectional",
  forecast_horizon = 1,
  comb = "shr",
  cs_n = 3,
  rfun = "csrec"
))
```
