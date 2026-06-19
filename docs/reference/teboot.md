# Temporal Joint Block Bootstrap

Joint block bootstrap for generating probabilistic base forecasts that
take into account the correlation between different temporal aggregation
orders (Girolimetto et al. 2023).

## Usage

``` r
teboot(model_list, boot_size, agg_order, block_size = 1, seed = NULL,
       xreg = NULL, ...)
```

## Arguments

- model_list:

  A list of all the \\p\\ base forecasts models ordered from the lowest
  frequency (most temporally aggregated) to the highest frequency. A
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) function for
  each model has to be available and implemented according to the
  package [forecast](https://CRAN.R-project.org/package=forecast), with
  the following mandatory parameters: *object*, *innov*, *future*, and
  *nsim*.

- boot_size:

  The number of bootstrap replicates.

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

- block_size:

  Block size of the bootstrap, which is typically equivalent to the
  forecast horizon for the most temporally aggregated series.

- seed:

  An integer seed.

- xreg:

  A (\\\text{block\\size}(k^\ast+m) \times N\_{xreg}\\) optional numeric
  matrix containing the new values of `xreg` to be used for forecasting
  ordered from the lowest frequency to the highest frequency. It can
  contains `NA`s.

- ...:

  Additional arguments for the
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) function.

## Value

A (\\\text{boot\\size}\times (k^\ast+m)\text{block\\size}\\) matrix.

## References

Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J.
(2024), Cross-temporal probabilistic forecast reconciliation:
Methodological and practical issues. *International Journal of
Forecasting*, 40, 3, 1134-1151.
[doi:10.1016/j.ijforecast.2023.10.003](https://doi.org/10.1016/j.ijforecast.2023.10.003)

## See also

Bootstrap samples:
[`csboot()`](https://danigiro.github.io/FoReco/reference/csboot.md),
[`ctboot()`](https://danigiro.github.io/FoReco/reference/ctboot.md)

Temporal framework:
[`tebu()`](https://danigiro.github.io/FoReco/reference/tebu.md),
[`tecov()`](https://danigiro.github.io/FoReco/reference/tecov.md),
[`telcc()`](https://danigiro.github.io/FoReco/reference/telcc.md),
[`temo()`](https://danigiro.github.io/FoReco/reference/temo.md),
[`temvn()`](https://danigiro.github.io/FoReco/reference/temvn.md),
[`terec()`](https://danigiro.github.io/FoReco/reference/terec.md),
[`tesmp()`](https://danigiro.github.io/FoReco/reference/tesmp.md),
[`tetd()`](https://danigiro.github.io/FoReco/reference/tetd.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md)

## Examples

``` r
set.seed(123)

# Minimal example functions: each "model" stores Gaussian residuals;
# simulate() draws new innovations (innov=NULL) or uses the supplied
# ones (innov given).
simple_model <- function(res) {
  structure(list(residuals = res, sigma = sd(res)), class = "simple_model")
}
simulate.simple_model <- function(object, nsim = 1, innov = NULL,
                                  future = TRUE, ...) {
  if (is.null(innov)) {
    rnorm(nsim, mean = 0, sd = object$sigma)
  } else {
    as.numeric(innov)[seq_len(nsim)]
  }
}

# Temporal hierarchy: annual-quarterly, m = 4
# Aggregation orders k = 4, 2, 1 => k* + m = 1 + 2 + 4 = 7 models,
# ordered from the lowest frequency (annual) to the highest (quarterly).
m <- 4
kset <- c(4, 2, 1)              # k = 4 (annual), 2 (semi), 1 (quarterly)
n_obs_per_k <- 40 / kset        # residuals available at each frequency

model_list <- lapply(seq_along(kset), function(i) {
  simple_model(rnorm(n_obs_per_k[i]))
})

# Joint block bootstrap: 100 replicates, block_size = 1 (one annual step)
boot <- teboot(model_list = model_list, boot_size = 100, agg_order = m,
               block_size = 1, seed = 1)
#> Error in UseMethod("simulate"): no applicable method for 'simulate' applied to an object of class "simple_model"
```
