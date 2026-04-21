# Cross-sectional forecast reconciliation

## Introduction

This vignette demonstrates the use of the `FoReco` package for
cross-sectional forecast reconciliation. We will work through examples
using grouped and general linearly constrained time series, showing how
to obtain base forecasts, reconcile these forecasts, and address
practical challenges such as non-negativity constraints and immutable
forecasts. We will also explore probabilistic forecast reconciliation.

### Packages

First, we load the necessary packages.

``` r

library(FoReco)   # -> To perform reconciliation
library(forecast)  # -> To obtain base forecasts
```

## `vndata`: Groupped time series

We will use the `vndata` dataset ([Wickramasuriya et al.,
2018](#ref-Wickramasuriya2019)), which contains grouped time series
data, and `vnaggmat`, which is the corresponding aggregation matrix. See
the [dataset
vignette](https://danigiro.github.io/FoReco/articles/Dataset-vndata-and-itagdp.md)
for more details.

``` r

data(vndata)      # dataset
data(vnaggmat)    # Agg mat matrix
```

### Base forecast

To obtained the base forecasts, we fit an ETS model with log
transformation to each series. We handle zeros by replacing them with
half the minimum non-zero value in the series ([Wickramasuriya et al.,
2020](#ref-Wickramasuriya2020-zk)), then fit the ETS model and generate
forecasts.

``` r

model <- setNames(vector(mode='list', length=NCOL(vndata)), colnames(vndata))
fc_obj <- setNames(vector(mode='list', length=NCOL(vndata)), colnames(vndata))

# ETS model with log transformation
ets_log <- function(x, ...){
  x[x==0] <- min(x[x!=0])/2
  ets(x, lambda = 0, ...)
}

for(i in 1:NCOL(vndata)){
  model[[i]] <- ets_log(vndata[, i])
  fc_obj[[i]] <- forecast(model[[i]], h = 12)
}
```

We extract the point forecasts and residuals from the fitted models.

``` r

# Point forecasts
base <- do.call(cbind, lapply(fc_obj, function(x) x$mean))
str(base, give.attr = FALSE)
#>  Time-Series [1:12, 1:525] from 2017 to 2018: 50651 21336 24567 29800 22846 ...
# Residuals
res <- do.call(cbind, lapply(fc_obj, residuals, type = "response"))
str(res, give.attr = FALSE)
#>  Time-Series [1:228, 1:525] from 1998 to 2017: 2143 -970 -115 133 951 ...
```

### Point forecast reconciliation

We apply various reconciliation methods to the base forecasts to ensure
they add up correctly across different levels of aggregation.

Bottom-up reconciliation ([Dunn et al., 1976](#ref-Dunn1976-kv))
aggregates forecasts from the lowest level to higher levels.

``` r

fc_bts <- base[, colnames(vnaggmat)]
rf_bu <- csbu(fc_bts, agg_mat = vnaggmat)
str(rf_bu, give.attr = FALSE)
#>  num [1:12, 1:525] 44094 18218 20585 25563 19385 ...
```

In the top-down reconciliation for genuine hierarchical/grouped time
series ([Gross & Sohl, 1990](#ref-Gross1990-uf)), the forecast of
`Total` (top-level series, expected to be positive) is disaggregated
according to a proportional scheme (weights) such that:

- the top-level value remains unchanged;

- all the bottom time series reconciled forecasts are non-negative.

``` r

bts <- vndata[, colnames(vnaggmat)]
total <- vndata[, "Total"]
fc_total <- base[, "Total"]

# Average historical proportions - Gross-Sohl method A
p_gsa <- colMeans(apply(bts, 2, function(x) x/total))
rf_td_gsa <- cstd(fc_total, agg_mat = vnaggmat, weights = p_gsa)
str(rf_td_gsa, give.attr = FALSE)
#>  num [1:12, 1:525] 50651 21336 24567 29800 22846 ...

# Proportions of the historical averages - Gross-Sohl method F
p_gsf <- colMeans(bts)/mean(total)
rf_td_gsf <- cstd(fc_total, agg_mat = vnaggmat, weights = p_gsf)
str(rf_td_gsf, give.attr = FALSE)
#>  num [1:12, 1:525] 50651 21336 24567 29800 22846 ...
```

The level conditional coherent reconciliation (LCC) is a generalization
of the the original proposal by Hollyman et al.
([2021](#ref-Hollyman2021-zq)) where the reconciled forecasts are
conditional to (i.e., constrained by) the base forecasts of a specific
upper level of the hierarchy ([Di Fonzo & Girolimetto,
2024](#ref-DiFonzo2024-ijf)).

``` r

rf_lcc <- cslcc(base = base, agg_mat = vnaggmat,
                res = res, comb = "wls")
str(rf_lcc, give.attr = FALSE)
#>  num [1:12, 1:525] 47896 20405 23276 28275 21938 ...
```

Finally we can obtained the optimal (in least squares sense) combination
cross-sectional reconciled forecast ([Girolimetto & Di Fonzo,
2023](#ref-Girolimetto2023-ft); [Panagiotelis et al.,
2021](#ref-Panagiotelis2021-rh); [Wickramasuriya et al.,
2018](#ref-Wickramasuriya2019)).

``` r

rf_opt <- csrec(base = base, agg_mat = vnaggmat, res = res, comb = "shr")
str(rf_opt, give.attr = FALSE)
#>  num [1:12, 1:525] 49160 21622 24815 29433 23260 ...
```

The following table shows some options for the optimal combination
cross-sectional reconciliation function
[`csrec()`](https://danigiro.github.io/FoReco/reference/csrec.md).

[TABLE]

### Practical challenges

#### Non negativity issues

Unfortunately, our reconciled forecasts contain negative values, even
though we used non-negative base forecasts during the reconciliation.

``` r

recoinfo(rf_opt)
#> ✔ Optimal Cross-sectional Forecast Reconciliation
#> ℹ FoReco function: `csrec`
#> ℹ Covariance approximation: `shr`
#> ℹ Non-negative forecasts: `FALSE`
```

To address this issue, we can use two approaches:

- State-of-the-art numerical optimization procedure, **osqp** ([Stellato
  et al., 2020](#ref-Stellato2020-jd)).

``` r

rf_osqp <- csrec(base = base, agg_mat = vnaggmat, res = res, comb = "shr", nn = "osqp")
tmp <- recoinfo(rf_osqp)
#> ✔ Optimal Cross-sectional Forecast Reconciliation
#> ℹ FoReco function: `csrec`
#> ℹ Covariance approximation: `shr`
#> ℹ Non-negative forecasts: `TRUE`
tmp$info # OSQP information matrix 
#>     obj_val run_time iter     prim_res status status_polish
#> 1 -3197.129 0.195952  500 3.289149e-11      1             1
```

- Simple heuristic strategy: set-negative-to-zero, **sntz** ([Di Fonzo &
  Girolimetto, 2023](#ref-Di_Fonzo2023-ae)).

``` r

rf_sntz <- csrec(base = base, agg_mat = vnaggmat, res = res, comb = "shr", 
                     nn = "sntz")
recoinfo(rf_sntz)
#> ✔ Optimal Cross-sectional Forecast Reconciliation
#> ℹ FoReco function: `csrec`
#> ℹ Covariance approximation: `shr`
#> ℹ Non-negative forecasts: `TRUE`
```

#### A priori constrained (immutable) forecasts

Sometimes we may wish to incorporate a priori knowledge during the
reconciliation process ([Zhang et al., 2023](#ref-Zhang2023-yj)) in
order to improve the accuracy of the reconciled forecasts. For example,
suppose we want to fix the forecasts of the states level series at the
base forecasts values.

``` r

rf_imm <- csrec(base = base, agg_mat = vnaggmat, res = res, comb = "shr", immutable = c(2:8))
str(rf_imm, give.attr = FALSE)
#>  num [1:12, 1:525] 50403 21300 24398 29664 22962 ...
```

``` r

round(rf_imm[, 2:8] - base[, 2:8], 6)
#>          A B C D E F G
#> Jan 2017 0 0 0 0 0 0 0
#> Feb 2017 0 0 0 0 0 0 0
#> Mar 2017 0 0 0 0 0 0 0
#> Apr 2017 0 0 0 0 0 0 0
#> May 2017 0 0 0 0 0 0 0
#> Jun 2017 0 0 0 0 0 0 0
#> Jul 2017 0 0 0 0 0 0 0
#> Aug 2017 0 0 0 0 0 0 0
#> Sep 2017 0 0 0 0 0 0 0
#> Oct 2017 0 0 0 0 0 0 0
#> Nov 2017 0 0 0 0 0 0 0
#> Dec 2017 0 0 0 0 0 0 0
```

### Probabilistic forecast reconciliation

Panagiotelis et al. ([2023](#ref-Panagiotelis2023-se)) shows that a
sample from the reconciled distribution can be obtained by reconciling a
sample from the incoherent distribution. This distinction between the
incoherent sample and the reconciliation allows us to separate the two
steps.

We can use a non-parametric method, the joint block bootstrap to
simulate \\B\\ samples and then reconciled them.

``` r

# Base forecasts' sample:
# we simulate from the base models by sampling errors 
# while keeping the cross-sectional dimension fixed.
B <- 100
base_csjb <- csboot(model, B, 12)$sample
reco_csjb <- cssmp(base_csjb, agg_mat = vnaggmat, res = res, nn = "sntz",
                   comb = "shr")
class(reco_csjb) # distribution object
#> [1] "distribution" "vctrs_vctr"   "list"

# Extracts mean:
str(mean(reco_csjb), give.attr = FALSE)
#>  num [1:12, 1:525] 49489 21762 25173 29922 23554 ...
```

A parametric method assumes a normal distribution (Gaussian), to
generate the incoherent sample set of forecasts.

``` r

# Gaussian reconciled distribution with fixed base covariance matrix for different forecast horizon
reco_csg_1 <- csmvn(base = base, agg_mat = vnaggmat, comb = "shr", res = res)

# Gaussian reconciled distribution with different base covariance matrix for each forecast horizon
## Multi-step residuals
hres <- lapply(1:12, function(h) 
  sapply(model, residuals, type='response', h = h))
## List of H=12 covariance matrix (one for each forecast horizon)
cov_shr <- lapply(hres, function(r) cscov(comb = "shr", res = r)) 
reco_csg_h <- sapply(1:NROW(base), function(h){
  csmvn(base = base[h, ], agg_mat = vnaggmat, comb = "shr", res = res, comb_base = cov_shr[[h]])
})
class(reco_csg_h) <- class(reco_csg_1)
names(reco_csg_h) <- paste0("h-", 1:length(reco_csg_h))

# Reconciled sample distribution starting from gaussian base forecasts
base_css <- lapply(1:NROW(base), function(h) MASS::mvrnorm(n = B, mu = base[h, ], Sigma = cov_shr[[h]]))
base_css <- aperm(simplify2array(base_css), c(3,2,1))
dimnames(base_css) <- dimnames(base_csjb)

reco_css <- cssmp(base_css, agg_mat = vnaggmat, res = res, nn = "sntz",
                   comb = "shr")
class(reco_css) # distribution object
#> [1] "distribution" "vctrs_vctr"   "list"
```

## `itagdp`: general linearly constrained multiple time series

In this section, we work with the `itagdp` dataset and the corresponding
zero-constrained matrix `gdpconsmat`. This dataset illustrates
reconciliation under more complex linear constraints where an unique
aggregation is not available. See the [dataset
vignette](https://danigiro.github.io/FoReco/articles/Dataset-vndata-and-itagdp.md)
for more details.

``` r

data(itagdp)      # dataset
data(gdpconsmat)  # Zero-constrained matrix
```

### Base forecasts

We fit ARIMA models to each series and generate base forecasts.

``` r

model <- setNames(vector(mode='list', length=NCOL(itagdp)), colnames(itagdp))
fc_obj <- setNames(vector(mode='list', length=NCOL(itagdp)), colnames(itagdp))
for(i in 1:NCOL(itagdp)){
  model[[i]] <- auto.arima(itagdp[, i])
  fc_obj[[i]] <- forecast(model[[i]], h = 4)
}
```

``` r

# Point forecasts
base <- do.call(cbind, lapply(fc_obj, function(x) x$mean))
str(base, give.attr = FALSE)
#>  Time-Series [1:4, 1:21] from 2020 to 2021: 435117 454372 451935 483302 169348 ...
# Residuals
res <- do.call(cbind, lapply(fc_obj, function(x) x$residuals))
str(res, give.attr = FALSE)
#>  Time-Series [1:80, 1:21] from 2000 to 2020: 167.9 89.7 51.4 -43.6 -630.9 ...
```

### Point forecast reconciliation

We apply the optimal reconciliation method to the base forecasts,
considering the linear constraints defined by `gdpconsmat`.

``` r

rf_opt <- csrec(base = base, cons_mat = gdpconsmat, res = res, comb = "wls")
str(rf_opt, give.attr = FALSE)
#>  num [1:4, 1:21] 434329 453375 451278 482269 169366 ...
```

#### Practical challenge: immutable forecast

In this case, we want to fix the forecasts of the top level series
(\\GDP\\) at the base forecasts values.

``` r

rf_imm <- csrec(base = base, cons_mat = gdpconsmat, res = res, comb = "wls", immutable = c(1))
str(rf_imm, give.attr = FALSE)
#>  num [1:4, 1:21] 435117 454372 451935 483302 169424 ...
```

``` r

rf_imm[,1]-base[,1]
#>      Qtr1 Qtr2 Qtr3 Qtr4
#> 2020    0    0    0    0
```

### Probabilistic forecast reconciliation

We can use a non-parametric method, the joint block bootstrap to
simulate \\B\\ samples and then reconciled them.

``` r

# Base forecasts' sample:
# we simulate from the base models by sampling errors 
# while keeping the cross-sectional dimension fixed.
B <- 100
base_csjb <- csboot(model, B, 4)$sample
reco_csjb <- cssmp(base_csjb, cons_mat = gdpconsmat, res = res, comb = "shr")
reco_csjb # distribution object
#> <distribution[4]>
#>         h-1         h-2         h-3         h-4 
#> sample[100] sample[100] sample[100] sample[100]

# Extracts mean:
str(mean(reco_csjb), give.attr = FALSE)
#>  num [1:4, 1:21] 433771 452039 450545 481584 168979 ...
```

Alternatively, we can use a parametric method.

``` r

# Gaussian reconciled distribution with fixed base covariance matrix for different forecast horizon
reco_csg_1 <- csmvn(base = base, cons_mat = gdpconsmat, comb = "shr", res = res)

# Gaussian reconciled distribution with different base covariance matrix for each forecast horizon
## Multi-step residuals
hres <- lapply(1:4, function(h) 
  sapply(model, residuals, type='response', h = h))
## List of H=12 covariance matrix (one for each forecast horizon)
cov_shr <- lapply(hres, function(r) cscov(comb = "shr", res = r)) 

reco_csg_h <- sapply(1:NROW(base), function(h){
  csmvn(base = base[h, ], cons_mat = gdpconsmat, comb = "shr", res = res, comb_base = cov_shr[[h]])
})
class(reco_csg_h) <- class(reco_csg_1)
names(reco_csg_h) <- paste0("h-", 1:length(reco_csg_h))

# Reconciled sample distribution starting from gaussian base forecasts
base_css <- lapply(1:NROW(base), function(h) MASS::mvrnorm(n = B, mu = base[h, ], Sigma = cov_shr[[h]]))
base_css <- aperm(simplify2array(base_css), c(3,2,1))
dimnames(base_css) <- dimnames(base_csjb)

reco_css <- cssmp(base_css, cons_mat = gdpconsmat, res = res, comb = "shr")
class(reco_css) # distribution object
#> [1] "distribution" "vctrs_vctr"   "list"
```

## References

Di Fonzo, T., & Girolimetto, D. (2023). Spatio-temporal reconciliation
of solar forecasts. *Solar Energy*, *251*, 13–29.
<https://doi.org/10.1016/j.solener.2023.01.003>

Di Fonzo, T., & Girolimetto, D. (2024). Forecast combination-based
forecast reconciliation: Insights and extensions. *International Journal
of Forecasting*, *40*(2), 490–514.
<https://doi.org/10.1016/j.ijforecast.2022.07.001>

Dunn, D. M., Williams, W. H., & Dechaine, T. L. (1976). Aggregate versus
subaggregate models in local area forecasting. *Journal of the American
Statistical Association*, *71*(353), 68–71.
<https://doi.org/10.1080/01621459.1976.10481478>

Girolimetto, D., & Di Fonzo, T. (2023). Point and probabilistic forecast
reconciliation for general linearly constrained multiple time series.
*Statistical Methods & Applications*.
<https://doi.org/10.1007/s10260-023-00738-6>

Gross, C. W., & Sohl, J. E. (1990). Disaggregation methods to expedite
product line forecasting. *Journal of Forecasting*, *9*(3), 233–254.
<https://doi.org/10.1002/for.3980090304>

Hollyman, R., Petropoulos, F., & Tipping, M. E. (2021). Understanding
forecast reconciliation. *European Journal of Operational Research*,
*294*(1), 149–160. <https://doi.org/10.1016/j.ejor.2021.01.017>

Panagiotelis, A., Athanasopoulos, G., Gamakumara, P., & Hyndman, R. J.
(2021). Forecast reconciliation: A geometric view with new insights on
bias correction. *International Journal of Forecasting*, *37*(1),
343–359. <https://doi.org/10.1016/j.ijforecast.2020.06.004>

Panagiotelis, A., Gamakumara, P., Athanasopoulos, G., & Hyndman, R. J.
(2023). Probabilistic forecast reconciliation: Properties, evaluation
and score optimisation. *European Journal of Operational Research*,
*306*(2), 693–706. <https://doi.org/10.1016/j.ejor.2022.07.040>

Stellato, B., Banjac, G., Goulart, P., Bemporad, A., & Boyd, S. (2020).
OSQP: An operator splitting solver for quadratic programs. *Mathematical
Programming Computation*, *12*(4), 637–672.
<https://doi.org/10.1007/s12532-020-00179-2>

Wickramasuriya, S. L., Athanasopoulos, G., & Hyndman, R. J. (2018).
Optimal Forecast Reconciliation for Hierarchical and Grouped Time Series
Through Trace Minimization. *Journal of the American Statistical
Association*, *114*(526), 804–819.
<https://doi.org/10.1080/01621459.2018.1448825>

Wickramasuriya, S. L., Turlach, B. A., & Hyndman, R. J. (2020). Optimal
non-negative forecast reconciliation. *Statistics and Computing*,
*30*(5), 1167–1182. <https://doi.org/10.1007/s11222-020-09930-0>

Zhang, B., Kang, Y., Panagiotelis, A., & Li, F. (2023). Optimal
reconciliation with immutable forecasts. *European Journal of
Operational Research*, *308*(2), 650–660.
<https://doi.org/10.1016/j.ejor.2022.11.035>
