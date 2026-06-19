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
summary(rf_bu)
#> ✔ Cross-sectional point forecast reconciliation
#> 
#> ── Method
#> • Function used: `csbu`
#> • Output: (12 x 525) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Forecast horizons (h): 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>        Total         A         B        C        D        E
#> h-1 44094.08 15090.110 10589.205 8964.637 2855.082 4821.167
#> h-2 18218.45  5881.055  3793.292 3994.225 1106.625 2516.225
#> h-3 20585.12  6592.540  4482.910 4183.151 1301.891 2849.777
#> h-4 25563.48  8609.232  5174.437 5631.049 1729.393 3352.184
#> ... (8 more rows, 519 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
#> Class structure: foreco > matrix
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
summary(rf_td_gsa)
#> ✔ Cross-sectional point forecast reconciliation
#> 
#> ── Method
#> • Function used: `cstd`
#> • Output: (12 x 525) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Forecast horizons (h): 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>        Total         A        B         C        D        E
#> h-1 50651.48 16251.582 9659.413 13361.324 3417.451 5239.247
#> h-2 21335.64  6845.562 4068.780  5628.115 1439.514 2206.899
#> h-3 24566.98  7882.342 4685.009  6480.509 1657.532 2541.140
#> h-4 29799.97  9561.352 5682.958  7860.916 2010.601 3082.425
#> ... (8 more rows, 519 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
#> Class structure: foreco > foreco

# Proportions of the historical averages - Gross-Sohl method F
p_gsf <- colMeans(bts)/mean(total)
rf_td_gsf <- cstd(fc_total, agg_mat = vnaggmat, weights = p_gsf)
summary(rf_td_gsf)
#> ✔ Cross-sectional point forecast reconciliation
#> 
#> ── Method
#> • Function used: `cstd`
#> • Output: (12 x 525) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Forecast horizons (h): 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>        Total         A        B         C        D        E
#> h-1 50651.48 16356.898 9782.335 13203.994 3430.735 5201.641
#> h-2 21335.64  6889.924 4120.558  5561.844 1445.109 2191.058
#> h-3 24566.98  7933.423 4744.628  6404.201 1663.975 2522.900
#> h-4 29799.97  9623.313 5755.277  7768.354 2018.417 3060.301
#> ... (8 more rows, 519 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
#> Class structure: foreco > foreco
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
summary(rf_lcc)
#> ✔ Cross-sectional point forecast reconciliation
#> 
#> ── Method
#> • Function used: `cslcc`
#> • Covariance approximation approach: `wls`
#> • Output: (12 x 525) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 555
#> • Forecast horizons (h): 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>        Total         A         B        C        D        E
#> h-1 47896.34 16171.607 11489.265 9986.096 3107.838 5106.384
#> h-2 20404.61  6464.813  4265.185 4642.466 1270.962 2668.350
#> h-3 23275.65  7339.439  5072.770 4895.165 1509.840 3108.643
#> h-4 28275.20  9332.229  5714.468 6388.531 1894.417 3671.172
#> ... (8 more rows, 519 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
#> Class structure: foreco > matrix
```

Finally we can obtained the optimal (in least squares sense) combination
cross-sectional reconciled forecast ([Girolimetto & Di Fonzo,
2023](#ref-Girolimetto2023-ft); [Panagiotelis et al.,
2021](#ref-Panagiotelis2021-rh); [Wickramasuriya et al.,
2018](#ref-Wickramasuriya2019)).

``` r

rf_opt <- csrec(base = base, agg_mat = vnaggmat, res = res, comb = "shr")
summary(rf_opt)
#> ✔ Cross-sectional point forecast reconciliation
#> 
#> ── Method
#> • Function used: `csrec`
#> • Covariance approximation approach: `shr`
#> • Output: (12 x 525) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Forecast horizons (h): 12
#> • Non-negative forecasts (check): `FALSE`
#> 
#> ── Reconciled forecasts
#>        Total         A         B         C        D        E
#> h-1 49160.00 16147.271 11780.774 10171.782 3331.693 5378.219
#> h-2 21621.71  6665.539  4548.175  4876.768 1415.294 2809.422
#> h-3 24815.28  7597.319  5422.476  5169.599 1691.818 3334.412
#> h-4 29432.98  9453.527  5950.275  6592.491 2040.414 3892.674
#> ... (8 more rows, 519 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
#> Class structure: foreco > matrix
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

summary(rf_opt)
#> ✔ Cross-sectional point forecast reconciliation
#> 
#> ── Method
#> • Function used: `csrec`
#> • Covariance approximation approach: `shr`
#> • Output: (12 x 525) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Forecast horizons (h): 12
#> • Non-negative forecasts (check): `FALSE`
#> 
#> ── Reconciled forecasts
#>        Total         A         B         C        D        E
#> h-1 49160.00 16147.271 11780.774 10171.782 3331.693 5378.219
#> h-2 21621.71  6665.539  4548.175  4876.768 1415.294 2809.422
#> h-3 24815.28  7597.319  5422.476  5169.599 1691.818 3334.412
#> h-4 29432.98  9453.527  5950.275  6592.491 2040.414 3892.674
#> ... (8 more rows, 519 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
#> Class structure: foreco > matrix
```

To address this issue, we can use two approaches:

- State-of-the-art numerical optimization procedure, **osqp** ([Stellato
  et al., 2020](#ref-Stellato2020-jd)).

``` r

rf_osqp <- csrec(base = base, agg_mat = vnaggmat, res = res, comb = "shr", nn = "osqp")
summary(rf_osqp) # OSQP information matrix 
#> ✔ Cross-sectional point forecast reconciliation
#> 
#> ── Method
#> • Function used: `csrec`
#> • Covariance approximation approach: `shr`
#> • Output: (12 x 525) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Forecast horizons (h): 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Non-negative reconciliation diagnostics
#>     obj_val  run_time iter     prim_res status status_polish
#> 1 -3197.129 0.1774314  500 2.073845e-11      1             1
#> 
#> ── Reconciled forecasts
#>        Total         A         B         C        D        E
#> h-1 49164.22 16147.024 11781.704 10175.049 3331.750 5378.725
#> h-2 21621.71  6665.539  4548.175  4876.768 1415.294 2809.422
#> h-3 24815.28  7597.319  5422.476  5169.599 1691.818 3334.412
#> h-4 29432.98  9453.527  5950.275  6592.491 2040.414 3892.674
#> ... (8 more rows, 519 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
#> Class structure: foreco > matrix
```

- Simple heuristic strategy: set-negative-to-zero, **sntz** ([Di Fonzo &
  Girolimetto, 2023](#ref-Di_Fonzo2023-ae)).

``` r

rf_sntz <- csrec(base = base, agg_mat = vnaggmat, res = res, comb = "shr", 
                 nn = "sntz")
summary(rf_sntz)
#> ✔ Cross-sectional point forecast reconciliation
#> 
#> ── Method
#> • Function used: `csrec`
#> • Covariance approximation approach: `shr`
#> • Output: (12 x 525) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Forecast horizons (h): 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Non-negative reconciliation diagnostics
#>       run_time          tol iter
#> 1 8.821487e-06 1.490116e-08    0
#> 
#> ── Reconciled forecasts
#>        Total         A         B         C        D        E
#> h-1 49163.32 16147.271 11780.774 10175.046 3331.693 5378.219
#> h-2 21621.71  6665.539  4548.175  4876.768 1415.294 2809.422
#> h-3 24815.28  7597.319  5422.476  5169.599 1691.818 3334.412
#> h-4 29432.98  9453.527  5950.275  6592.491 2040.414 3892.674
#> ... (8 more rows, 519 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
#> Class structure: foreco > matrix
```

#### A priori constrained (immutable) forecasts

Sometimes we may wish to incorporate a priori knowledge during the
reconciliation process ([Zhang et al., 2023](#ref-Zhang2023-yj)) in
order to improve the accuracy of the reconciled forecasts. For example,
suppose we want to fix the forecasts of the states level series at the
base forecasts values.

``` r

rf_imm <- csrec(base = base, agg_mat = vnaggmat, res = res, comb = "shr", immutable = c(2:8))
summary(rf_imm)
#> ✔ Cross-sectional point forecast reconciliation
#> 
#> ── Method
#> • Function used: `csrec`
#> • Covariance approximation approach: `shr`
#> • Output: (12 x 525) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Forecast horizons (h): 12
#> • Non-negative forecasts (check): `FALSE`
#> 
#> ── Reconciled forecasts
#>        Total         A         B         C        D        E
#> h-1 50402.94 16970.853 12179.676 10128.619 3370.264 5486.059
#> h-2 21299.59  6726.861  4511.988  4768.206 1387.964 2640.546
#> h-3 24397.73  7685.102  5426.145  4932.492 1679.915 3188.778
#> h-4 29663.68  9787.558  6028.565  6497.439 2021.066 3872.930
#> ... (8 more rows, 519 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
#> Class structure: foreco > matrix
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
base_csjb <- csboot(model, B, 12)
reco_csjb <- cssmp(base_csjb, agg_mat = vnaggmat, res = res, nn = "sntz",
                   comb = "shr")
summary(reco_csjb)
#> ✔ Cross-sectional probabilistic forecast reconciliation
#> 
#> ── Method
#> • Function used: `cssmp(csrec)`
#> • Covariance approximation approach: `shr`
#> • Output: distributional object
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Forecast horizons (h): 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Non-negative reconciliation diagnostics
#>      run_time          tol iter
#> 1 0.002269983 1.490116e-08    0
#> 2 0.002269983 1.490116e-08    0
#> 3 0.002269983 1.490116e-08    0
#> 4 0.002269983 1.490116e-08    0
#> 5 0.002269983 1.490116e-08    0
#> ℹ Showing the first 5 rows of the non-negativity diagnostics info matrix.
#> 
#> ── Reconciled forecasts
#> <distribution[4]>
#>         h-1         h-2         h-3         h-4 
#> sample[100] sample[100] sample[100] sample[100] 
#> ... (8 more elements)
#> Use `print(x, n_row)` to see more elements.
#> Class structure: foreco > distribution

# Extracts mean:
str(mean(reco_csjb), give.attr = FALSE)
#>  num [1:12, 1:525] 50024 21686 25301 29973 23541 ...
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
summary(reco_css)
#> ✔ Cross-sectional probabilistic forecast reconciliation
#> 
#> ── Method
#> • Function used: `cssmp(csrec)`
#> • Covariance approximation approach: `shr`
#> • Output: distributional object
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Forecast horizons (h): 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Non-negative reconciliation diagnostics
#>      run_time          tol iter
#> 1 0.001087904 1.490116e-08    0
#> 2 0.001087904 1.490116e-08    0
#> 3 0.001087904 1.490116e-08    0
#> 4 0.001087904 1.490116e-08    0
#> 5 0.001087904 1.490116e-08    0
#> ℹ Showing the first 5 rows of the non-negativity diagnostics info matrix.
#> 
#> ── Reconciled forecasts
#> <distribution[4]>
#>         h-1         h-2         h-3         h-4 
#> sample[100] sample[100] sample[100] sample[100] 
#> ... (8 more elements)
#> Use `print(x, n_row)` to see more elements.
#> Class structure: foreco > distribution
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
summary(rf_opt)
#> ✔ Cross-sectional point forecast reconciliation
#> 
#> ── Method
#> • Function used: `csrec`
#> • Covariance approximation approach: `wls`
#> • Output: (4 x 21) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 21
#> • Forecast horizons (h): 4
#> • Non-negative forecasts (check): `FALSE`
#> 
#> ── Reconciled forecasts
#>          GDP       D1    P3_P5       P3      P5G P31_S14_S15
#> h-1 434329.3 169366.2 422259.7 346156.0 76103.62    265333.7
#> h-2 453374.5 185021.2 432796.4 351654.7 81141.71    269408.4
#> h-3 451278.0 172359.1 429331.7 354245.5 75086.25    276233.6
#> h-4 482268.9 208754.0 464192.2 369154.7 95037.46    271789.6
#> ... (0 more rows, 15 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
#> Class structure: foreco > matrix
```

#### Practical challenge: immutable forecast

In this case, we want to fix the forecasts of the top level series
(\\GDP\\) at the base forecasts values.

``` r

rf_imm <- csrec(base = base, cons_mat = gdpconsmat, res = res, comb = "wls", immutable = c(1))
summary(rf_imm)
#> ✔ Cross-sectional point forecast reconciliation
#> 
#> ── Method
#> • Function used: `csrec`
#> • Covariance approximation approach: `wls`
#> • Output: (4 x 21) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 21
#> • Forecast horizons (h): 4
#> • Non-negative forecasts (check): `FALSE`
#> 
#> ── Reconciled forecasts
#>          GDP       D1    P3_P5       P3      P5G P31_S14_S15
#> h-1 435117.1 169424.0 422702.7 346277.9 76424.76    265406.0
#> h-2 454372.1 185094.4 433357.3 351809.0 81548.34    269499.9
#> h-3 451935.2 172407.3 429701.2 354347.1 75354.13    276293.9
#> h-4 483302.0 208829.7 464773.1 369314.5 95458.58    271884.3
#> ... (0 more rows, 15 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
#> Class structure: foreco > matrix
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
base_csjb <- csboot(model, B, 4)
reco_csjb <- cssmp(base_csjb, cons_mat = gdpconsmat, res = res, comb = "shr")
summary(reco_csjb)
#> ✔ Cross-sectional probabilistic forecast reconciliation
#> 
#> ── Method
#> • Function used: `cssmp(csrec)`
#> • Covariance approximation approach: `shr`
#> • Output: distributional object
#> 
#> ── Structure
#> • Number of cross-sectional series: 21
#> • Forecast horizons (h): 4
#> • Non-negative forecasts (check): `FALSE`
#> 
#> ── Reconciled forecasts
#> <distribution[4]>
#>         h-1         h-2         h-3         h-4 
#> sample[100] sample[100] sample[100] sample[100] 
#> Class structure: foreco > distribution

# Extracts mean:
str(mean(reco_csjb), give.attr = FALSE)
#>  num [1:4, 1:21] 434228 452963 451142 481849 169111 ...
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
summary(reco_css)
#> ✔ Cross-sectional probabilistic forecast reconciliation
#> 
#> ── Method
#> • Function used: `cssmp(csrec)`
#> • Covariance approximation approach: `shr`
#> • Output: distributional object
#> 
#> ── Structure
#> • Number of cross-sectional series: 21
#> • Forecast horizons (h): 4
#> • Non-negative forecasts (check): `FALSE`
#> 
#> ── Reconciled forecasts
#> <distribution[4]>
#>         h-1         h-2         h-3         h-4 
#> sample[100] sample[100] sample[100] sample[100] 
#> Class structure: foreco > distribution
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
