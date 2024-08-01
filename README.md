
# FoReco <img src="man/figures/logo.svg" alt="logo" align="right" width="150" style="border: none; float: right;"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/daniGiro/FoReco/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/daniGiro/FoReco/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/FoReco)](https://CRAN.R-project.org/package=FoReco)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![devel
version](https://img.shields.io/badge/devel%20version-1.0.0-blue.svg)](https://github.com/daniGiro/FoReco)
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-forestgreen.svg)](https://cran.r-project.org/web/licenses/GPL-3)

<!-- badges: end -->

**Fo**recast **Reco**nciliation is a a post-forecasting process aimed to
improve the accuracy and align forecasts for a system of linearly
constrained (e.g. hierarchical/grouped) time series. The **FoReco**
package provides a comprehensive set of classical (bottom-up, top-down
and middle-out), and modern (optimal and heuristic combination) forecast
reconciliation procedures in different frameworks including
cross-sectional, temporal, or cross-temporal settings.

The core functions for reconciliation categorized by framework are as
follows:

| **Reconciliation**                |                         **Cross-sectional**                         |                            **Temporal**                             |                         **Cross-Temporal**                          |
|-----------------------------------|:-------------------------------------------------------------------:|:-------------------------------------------------------------------:|:-------------------------------------------------------------------:|
| *Classical reconciliation*        |                                                                     |                                                                     |                                                                     |
| Top-down: `*td()`                 |  [`cstd()`](https://danigiro.github.io/FoReco/reference/cstd.html)  |  [`tetd()`](https://danigiro.github.io/FoReco/reference/tetd.html)  |  [`cttd()`](https://danigiro.github.io/FoReco/reference/cttd.html)  |
| Bottom-up: `*bu()`                |  [`csbu()`](https://danigiro.github.io/FoReco/reference/csbu.html)  |  [`tebu()`](https://danigiro.github.io/FoReco/reference/tebu.html)  |  [`ctbu()`](https://danigiro.github.io/FoReco/reference/ctbu.html)  |
| Middle-out: `*mo()`               |  [`csmo()`](https://danigiro.github.io/FoReco/reference/csmo.html)  |  [`temo()`](https://danigiro.github.io/FoReco/reference/temo.html)  |  [`ctmo()`](https://danigiro.github.io/FoReco/reference/ctmo.html)  |
| *Regression‑based reconciliation* |                                                                     |                                                                     |                                                                     |
| Least squares: `*rec()`           | [`csrec()`](https://danigiro.github.io/FoReco/reference/csrec.html) | [`terec()`](https://danigiro.github.io/FoReco/reference/terec.html) | [`ctrec()`](https://danigiro.github.io/FoReco/reference/ctrec.html) |
| LCC: `*lcc()`                     | [`cslcc()`](https://danigiro.github.io/FoReco/reference/cslcc.html) | [`telcc()`](https://danigiro.github.io/FoReco/reference/telcc.html) | [`ctlcc()`](https://danigiro.github.io/FoReco/reference/ctlcc.html) |

Additionally, **FoReco** provides various functions for different
aspects of forecast reconciliation, including aggregating time series,
balancing hierarchies, computing projection and covariance matrices, and
more.

## Installation

You can install the **stable** version on [R
CRAN](https://cran.r-project.org/)

``` r
install.packages("FoReco")
```

You can also install the **development** version from
[Github](https://github.com/daniGiro/FoReco)

``` r
# install.packages("devtools")
devtools::install_github("danigiro/FoReco")
```

## Getting Started

To get started with using the FoReco package, refer to the
[documentation](https://danigiro.github.io/FoReco/) for detailed
information on how to apply different forecast reconciliation procedures
to your data.

<!-- - [Introduction to `FoReco`](https://danigiro.github.io/FoReco/articles/Introduction-to-FoReco.html) -->

- [The `vndata` and `itagdp`
  dataset](https://danigiro.github.io/FoReco/articles/Dataset-vndata-and-itagdp.html)

- [Cross-sectional forecast
  reconciliation](https://danigiro.github.io/FoReco/articles/Cross-sectional-forecast-reconciliation.html)

- [Temporal forecast
  reconciliation](https://danigiro.github.io/FoReco/articles/Temporal-forecast-reconciliation.html)

- [Cross-temporal forecast
  reconciliation](https://danigiro.github.io/FoReco/articles/Cross-temporal-forecast-reconciliation.html)

- [Replicate the `hts`
  package](https://danigiro.github.io/FoReco/articles/Replicate-the-hts-package.html)

- [Replicate the `thief`
  package](https://danigiro.github.io/FoReco/articles/Replicate-the-thief-package.html)

## Issues and Contributions

If you encounter any bugs or have suggestions for improvements, please
feel free to report them on [GitHub Issues
page](https://github.com/daniGiro/FoReco/issues). Contributions are also
welcome!
