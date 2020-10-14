
<!-- README.md is generated from README.Rmd. Please edit that file -->

FoReco <img src="man/figures/logo.svg" align="right" alt="logo" width="150" style = "border: none; float: right;">
==================================================================================================================

<!-- badges: start -->

[![R build
status](https://github.com/daniGiro/FoReco/workflows/R-CMD-check/badge.svg)](https://github.com/daniGiro/FoReco/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/FoReco)](https://CRAN.R-project.org/package=FoReco)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![devel
version](https://img.shields.io/badge/devel%20version-0.1.1-blue.svg)](https://github.com/daniGiro/FoReco)
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-forestgreen.svg)](https://cran.r-project.org/web/licenses/GPL-3)
<!-- badges: end -->

The **FoReco** (**Fo**recast **Reco**nciliation) package is designed for
point forecast reconciliation, a **post-forecasting** process aimed to
improve the quality of the base forecasts for a system of linearly
constrained (e.g. hierarchical/grouped) time series.

It offers classical (bottom-up), optimal and heuristic combination
forecast reconciliation procedures by exploiting cross-sectional,
temporal, and cross-temporal relationships linking the time series.

The main functions are:

-   `htsrec()`: cross-sectional (contemporaneous) forecast
    reconciliation.
-   `thfrec()`: forecast reconciliation for a single time series through
    temporal hierarchies.
-   `tcsrec()`: heuristic first-temporal-then-cross-sectional
    cross-temporal forecast reconciliation.
-   `cstrec()`: heuristic first-cross-sectional-then-temporal
    cross-temporal forecast reconciliation.
-   `iterec()`: heuristic iterative cross-temporal forecast
    reconciliation.
-   `octrec()`: optimal combination cross-temporal forecast
    reconciliation.

Installation
------------

You can install the **stable** version on [R
CRAN](https://cran.r-project.org/)

    install.packages('FoReco', dependencies = TRUE)

You can also install the **development** version from
[Github](https://github.com/daniGiro/FoReco)

    # install.packages("devtools")
    devtools::install_github("daniGiro/FoReco", build_vignettes = TRUE)

Links
-----

-   Source code:
    <a href="https://github.com/daniGiro/FoReco" class="uri">https://github.com/daniGiro/FoReco</a>
-   Site documentation:
    <a href="https://danigiro.github.io/FoReco/" class="uri">https://danigiro.github.io/FoReco/</a>

Getting help
------------

If you encounter a clear bug, please file a minimal reproducible example
on [GitHub](https://github.com/daniGiro/FoReco/issues).
