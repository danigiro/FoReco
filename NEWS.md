# FoReco 0.2.6

##### Major changes: [probabilistic forecast reconciliation](https://danigiro.github.io/FoReco/articles/FoReco_prob.html)
* Added `boot_cs()`, `boot_te()` and `boot_ct()` to draw samples from, respectively, cross-sectional, temporal and cross-temporal joint (block) bootstrap.

##### Minor changes
* Fixed deprecation warnings in [Matrix](https://CRAN.R-project.org/package=Matrix) (v. 1.5-0);
* Improved docs and bug fixes;
* Fixed `ctbu()` inputs;
* Added `FoReco2matrix()` to trasform FoReco forecasts input and output in a list of matrix/vector class;
* Added `agg_ts()`: non-overlapping temporal aggregation of a time series according to a specific aggregation order.
* Added `arrange_hres()` and `residuals_matrix()` functions to arrange residuals for the covariance matrix under Gaussianity;

# FoReco 0.2.5

##### Minor changes
* Fixed the not negative reconciliation "**sntz**" for `octrec()`;
* Fixed documentation.

# FoReco 0.2.4

##### Major changes
* Added `lcmat()` function.

##### Minor changes
* Fixed BU approach when the number of columns of basef is equal to the number of bottom time series `htsrec()`;
* Fixed `score_index()`;
* Fixed the `bounds` param when `type = "S"` in `htsrec()`, `thfrec()` and `octrec()`;
* Add the possibility to fix base forecasts through the `v` param in `htsrec()`, `thfrec()` and `octrec()` - experimental;
* Add two new type of optimal cross-temporal reconciliation (**cs_struc** and **t_struc**);
* Improved docs and bug fixes.

# FoReco 0.2.2

##### Minor changes
* Fixed documentation;
* Removed `ut2c()` and `srref()`.

# FoReco 0.2.1

##### Minor changes
* Fixed a bug in the output of `lccrec()` (now the function returns the Level Conditional Coherent or the Combined Conditional Coherent forecasts);
* Fixed the not negative reconciliation "**KAnn**" when `keep = "list"`.

# FoReco 0.2.0

##### Major changes
* It's possible to use a subset of factors of *m* (max. order of temporal aggregation);
* Added the possibility for `htsrec()`, `thfrec()` and `octrec()` to introduce a list of *h* covariance matrices in the parameters `W` and `Omega`, where *h* stands for the forecast horizon (note that for `thfrec()` and `octrec()` this is the forecast horizon of the entire cycle);
* Param `Sstruc` is no more avaible in `octrec()` and `ctf_tools()`. FoReco uses a fast algorithm to compute **Scheck**, so no external input is needed;
* Modified output of `ctf_tools()` (added `Ccheck`, `Htcheck`, `Scheck`, removed `Cstruc`, `Sstruc`), `hts_tools()` (added `C`) and `thf_tools()` (added `m`);
* Added two new not negative reconciliation techniques ("**KAnn**" and "**sntz**") with a new parameter (`nn_type`) in `htsrec()`, `thfrec()` and `octrec()`;
* Added the top-down reconciliation function `tdrec()`;
* Added the level conditional forecast reconciliation (with and without not-negative constraints) for genuine hierarchical/grouped time series `levrec()` (cross-sectional, temporal and cross-temporal).

##### Minor changes
* Now in `octrec()` it is also possible to introduce the **Ω** covariance matrix variant through the `Omega` parameter and not only the **W** variant with the `W` parameter;
* Updated `tcsrec()`, `cstrec()` and `iterec()`. In the `iterec()` function the `maxit` parameter has been replaced by `itmax`, however for the moment `maxit` is still supported;
* Now FoReco removes null rows from the cross-sectional aggregation matrix **C** and it warns the user if the balanced version of an unbalanced hierarchy is considering duplicated variables;
* Redesigned the console output and added a new convergence norm as *default* for `iterec()` (`norm` parameter).

##### Experimental
* Add the possibility to introduce constraints through the `bounds` param in `htsrec()`, `thfrec()` and `octrec()`;
* Add a function `oct_bounds()` to organize the bounds on a specific dimension (i.e. only cross-sectional or only temporal) in a cross-temporal framework;
* Added `ut2c()` and `srref()` to develop a cross-sectional structural representation starting from a zero constraints kernel matrix;
* Added in `score_index()` the calculation of multiple forecast horizons index (like 1:6) and multiple cross-sectional levels for a forecasting experiment.


# FoReco 0.1.1
Minore release, fixing some bugs and the documentation

* Fixed a bug in `iterec()` when calculating the incoherence
* Fixed documentation 
* Changed the contact mail (now it's daniele.girolimetto@phd.unipd.it)
* Corrected the second section of the vignette "`Average relative accuracy indices`"

# FoReco 0.1.0

* Release on github
