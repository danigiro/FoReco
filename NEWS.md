# FoReco 0.1.2
**Note**: *Documentation for this version is still under development*

##### Major changes
* Add the possibility to use part of factors of *m* (max. order of temporal aggregation)
* Add the possibility for `htsrec()`, `thfrec()` and `octrec()` to introduce a list of *h* covariance matrices in the parameters `W` and `Omega`, where *h* stands for the forecast horizon (note that for `thfrec()` and `octrec()` this is the forecast horizon of the entire cycle)

##### Minor changes
* Now in `octrec()` it is also possible to introduce the **Ω** covariance matrix variant through the `Omega` parameter and not only the **W** variant with the`W` parameter
* Renew `tcsrec()`, `cstrec()` and `iterec()`. Also in the iterec function the `maxit` parameter has been replaced by `itmax`, however for the moment `maxit` is still supported

##### Experimental
* Add the possibility to introduce constraints through the `bounds` param in `htsrec()`, `thfrec()` and `octrec()`.
* Add a (hidden) function `oct_bounds()` to organize the bounds on a singular dimension (i.e. only cross-sectional or only temporal) in a cross-temporal framework

# FoReco 0.1.1
This is a small release focusing on fixing some bugs and the documentation

* Fixed a bug in `iterec()` when calculating incoherence
* Fixed documentation 
* Changed the contact mail (now it's daniele.girolimetto@phd.unipd.it)
* Corrected the second section of the vignette "`Average relative accuracy indices`"

# FoReco 0.1.0

* Release on github
