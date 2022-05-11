# CHANGES in rsae VERSION 0.1-6 (2022-05-11)

## BUG FIXES

* fixed namespace issues
* fixed DESCRIPTION file
* replaced `require()` in function `.initmethod()` with
    `requireNamespace()`

## CHANGES

* `coef()` method for objects of class `fitsaemodel` only returns the coefficients but does print to the console anymore 
* `fitsaemodel.control()` gained the additional argument `k_Inf` to specify the robustness tuning constant k that represents infinity (used to define the maximum likelihood estimator)

## MISC

* updated CITATION file
* cleaned documentation and vignette
* combined several R source files into one source file
* source repository moved from R-Forge to GitHub

# CHANGES in rsae VERSION 0.1-5 (2014-02-13)
added vignette folder

# CHANGES in rsae VERSION 0.1-4 (2012-01-08)

## NEW FEATURES

* added status bar for parametric bootstrap
* added arguments `dec` and `decor` in `fitsaemodel.control()`
    to choose between different types of decorrelation of the
    residuals and the type of matrix square roots (Cholesky vs. SVD)

## MISC

* Modified documentation
* implemented better checks for the robustness tuning constant
    in the internal function `.fitsaemodel.huberm()`

## BUG FIXES

* fixed the names of coefficients returned by `makedata()`
* fixed print method `fitsaemodel()` when different robustness
    tuning are used

# CHANGES in rsae VERSION 0.1-3 (2011-07-26)

## NEW FEATURES

* added function `convergence()` to learn more about the convergence of the
    methods
* added argument `areameans` in function `robpredict()`
* added argument `type` in plot method for class `meanssaemodel`


## MISC

removed argument `full` in the summary method for objects of class `fitmodel`

# CHANGES in rsae VERSION 0.1-2 (2011-07-23)

Initial release
