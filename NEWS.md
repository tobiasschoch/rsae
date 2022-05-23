# CHANGES in rsae VERSION 0.2 (2022-05-23)

## BUG FIXES

* fixed name space issues
* fixed DESCRIPTION file
* The call of `require()` in function `.initmethod()` has been replaced with
    `requireNamespace()` in order to load functions from the `robustbase` package.

## CHANGES

* Function `fitsaemodel()`
  * It returns now an object of class`fit_model_b` (in earlier versions, the class was called `fitsaemodel`).
  * If the algorithm does not converge, `fitsaemodel()` returns an object with the slots `beta` (estimated fixed-effects coefficients) and `theta` (variance components) equal to a vector of `NA`'s.
  * The covariance matrix of the fixed effects (`vcovbeta`) is now computed in function `summary()` and not `fitsaemodel()`. This goes unnoticed by the user. 
* The  function `coef.fit_model_b()` only returns the coefficients but does not print to the console anymore. 
* Function `fitsaemodel.control()` gains the additional argument `k_Inf` to specify the robustness tuning constant k that represents infinity (used to define the maximum likelihood estimator)
* Argument `digit` of the function `print.fit_model_b()` depends now on `getOption("digits")`.
* Function `summary.fit_model_b()` returns now the model summary.
* Function `robpredict()` 
  * The function gains the additional argument `seed`  to specify the random seed used in the parametric bootstrap.
  * The effective number of bootstrap replicates is now returned and printed to the console. In earlier package versions, the nominal number of replicates was reported, which did not account for bootstrap model estimates that failed to converge. 
  * The progress bar in function `robpredict()` is now implemented with `utils:txtProgressBar()` and can be switched off using  `progress_bar = FALSE`.
* The column of ones of the matrix of independent variables (if there is an intercept), see slot `X` of an `saemodel` instance, is now called `(Intercept)` not `X.Intercept` .


## MISC

* updated CITATION file
* cleaned documentation and vignette
* combined several R source files into one source file and added `init.c`
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
