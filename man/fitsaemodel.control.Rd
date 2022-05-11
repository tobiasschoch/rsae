\name{fitsaemodel.control}
\alias{fitsaemodel.control}
\title{
Tuning Parameters of \code{fitsaemodel}
}
\description{
This function is used to define global settings and parameters that
are used by \code{fitsaemodel}.
}
\usage{
fitsaemodel.control(niter = 40, iter = c(200, 200), acc = 1e-05,
    dec = 0, decorr = 0, init = "default", k_Inf = 20000, ...)
}
\arguments{
    \item{niter}{integer that defines the maximum number of outer-loop
        iterations (default: \code{niter = 40}).}
    \item{iter}{integer or vector of integers of size 2 that defines the
        maximum number of inner-loop iterations. The first element of the
        vector \code{iter} refers to the estimation of the regression
        coefficients \eqn{\beta}{beta}; the second element refers to the
        estimation of the variance of the unit-level errors, \eqn{v}; the
        maximum number of iterations used to compute the ratio of variances,
        \eqn{d}, cannot be modified (default: \code{iter = c(200, 200)}).}
    \item{acc}{scalar or numeric vector of size 4 that define the numeric
        tolerance used in the termination rule of the iterative updating
        algorithms. The positions \code{1:4} of the vector \code{acc}
        refer to \code{1}: (overall) outer-loop, \code{2}: regression
        coefficients, \eqn{\beta}{beta}, \code{3}: variance component,
        \eqn{v}, \code{4}: ratio of variances \eqn{d}; default:
        \code{acc = 1e-05}.}
    \item{dec}{integer to toggle the type of the matrix square root
        (decomposition); \code{dec = 0} for SVD (default) or \code{dec = 1}
        for Cholesky decomposition.}
    \item{decorr}{integer to toggle the type of decorrelation of the
        residuals (default \code{decorr = 0}).}
    \item{init}{a character string that specifies by what method the main
        algorithm is initialized. By default: \code{init = "default"};
        alternatively, (and provided that the \pkg{robustbase} package is
        installed) one may choose a high-breakdown-point initial estimate:
        either "lts" (fast LTS regression) or "s" (S-estimate of regression).
        For datasets with more than 100,000 observations, the former is
        rather slow. For more details on the initializing methods, see the
        documentation of \pkg{robustbase} (\code{"ltsReg"} and
        \code{"lmrob.S"}).}
    \item{k_Inf}{tuning constant that represents infinity (default:
        \code{k_Inf = 20000}).}
    \item{\dots}{(may be used in the future).}
}
\details{
Changing the default values of the parameters may result in failure of
convergence or loss of convergence speed.
}
\value{
(an object used by the robust methods)
}
\seealso{
\code{\link{fitsaemodel}}
}
