\name{fitsaemodel.control}
\alias{fitsaemodel.control}
\title{Tuning Parameters of \code{fitsaemodel}}
\description{
    This function is used to define global settings and parameters that
    are used by \code{\link[=fitsaemodel]{fitsaemodel()}}.
}
\usage{
fitsaemodel.control(niter = 40, iter = c(200, 200), acc = 1e-05,
    dec = 0, decorr = 0, init = "default", k_Inf = 20000, ...)
}
\arguments{
    \item{niter}{\code{[integer]} the maximum number of outer-loop
        iterations (default: \code{niter = 40}).}
    \item{iter}{\code{[integer]} the maximum number of inner-loop
        iterations. It can be a vector of size 2. The first element of the
        vector refers to the estimation of the regression
        coefficients \eqn{\beta}{beta}; the second element refers to the
        estimation of the variance of the unit-level errors, \eqn{v}; the
        maximum number of iterations used to compute the ratio of variances,
        \eqn{d}, cannot be modified (default: \code{iter = c(200, 200)}).}
    \item{acc}{\code{[numeric]} numeric tolerance used in the
        termination rule of the iterative updating algorithms.
        It can be a vector of size 4. The positions \code{1:4} of the
        vector \code{acc} refer to \code{1}: (overall) outer-loop,
        \code{2}: regression coefficients, \eqn{\beta}{beta}, \code{3}:
        variance component, \eqn{v}, \code{4}: ratio of variances
        \eqn{d}; default: \code{acc = 1e-05}.}
    \item{dec}{\code{[integer]} type of matrix square root (decomposition);
        \code{dec = 0} for eigenvalue decomposition (default)
        or \code{dec = 1} for Cholesky decomposition.}
    \item{decorr}{\code{[integer]} type of decorrelation of the
        residuals; \code{decorr = 0}: no robust decorrelation (default);
        \code{decorr = 1}: means are replaced by medians.}
    \item{init}{\code{[character]} method by which the main
        algorithm is initialized. By default, \code{init = "default"}
        the algorithm is initialized by a robust fixed-effects estimator;
        alternatively, (and provided that the \pkg{robustbase} package is
        installed) one may choose one of the high-breakdown-point initial
        estimators: "lts" (fast least-trimmed squares, LTS, regression) or
        "s" (regression S-estimator). For more details
        on the initialization methods, see documentation of
        \code{\link[=fitsaemodel]{fitsaemodel()}}.}
    \item{k_Inf}{\code{[numeric]} tuning constant of the robust estimator
        that represents infinity (default: \code{k_Inf = 20000}).}
    \item{\dots}{additional arguments (not used).}
}
\details{
    Changing the default values of the parameters may result in failure of
    convergence or loss of convergence speed.
}
\value{
    A list with entries
    \itemize{
        \item \code{niter}
        \item \code{iter}
        \item \code{acc}
        \item \code{k_Inf}
        \item \code{init}
        \item \code{dec}
        \item \code{decorr}
        \item \code{add}
    }
}
\seealso{
    \code{\link[=fitsaemodel]{fitsaemodel()}}
}
\examples{
# use the landsat data
head(landsat)

# define the saemodel using the landsat data
model <- saemodel(formula = HACorn ~ PixelsCorn + PixelsSoybeans,
    area = ~CountyName,
    data = subset(landsat, subset = (outlier == FALSE)))

# summary of the model
summary(model)

# obtain the maximum likelihood estimates with, for instance, 'niter = 50'
# number of outer-loop iterations (by default: niter = 40). Here, we use
# 'niter = 50' for the sake of demonstration, not because it is needed.
fitsaemodel("ml", model, niter = 50)
}
