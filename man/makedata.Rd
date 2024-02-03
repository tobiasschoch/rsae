\name{makedata}
\alias{makedata}
\title{Synthetic Data Generation for the Basic Unit-Level SAE Model}
\description{
    This function generates synthetic data (possibly contaminated by outliers)
    for the basic unit-level SAE model.
}
\usage{
makedata(seed = 1024, intercept = 1, beta = 1, n = 4, g = 20, areaID = NULL,
         ve = 1, ve.contam = 41, ve.epsilon = 0, vu = 1, vu.contam = 41,
         vu.epsilon = 0)
}
\arguments{
    \item{seed}{\code{[integer]} seed value used in \code{set.seed} (default
        \code{seed = 1024}).}
    \item{intercept}{\code{[numeric]} or \code{[NULL]} value of the intercept
        of the fixed-effects model or \code{NULL} for a model without
        intercept (default: \code{intercept = 1}).}
    \item{beta}{\code{[numeric vector]} value of the fixed-effect coefficients
        (without intercept; default: \code{beta = 1}). For each given
        coefficient, a vector of realizations is drawn from the standard
        normal distribution.}
    \item{n}{\code{[integer]} number of units per area in balanced-data
        setups (default: \code{n = 4}).}
    \item{g}{\code{[integer]} number of areas (default: \code{g = 20}).}
    \item{areaID}{\code{[integer vector]} or \code{[NULL]}. If one attempts
        to generate synthetic unbalanced data, one calls \code{makedata} with
        a vector, the elements of which area identifiers. This vector
        should contain a series of (integer valued) area IDs. The number
        of areas is set equal to the number unique IDs.}
    \item{ve}{\code{[numeric]} nonnegative value of model/ residual variance.}
    \item{ve.contam}{\code{[numeric]} nonnegative value of model variance of
        the outlier part in a mixture distribution (Tukey-Huber-type
        contamination model) \eqn{e = (1-h)N(0, ve) + hN(0, ve.contam)
            }{e = (1-h)*N(0, ve) + h*N(0, ve.contam)}.}
    \item{ve.epsilon}{\code{[numeric]} value in \eqn{[0,1]}{[0,1]} that
        defines the relative number of outliers (i.e., epsilon or h in
        the contamination mixture distribution). Typically, it takes
        values between 0 and 0.5 (but it is not restricted to this interval).}
    \item{vu}{\code{[numeric]} value of the (area-level) random-effect
        variance.}
    \item{vu.contam}{\code{[numeric]} nonnegative value of the (area-level)
        random-effect variance of the outlier part in the contamination
        mixture distribution.}
    \item{vu.epsilon}{\code{[numeric]} value in \eqn{[0,1]}{[0,1]} that
        defines the relative number of outliers in the contamination
        mixture distribution of the (area-level) random effects.}
}
\details{
    Let \eqn{y_i}{e[i]} denote an area-specific \eqn{n_i}{n[i]}-vector of
    the response variable for the areas \eqn{i = 1,..., g}. Define a
    \eqn{(n_i \times p)}{(n[i] * p)}-matrix \eqn{X_i}{X[i]} of realizations
    from the std. normal distribution, \eqn{N(0,1)}{N(0,1)}, and let
    \eqn{\beta} denote a \eqn{p}-vector of regression coefficients. Now, the
    \eqn{y_i}{y[i]} are drawn using the law
    \eqn{y_i \sim N(X_i\beta, v_e I_i + v_u J_i)}{y[i] ~ N(X[i]\beta,
    v[e] I[i] + v[u] J[i])} with \eqn{v_e}{v[e]} and
    \eqn{v_u}{v[u]} the variances of the model error and random-effect
    variance, respectively, and \eqn{I_i}{I[i]} and \eqn{J_i}{J[i]} denoting
    the identity matrix and matrix of ones, respectively.

    In addition, we allow the distribution of the model/residual and
    area-level random effect to be contaminated (cf. Stahel and Welsh, 1997).
    Notably, the laws of \eqn{e_{i,j}}{e[ij]} and \eqn{u_i}{u[i]} are replaced
    by the Tukey-Huber contamination mixture:
    \itemize{
        \item \eqn{e_{i,j} \sim (1-\epsilon^{ve})N(0,v_e) +
            \epsilon^{ve}N(0, v_e^{\epsilon})}{e[i,j] ~ (1 -
            \epsilon[ve]) N(0,v[e]) + \epsilon[ve] N(0, v[e,\epsilon])}

        \item \eqn{u_{i} \sim (1-\epsilon^{vu})N(0,v_u) +
            \epsilon^{vu}N(0, v_u^{\epsilon})}{u[i] ~ (1 -
            \epsilon[vu]) N(0,v[u]) + \epsilon[vu] N(0, v[u,\epsilon])}
    }
    where \eqn{\epsilon^{ve}}{\epsilon[ve]} and
    \eqn{\epsilon^{vu}}{\epsilon[vu]} regulate the degree of contamination;
    \eqn{v_e^{\epsilon}}{v[e,\epsilon]} and \eqn{v_u^{\epsilon}}{v[u,\epsilon]}
    define the variance of the contamination part of the mixture distribution.

    Four different contamination setups are possible:
        \itemize{
            \item no contamination (i.e., \code{ve.epsilon = vu.epsilon = 0}),
            \item contaminated model error (i.e., \code{ve.epsilon != 0} and
                \code{vu.epsilon = 0}),
            \item contaminated random effect (i.e., \code{ve.epsilon = 0} and
                \code{vu.epsilon != 0}),
            \item both are conaminated (i.e., \code{ve.epsilon != 0} and
                \code{vu.epsilon != 0}).
        }
}
\value{
    An instance of the class \code{saemodel}.
}
\references{
Schoch, T. (2012). Robust Unit-Level Small Area Estimation: A Fast Algorithm
    for Large Datasets. \emph{Austrian Journal of Statistics} \bold{41},
    243--265. \doi{https://doi.org/10.17713/ajs.v41i4.1548}

Stahel, W. A. and A. Welsh (1997).
    Approaches to robust estimation in the simplest variance components model.
    \emph{Journal of Statistical Planning and Inference} \bold{57}, 295--319.
    \doi{https://doi.org/10.1016/S0378-3758(96)00050-X}
}
\seealso{
    \code{\link[=saemodel]{saemodel()}},
    \code{\link[=fitsaemodel]{fitsaemodel()}}
}
\examples{
# generate a model with synthetic data
model <- makedata()
model

# summary of the model
summary(model)
}
