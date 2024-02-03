\name{saemodel}
\alias{saemodel}
\alias{print.saemodel}
\alias{summary.saemodel}
\alias{as.matrix.saemodel}
\title{Setting Up a SAE Model}
\description{
    Function \code{saemodel()} is used to specify a model. Once a model
    has been specified, it can be fitted using
    \code{\link[=fitsaemodel]{fitsaemodel()}} by different estimation methods.
}
\usage{
saemodel(formula, area, data, type = "b", na.omit = FALSE)

\method{print}{saemodel}(x, ...)
\method{summary}{saemodel}(object, ...)
\method{as.matrix}{saemodel}(x, ...)
}
\arguments{
    \item{formula}{a \code{formula} object of describing the fixed-effects
        part of the model, with the response on the RHS of the \code{~}
        operator and the terms or regressors, separated by \code{+}
        operators, on the LHS of the formula.}
    \item{area}{a one-sided \code{formula} object. A \code{~} operator
        followed by only one single term defining the area-specific
        random-effect part.}
    \item{data}{data.frame.}
    \item{type}{\code{[character]} \code{"a"} or \code{"b"} refering to
        J.N.K. Rao's definition of model type A (area-level model) or B
        (unit-level model); default is \code{"b"}.}
    \item{na.omit}{\code{[logical]} indicating whether \code{NA} values
        should be removed before the computation proceeds. Note that none
        of the algorithms can cope with missing values.}
    \item{x}{an object of class \code{"saemodel"}.}
    \item{object}{an object of the class \code{"saemodel"}.}
    \item{\dots}{additional arguments (not used).}
}
\details{
    Function \code{saemodel()} is used to specify a model.
    \itemize{
        \item \code{model} is a symbolic description (\code{formula} of the
            fixed-effects model to be fitted.

            A typical model has the form \code{response ~ terms} where
            \code{response} is the (numeric) response vector and
            \code{terms} is a series of terms which specifies a linear
            predictor for response (explanatory variables); see
            \code{\link{formula}}.

            A \code{formula} has an implied intercept term. To remove
            this use either \code{y ~ x - 1} or \code{y ~ 0 + x};
            see \code{\link{formula}} for more details of allowed formulae.
        \item \code{area} is a symbolic description (\code{formula}) of
            the random effects (nested error structure). It must be
            right-hand side only formula consisting of one term,
            e.g., \code{~ areaDefinition}.
    }
    The data must no contain missing values.

    The design matrix (i.e., matrix of the explanatory variables
    defined the right-hand side of \code{model}) must have full column
    rank; otherwise execution is terminated by an error.

    Once a model has been specified, it can be fitted by
    \code{\link[=fitsaemodel]{fitsaemodel()}}.
}
\value{
    An instance of the S3 class \code{"saemodel"}
}
\references{
Rao, J.N.K. (2003). \emph{Small Area Estimation}, New York: John Wiley and Sons.
}
\seealso{
    \code{\link[=makedata]{makedata()}},
    \code{\link[=fitsaemodel]{fitsaemodel()}}
}
\examples{
# use the landsat data
head(landsat)

# set up the model
model <- saemodel(formula = HACorn ~ PixelsCorn + PixelsSoybeans,
    area = ~CountyName,
    data = subset(landsat, subset = (outlier == FALSE)))

# summar of the model
summary(model)
}
