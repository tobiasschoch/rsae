\name{saemodel}
\alias{saemodel}
\alias{print.saemodel}
\alias{summary.saemodel}
\alias{as.matrix.saemodel}
\title{
Setting Up a SAE Model
}
\description{
\code{saemodel} is the workhorse function to set up a model (i.e., an
instance of the \code{"saemodel"} class). It is the starting point of
every model fitting exercise. Once a model has been initilized/ set up,
we consider estimating its parameter.
}
\usage{
saemodel(formula, area, data, type = "b", na.omit = FALSE)

\method{print}{saemodel}(x, ...)
\method{summary}{saemodel}(object, ...)
\method{as.matrix}{saemodel}(x, ...)
}
\arguments{
    \item{formula}{a two-sided linear formula object describing the
        fixed-effects part, with the response on the RHS of the \code{~}
        operator and the terms or regressors, separated by \code{+}
        operators, on the LHS of the formula.}
    \item{area}{a one-sided formula object. A \code{~} operator
        followed by only one single term defining the area-specific
        random-effect part.}
    \item{data}{data.frame.}
    \item{type}{either "a" or "b" refering to J.N.K. Rao's definition
        of model type A (area-level model) or B (unit-level model);
        default is \code{type="b"}.}
    \item{na.omit}{a logical indicating whether \code{NA} should be removed
        (default is \code{FALSE}). Note that none of the algorithms
        can cope with missing values.}
    \item{x}{an object of the class \code{"saemodel"} (this argument is
        implicitly used by the \code{print} and \code{as.matrix} methods).}
    \item{object}{an object of the class \code{"saemodel"} (this argument is
        implicitly used by the \code{summary} method).}
    \item{\dots}{not used.}
}
\details{
The step of setting up a SAE model is the starting point of any (robust)
SAE modeling exercise. (Use the makedata to generate a synthetic
dataset; see also, below). Here, we have to define the fixed-effects-
and random-effects part of the model, and to tell \code{R} what data
it shall use.

Once a model has been initilized/ set up, we consider estimating its
parameter; see \code{\link{fitsaemodel}}.
}
\value{
Instance of the S3 class \code{"saemodel"}.
}
\references{
Rao, J.N.K. (2003).
    \emph{Small Area Estimation}, New York: John Wiley and Sons.
}
\seealso{
\code{\link{makedata}}
}
