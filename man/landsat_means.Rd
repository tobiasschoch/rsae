\name{landsat_means}
\alias{landsat_means}
\docType{data}
\title{Means of the LANDSAT Data for Corn and Soybeans}
\description{
    The \code{\link{landsat}} data is a compilation of survey and satellite
    data from Battese et al. (1988). The county-specific population means of
    pixels of the segments under corn and soybeans, respectively, are available
    in the data.frame \code{landsat_means}.
}
\usage{data(landsat_means)}
\format{
    A data frame with 12 observations (counties) on the following 3 variables.
    \describe{
        \item{\code{(intercept)}}{all ones.}
        \item{\code{MeanPixelsCorn}}{county-specific mean of number of pixels
            classified as corn (LANDSAT readings).}
        \item{\code{MeanPixelsSoybeans}}{county-specific number of pixels
            classified as soybeans (LANDSAT readings).}
    }
}
\details{
    The data.frame \code{landsat_means} is an aggregation of the data.frame
    \code{\link{landsat}}.
}
\references{
Battese, G. E., Harter, R. M., and W.A. Fuller (1988).
    An error component model for prediction of county crop areas using.
    \emph{Journal of the American Statistical Association} \bold{83}, 28--36.
    \doi{https://doi.org/10.2307/2288915}
}
\examples{
head(landsat_means)
}
\keyword{datasets}
