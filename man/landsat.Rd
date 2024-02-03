\name{landsat}
\alias{landsat}
\docType{data}
\title{LANDSAT Data: Prediction of County Crop Areas Using Survey and
    Satellite Data}
\description{
    The \code{landsat} data is a compilation of survey and satellite
    data from Battese et al. (1988). It consists of data on segments
    (primary sampling unit; 1 segement approx. 250 hectares) under
    corn and soybeans for 12 counties in north-central Iowa.
}
\usage{data(landsat)}
\format{
    A data frame with 37 observations on the following 10 variables.
    \describe{
        \item{\code{SegmentsInCounty}}{number of segments per county.}
        \item{\code{SegementID}}{sample segment identifier (per county).}
        \item{\code{HACorn}}{hectares of corn for each sample segment
            (as reported in the June 1978 Enumerative Survey).}
        \item{\code{HASoybeans}}{hectares of soybeans for each sample
            segment (as reported in the June 1978 Enumerative Survey).}
        \item{\code{PixelsCorn}}{number of pixels classified
            as corn for each sample segment (LANDSAT readings).}
        \item{\code{PixelsSoybeans}}{number of pixels classified as
            soybeans for each sample segment (LANDSAT readings).}
        \item{\code{MeanPixelsCorn}}{county-specific mean number of
            pixels classified as corn.}
        \item{\code{MeanPixelsSoybeans}}{county-specific mean number
            of pixels classified as soybeans.}
        \item{\code{outlier}}{outlier indicator; observation number 33
           is flagged as outlier.}
        \item{\code{CountyName}}{county names (factor variable):
            Cerro Gordo, Hamilton, Worth, Humboldt, Franklin, Pocahontas,
            Winnebago, Wright, Webster, Hancock, Kossuth, Hardin.}
    }
}
\details{
    The \code{landsat} data in Battese et al. (1988) is a compilation
    of the LANDSAT satellite data from the U.S. Department of Agriculture
    (USDA) and the 1978 June Enumerative Survey.

    \subsection{Survey data}{
        The survey data on the areas under corn and soybeans (reported in
        hectares) in the 37 segments of the 12 counties (north-central Iowa)
        have been determined by USDA Statistical Reporting Service staff,
        who interviewed farm operators. A segment is about 250 hectares.
    }

    \subsection{Satellite data}{
        For the LANDSAT satellite data, information is recorded as "pixels".
        A pixel is about 0.45 hectares. The USDA has been engaged in
        research toward transforming satellite information into good
        estimates of crop areas at the individual pixel and segments
        level. The satellite (LANDSAT) readings were obtained during
        August and September 1978.
    }

    Data for more than one sample segment are available for several counties
    (i.e., unbalanced data).

    Observations No. 33 has been flaged as outlier; see Battese et al.
    (1988, p. 28).
}
\source{
    The \code{landsat} data is from Table 1 of Battese et al. (1988, p. 29).
}
\references{
Battese, G. E., Harter, R. M., and W.A. Fuller (1988).
    An error component model for prediction of county crop areas using.
    \emph{Journal of the American Statistical Association} \bold{83}, 28--36.
    \doi{https://doi.org/10.2307/2288915}
}
\seealso{
    \code{\link{landsat_means}}
}
\examples{
head(landsat)
}
\keyword{datasets}
