\name{GrowthData}
\docType{data}
\alias{GrowthData}
\title{
    Barro & Lee Growth Data Set
}
\description{
    Data set of growth compiled by Barro Lee, directly excerpted from pacakage \code{hdm}.
}
\usage{data("GrowthData")}
\format{
This dataset contains 90 observations with 63 variables.
    \describe{
    \item{\code{Outcome}}{national growth rates in GDP per capita for the periods 1965-1975 and 1975-1985;}
    \item{\code{...}}{covariates which might influence growth.}
    }
}
\details{
    The data set contains growth data of Barro-Lee. The Barro Lee data consists of a panel of 138 countries for the period 1960 to 1985. The dependent variable is national growth rates in GDP per capita for the periods 1965-1975 and 1975-1985. The growth rate in GDP over a period from t_1 to t_2 is commonly defined as \eqn{\log(GDP_{t_1}/GDP_{t_2}}). The number of covariates is p=62. The number of complete observations is 90.
}
\source{
The full data set and further details can be found at
\url{http://www.nber.org/pub/barro.lee}, \url{http://www.barrolee.com}, and,
\url{http://www.bristol.ac.uk//Depts//Economics//Growth//barlee.htm}.
}
\references{
    Barro, R. J. and Lee, J. W. (1994), Data set for a panel of 139
countries. NBER.

    Barro, R. J. and Sala-i-Martin, X. (1995), \emph{Economic Growth}, McGrwa-Hill, New York.
}
\note{
This dataset is directly excerpted from package \code{hdm}, for more information please check the webpages listed below.
}
\examples{
data(GrowthData)
}
\keyword{datasets}
