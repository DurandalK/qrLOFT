\name{Companies}
\alias{Companies}
\docType{data}
\title{
Forbes 500 Companies Dataset
}
\description{
Facts about companies selected from the Forbes 500 list
for 1986.
}
\usage{data("Companies")}
\format{
  This dataset contains 79 observations with 9 variables.
  \describe{
    \item{\code{Company}}{Company name.}
    \item{\code{Assets}}{Companies' assets.}
    \item{\code{Sales}}{Companies' sales.}
    \item{\code{Market.Value}}{Companies' estimated market value.}
    \item{\code{Profits}}{Companies' profits.}
    \item{\code{Cash.Flow}}{Companies' cash flow.}
    \item{\code{Employees}}{Numbers of companies' employees.}
    \item{\code{sector}}{Companies' industry sector.}
    \item{\code{Banks}}{Whether this company is a bank or not.}
  }
}
\details{
Facts about companies selected from the Forbes 500 list for 1986. This is a 1/10 systematic sample from the alphabetical list of companies. The Forbes 500 includes all companies in the top 500 on any of the criteria, and thus has almost 800 companies in the list. More details about this dataset can be found on The Data And Story Library, DASL.
}
\source{
\url{https://dasl.datadescription.com/datafile/companies/}
}
\references{
Forbes magazine
}
\examples{
data(Companies)
}
\keyword{datasets}
