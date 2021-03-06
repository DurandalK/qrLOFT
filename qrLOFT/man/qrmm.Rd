\name{qrmm}
\alias{qrmm}
\title{
    Funciton to estimate quantile regression coefficients
}
\description{
    Funciton to estimate quantile regression coefficients using MM-algorithm.
}
\usage{
qrmm(X, Y, tau)
}
\arguments{
    \item{X}{the design matrix}
    \item{Y}{the response variable}
    \item{tau}{current quantile level, should lie in (0, 1).}
}
\details{
    This function can only be used to estimate multivariate quantile regression, otherwise a singluar design matrix error might be returned.
}
\value{
    \item{beta_out}{quantile regression coefficients estimates.}
}
\references{
    Hunter, D. R. and Lange, K. (2000), Quantile Regression via an MM Algorithm, \emph{Journal of Computational and Graphical Statistics}, \bold{9}, 60-77.
}
\examples{
data("Companies")
qrmm(as.matrix(Companies[, c(5,6,7)]), Companies[, 2], 0.25)
}
\keyword{MM algorithm}
