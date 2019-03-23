\name{qrLOFT-package}
\alias{qrLOFT-package}
\alias{qrLOFT}
\docType{package}
\title{
  Lack of Fit Test for Multivariate and High Dimensional Quantile Regression Model
}
\description{
    This package provides a fast way to implement the lack of fit test method for both multivariate and high dimensional quantile regression model. More details about the test method can be found in Dong, Li and Feng (2019).
}
\details{
    Package \code{qrLOFT} constains several functions for estimating and testing quantile regression models.
    \describe{
    \item{qrmm}{MM-algorithm for estimating multivariate quantile regression.}
    \item{qrfit}{QICD algorithm for estimating high-dimensional quantile regression.}
    \item{qrloft}{Lack of fit test for multivariate or high dimensional quantile regression.}
    }
}
\author{
Chen Dong, Guodong Li, and Xingdong Feng
}
\references{
    Hunter, D. R. and Lange, K. (2000), Quantile Regression via an MM Algorithm, \emph{Journal of Computational and Graphical Statistics}, \bold{9}, 60-77.
    
    Peng, B. and Wang, L. (2015), An Iterative Coordinate Descent Algorithm for High-Dimensional Nonconvex Penalized Quantile Regression. \emph{Journal of Computational and Graphical Statistics}, \bold{24}, 676-694.
    
    Dong, C., Li, G. D. and Feng, X. D. (2019), Lack-of-fit tests for quantile regression models, \emph{Journal of Royal Statistical Society, Series B}, \bold{To appear}.
}
\keyword{qrLOFT}
\examples{
  library(qrLOFT)
}
