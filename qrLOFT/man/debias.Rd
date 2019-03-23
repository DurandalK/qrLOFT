\name{debias}
\alias{debias}
\title{
    Debias procedure for high dimensional quantile regression model.
}
\description{
    Function \code{debias} applies idea in Bradic and Kolar (2017) and performs partially debiase procedure to high dimensional quantile regression coefficients.
}
\usage{
debias(Y, X, beta_hat, tau = 0.5, intercept = TRUE, c_h = 0.01, 
       delta = 2.0, mat_thres = 1e-8, ...)
}
\arguments{
    \item{Y}{\emph{n} x 1 numeric vector, regression response.}
    \item{X}{\emph{n} x \emph{p} numeric matrix, \emph{p} dimensional predictors.}
    \item{beta_hat}{\emph{p} x 1 or (\emph{p}+1) x 1 numeric vector, pre-estimated quantile regression coefficients vector.}
    \item{tau}{numeric value, current quantile level, should lie in (0, 1) and the default value is 0.5.}
    \item{intercept}{logical value, intercept term indicator, the default value is \code{TRUE}.}
    \item{c_h}{numeric value, a constant used in bandwidth selection.}
    \item{delta}{numeric value, a constant used in precision matrix estimation, c.f. Cai, Liu and Xia (2014), the default value is 2.0.}
    \item{mat_thres}{numeric value, a threshold used for determining the number of non-zeron elements in \code{beta_hat}.}
    \item{...}{Further arguments to be passed to function \code{qrfit}, more info can be found in details.}
}
\details{
    To carry out the partial debias procedure, this function split the whole process into two steps, first is to estimate precision matrix c.f. Cai, Liu and Xia (2014) and second is to estiamte the \eqn{\\vartau}, c.f. Belloni, Chernozhukov and Kat (2013). 
    \code{c_h} is a constant used to control the bandwidth size. To catch the exception "Bandwidith h is too wide" a smaller \code{c_h} is preferred.
    More details about paramters in \code{qrfit} can be found in its document.
}
\value{
    \item{beta_amend}{\emph{p} x 1 or (\emph{p}+1) x 1 numeric vector, debaised high dimensional quantile regression coefficients.}
}
\references{
    Bradic, J. and Kolar, M. (2017), Uniform Inference for High-dimensional Quantile Regression: Linear Functionals and Regression Rank Scores. \emph{arXiv:1702.06209}.
    
    Belloni, A., Chernozhukov, V. and Kato, K. (2013), Robust inference in high-dimensional approximately sparse quantile regression models. arXiv preprint arXiv:1312.7186,
    
    Cai, T., Liu, W. and Xia, Y. (2014), Two-sample test of high dimensional means under dependence. \emph{Journal of Royal Statistical Society, Series B}, \bold{76}, 34-372.
}
\seealso{
    \code{\link{qrfit}}
}
\examples{
data("GrowthData")
out <- qrfit(GrowthData[, 1], as.matrix(GrowthData[, 25 : 40]), tau = 0.5, lambda = seq(0.01, 0.03, length.out = 100), criteria = 'BIC')
out_new <- debias(GrowthData[, 1], as.matrix(GrowthData[, 25 : 40]), out$beta_hat, tau = 0.5)
}
\keyword{debias}