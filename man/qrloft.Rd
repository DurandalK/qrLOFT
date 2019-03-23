\name{qrloft}
\alias{qrloft}
\title{
    Lack of Fit Test for High-dimensional Quantile Regression
}
\description{
    This function performs lack of fit test for either multivariate and high dimensional quantile regression.
}
\usage{
qrloft(Y, X, beta_in, tau = 0.5, intercept = TRUE, type = NULL, 
       B = 500, Cret = FALSE, ...)
}
\arguments{
    \item{Y}{\emph{n} x 1 numeric vector, regression response.}
    \item{X}{\emph{n} x \emph{p} numeric matrix, \emph{p} dimensional predictors.}
    \item{beta_in}{\emph{p} x 1 or (\emph{p}+1) x 1 numeric vector, quantile regression coefficients vector.}
    \item{tau}{numeric value, current quantile level, should lie in (0, 1) and the default value is 0.5.}
    \item{intercept}{logical value, intercept term indicator.}
    \item{type}{character string, type of test, can be "Multi" or "HighDim". if doesn't specify, then this function will choose the test type using (\emph{n}, \emph{p}) accordingly.}
    \item{B}{integer, bootstrap test sample size.}
    \item{Cret}{logical value, indicator of finite sample correction used in bootstrap.}
    \item{...}{Further arguments to be passed to function \code{.QRlof}, more info can be found in details.}
}
\details{
    The \code{beta_in} can be estimated by other functions such as \code{\link{qrfit}}, \code{\link{qrmm}} and functions in other packages, such as \code{quantreg}, \code{cobs}.
    Unless \code{type} is specified, the function will choose the test type based on the ratio of sample size \emph{n} and number of predctiors \emph{p}. If (\emph{n} / \emph{p} > 10) the function will perform a lack of fit test based on multivariate quantile regression case, same as specifing \code{type} as "\code{Multi}", and perform a lack of fit test based on high dimensional quantile regression case if (\emph{n} / \emph{p} <= 10), same as specifing \code{type} as "\code{HighDim}".
    In this multivariate test case, \code{qrloft} uses wild bootstrap method (Feng, He and Hu (2011)) to obtain the test p-value. Hence if \code{type = "Multi"}, one can specify the bootstrap sample size \code{B} and whether to perform a finte number correction setting parameter \code{Cret} to \code{FALSE} or \code{TRUE}. If neither is specified, \code{B = 500} and \code{Cret = FALSE} by default.
    In high-dimensional case, \code{beta_in} should be specified in the first place, if not, the test procedure will automatically estimate the high dimensional quantile regression and further parameter can be passed to \code{qrfit} through \code{...}.
}
\value{
    \item{test.stat}{numeric value, proposed test statistic.}
    \item{p.value}{numeric value, p value for proposed test.}
}
\references{
    Feng, X. D., He, X. M. and Hu, J. H. (2011), Wild bootstrap for quantile regression, \emph{Biometrika}, \bold{94}, 995-999.
    
    Dong, C., Li, G. D. and Feng, X. D. (2019), Lack-of-fit tests for quantile regression models, \emph{Journal of Royal Statistical Society, Series B}, \bold{to appear}.
}
\seealso{
    See \code{\link{qrfit}}, \code{\link{qicdalg}} for estimating quantile regression model.
}
\examples{
data("GrowthData")
y <- GrowthData[, 1]
x <- as.matrix(GrowthData[, 25 : 40])
out <- qrfit(y, x, tau = 0.65, lambda = seq(0.01, 0.03, length.out = 100), criteria = 'BIC')
t_out <- qrloft(y, x, out$beta_hat, tau = 0.65, type = "HighDim")
}
\keyword{Lack of Fit Test}
