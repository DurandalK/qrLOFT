\name{qrfit}
\alias{qrfit}
\title{
    High dimensional quantile regression estimation using QICD algorithm
}
\description{
    \code{qrfit} uses QICD in Peng and Wang (2015) to estimate high dimensional quantile regression model and apply different criterions to select the penalty coefficient.
}
\usage{
qrfit(Y, X, beta_in = NULL, tau = 0.5, lambda = seq(0.005, 0.3, length.out = 200),
      intercept = TRUE, criteria = 'Lambda', bic_c = 1, 
      nfold = 5, l1_c = 2, ...)
}
\arguments{
    \item{Y}{\emph{n} x 1 numeric vector, regression response.}
    \item{X}{\emph{n} x \emph{p} numeric matrix, \emph{p} dimensional predictors.}
    \item{beta_in}{\emph{p} x 1 or (\emph{p}+1) x 1 numeric vector, pre-estimated quantile regression coefficients vector.}
    \item{tau}{numeric value, current quantile level, should lie in (0, 1) and the default value is 0.5.}
    \item{lambda}{numeric vector, possible penalty level, should be (>0).}
    \item{intercept}{logical value, intercept term indicator, the default value is \code{TRUE}.}
    \item{criteria}{character string, one can choose one of "Lambda", "BIC", "CV" to determine the penalty coefficient selection procedure, the default value is "Lambda".}
    \item{bic_c}{numeric value, a constant used in BIC criterion computation.}
    \item{nfold}{integer, number of folds used in CV criterion.}
    \item{l1_c}{numeric value, a constant used in l1 penalty selection.}
    \item{...}{Further arguments to be passed to function \code{qicdalg}, more info can be found in details.}
}
\details{
    In penalty selction method, one can sepecify \code{criteria} to one of "Lambda", "BIC" and "CV". "Lambda" stands for using L1-penalty criteria in Chernozhukov et al. (2010); "BIC" stands for using BIC criteria to choose lambda; "CV" stands for using Cross validation to choose lambda. One should also provide a sequence of possible lambdas before estimation, the funtion will use the default sequence otherwise.
    bic_c is a controlling constant used in BIC evaluation, more details can be found in Peng and Wang (2015). l1_c is a controlling constant used in l1 penalty selection procedure, more details can be found in Belloni and Chernozhukov (2011).
    Further parameter used in \code{qicdalg} application can be found in \code{qicdalg} document.
}
\value{
    \item{beta_hat}{\emph{p} x 1 or (\emph{p}+1) x 1 numeric vector, penalized quantile regression coefficients estimates.}
    \item{lambda_min}{numeric value, corresponding penalty coefficient used in estimation.}
}
\references{
    Peng, B. and Wang, L. (2015), An Iterative Coordinate Descent Algorithm for High-Dimensional Nonconvex Penalized Quantile Regression. \emph{Journal of Computational and Graphical Statistics}, \bold{24}, 676-694.
    
    Belloni, A. and Chernozhukov V. (2011) L1-penalized Quantile Regression in High-dimensional Sparse Models. \emph{The Annals of Statistics}, \bold{39}, 82-130.
    
    Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery, B. P. (2007) \emph{Numerical Recipes: The Art of Scientific Computing, 3rd ed.}, Cambridge University Press, Cambridge.
}
\seealso{
    \code{\link{qicdalg}}
}
\examples{
data("GrowthData")
qrfit(GrowthData[, 1], as.matrix(GrowthData[, 25 : 40]), tau = 0.5, lambda = seq(0.01, 0.03, length.out = 100), criteria = 'BIC')
}
\keyword{QICD}
\keyword{High dimensional quantile regression}
