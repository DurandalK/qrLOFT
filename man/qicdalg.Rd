\name{qicdalg}
\alias{qicdalg}
\title{
    QICD algorithm for solving quantile regression coefficients with different penalties.
}
\description{
    This function \code{qicdalg} can be used to obtain quantile regression coefficient estimates with different penalty term, including LASSO, SCAD, MCP and L1-Penalty. More details about the algorithm can be found in Peng and Wang (2015).
}
\usage{
qicdalg(Y, X, beta_in = NULL, tau = 0.5, lambda = NULL, intercept = TRUE, 
        penalty = "L1-Penalty", a = 3.7, threshod = 1e-08, maxint = 100, 
        post.est = TRUE)
}
\arguments{
  \item{Y}{\emph{n} x 1 numeric vector, regression response.}
  \item{X}{\emph{n} x \emph{p} numeric matrix, \emph{p} dimensional predictors.}
  \item{beta_in}{\emph{p} x 1 or (\emph{p}+1) x 1 numeric vector, quantile regression coefficients vector.}
  \item{tau}{numeric value, current quantile level, should lie in (0, 1) and the default value is 0.5.}
  \item{lambda}{numeric value, current penalty level, should be (>0).}
  \item{intercept}{logical value, intercept term indicator, the default value is \code{TRUE}.}
  \item{penalty}{character string, penalty model indicator, one of "\code{SCAD}", "\code{MCP}", "\code{LASSO}" and "\code{Penalized QR}", the default value is "\code{L1-Penalty}".}
  \item{a}{numeric value, a pre-set constant used in "\code{SCAD}" and "\code{MCP}".}
  \item{threshod}{numeric value, a threshold value used for MM-algorithm iteration controlling.}
  \item{maxint}{integer, maximum iteration number for greedy coordinate descent and MM-algorithm.}
  \item{post.est}{logical value, post re-estimation indicator, the default value is \code{TRUE}.}
}
\details{
    We use greedy coordinate descent method to implement the QICD algorithm and perform post re-estimation suggested in Belloni and Chernozhukov (2011) by default. The post re-estimation stands for re-estimating quantile regression mode using only selected predictors, this re-estimation is performed without penalty.
    The penalty coefficient should be provided in the first place, if not, \code{lambda} is set to 0.
    The length of \code{beta_in} should be equal to the number of columns of \code{x}, denoted as \code{p}, if \code{intercept = 0}, or \code{length(beta_in)} is equal to \code{p + 1} if \code{intercept = 1}. A vector of length \code{length(beta_in)} is returned. 
    \code{a} is vaild only in the "\code{SCAD}" and "\code{MCP}" case.
}
\value{
  \item{beta_out}{\emph{p} x 1 or (\emph{p}+1) x 1 numeric vector, penalized quantile regression coefficients estimates.}
}
\references{
    Peng, B. and Wang, L. (2015), An Iterative Coordinate Descent Algorithm for High-Dimensional Nonconvex Penalized Quantile Regression; \emph{Journal of Computational and Graphical Statistics}, \bold{24}, 676-694.
    
    Fan, J. and Lv, J. (2001), A Selective Overview of Variable Selection in High Dimensional Feature Space, \emph{Statistica Sinica}, \bold{20}, 101-148.
    
    Zhang, C. H. (2010), Nearly Unbiased Variable Selection Under Minimax Concave Penalty, \emph{The Annals of Statistics}, \bold{38}, 894-942.
    
    Belloni, A. and Chernozhukov V. (2011) L1-penalized Quantile Regression in High-dimensional Sparse Models; \emph{The Annals of Statistics}, \bold{39}, 82-130.
    
    Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery, B. P. (2007) \emph{Numerical Recipes: The Art of Scientific Computing, 3rd ed.}, Cambridge University Press, Cambridge.
}
\seealso{
    \code{\link{qrfit}}
}
\examples{
data("GrowthData")
qicdalg(GrowthData[, 1], as.matrix(GrowthData[, 40 : 63]), tau = 0.75, lambda = 0.05)
}
\keyword{QICD}