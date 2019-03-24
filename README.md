# qrLOFT
This package provides a fast way to implement the lack of fit tests for both multivariate and high dimensional quantile regression models. More details about these test methods can be found in Dong, Li and Feng (2019).

# Installation
```
#install.packages("devtools")
library(devtools)
install_github("DurandalK/qrLOFT")
```  
# Example
```
library(qrLOFT)
data("GrowthData")

y <- GrowthData[, 1]
x <- as.matrix(GrowthData[, 25 : 40])
out <- qrfit(y, x, tau = 0.75, lambda = seq(0.01, 0.03, length.out = 100), criteria = 'BIC')
beta_hat <- debias(y, x, out$beta_hat, tau = 0.75, c_h = 5e-3)
t_out <- qrloft(y, x, beta_hat, tau = 0.75, type = "HighDim")
```
# References
Peng, B. and Wang, L. (2015), An Iterative Coordinate Descent Algorithm for High-Dimensional Nonconvex Penalized Quantile Regression, Journal of Computational and Graphical Statistics, 24, 676-694.
    
Belloni, A. and Chernozhukov V. (2011) L1-penalized Quantile Regression in High-dimensional Sparse Models, The Annals of Statistics, 39, 82-130.

Hunter, D. R. and Lange, K. (2000), Quantile Regression via an MM Algorithm, Journal of Computational and Graphical Statistics, 9, 60-77.

Bradic, J. and Kolar, M. (2017), Uniform Inference for High-dimensional Quantile Regression: Linear Functionals and Regression Rank Scores, arXiv:1702.06209.

Cai, T., Liu, W. and Xia, Y. (2013), Two-Sample Covariance Matrix Testing and Support Recovery in High-Dimensional and Sparse Settings, Journal of the American Statistical Association, 501, 265-277.
    
Cai, T., Liu, W. and Xia, Y. (2014), Two-sample test of high dimensional means under dependence, Journal of Royal Statistical Society, Series B, 76, 34-372.

Feng, X. D., He, X. M. and Hu, J. H. (2011), Wild bootstrap for quantile regression, Biometrika, 94, 995-999.
    
Dong, C., Li, G. D. and Feng, X. D. (2019), Lack-of-fit tests for quantile regression models, Journal of Royal Statistical Society, Series B, to appear.

Barro, R. J. and Lee, J. W. (1994), Data set for a panel of 139 countries. NBER.

Barro, R. J. and Sala-i-Martin, X. (1995), Economic Growth. McGrwa-Hill, New York.

Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery, B. P. (2007) Numerical Recipes: The Art of Scientific Computing, 3rd ed, Cambridge University Press, Cambridge.

# Development
This R package is maintained by Chen Dong (dongchen39@hotmail.com).