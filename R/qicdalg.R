qicdalg <- function(Y, X, beta_in = NULL, tau = 0.5, lambda = NULL, intercept = TRUE, 
                    penalty = "L1-Penalty", a = 3.7, threshod = 1e-08, maxint = 100, 
                    post.est = TRUE){
    
    # Input: 
    #       Y, (n x 1) vector,
    #           response;
    #       X, (n x p) matrix,
    #           predictors;
    #       beta_in, (p x 1 or (p + 1) x 1) vector,
    #           pre-estimated coefficients vector;
    #       tau, numeric,
    #           current quantile level, should lie in (0, 1);
    #       lambda, numeric,
    #           penalty coefficient;
    #       intercept, logical,
    #           intercept term indicator;
    #       penalty, char,
    #           penalty type indicator;
    #       a, numeric,
    #           a parameter used in SCAD;
    #       threshod, numeric,
    #           threshold value for QICD iteration convergence;
    #       maxint, integer,
    #           maximum iteration number for QICD algorithm.
    #       post.est, logical,
    #           indicator for whether do post-estimation.
    #
    # Output:
    #       beta_hat, (p x 1 or (p + 1) x 1) vector,
    #           QICD quantile regression coefficients estimates.
    #
    # Refer:
    #       An Iterative Coordinate Descent Algorithm for High-dimensional Nonconvex 
    #           Penalized Quantile Regression (2014)
    #
    
    p <- ncol(X)
    n <- nrow(X)
    
    if(is.null(beta_in)){
        if(intercept)
            beta_in <- rep(0, p + 1)
        else
            beta_in <- rep(0, p)
    }
    
    if(penalty == "SCAD"){
        model <- 0
    }else if(penalty == "MCP"){
        model <- 1
    }else if(penalty == "LASSO"){
        model <- 2
    }else if(penalty == "L1-Penalty"){
        model <- 3
    }else{
        stop("Wrong Penalty Function Specificed")
    }
    
    if(is.null(lambda))
        lambda <- 0
    
    beta_hat <- .Cqicd(Y, X, beta_in, tau, intercept, lambda, model, a, threshod, maxint)
    
    # re-estimate the non-zero coefficients without penalty term by author's suggestion
    if(post.est){
        if(intercept){
            idx1 <- which(beta_hat[-1] != 0)
            idx2 <- c(1, idx1 + 1)
            if(length(idx1) > 0){
                beta_tmp <- qrmm(cbind(rep(1, n), as.matrix(X[, idx1])), Y, tau)
                beta_hat[idx2] <- beta_tmp
            }
        }else{
            idx1 <- which(beta_hat != 0)
            if(length(idx1) > 0){
                beta_tmp <- qrmm(as.matrix(X[, idx1]), Y, tau)
                beta_hat[idx1] <- beta_tmp
            }
        }
    }
    
    return(beta_hat)
}