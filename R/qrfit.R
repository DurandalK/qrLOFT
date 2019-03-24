qrfit <- function(Y, X, beta_in = NULL, tau = 0.5, lambda = seq(0.005, 0.3, length.out = 200), 
                  intercept = TRUE, criteria = 'Lambda', bic_c = 1, nfold = 5, 
                  l1_c = 2, ...){
    
    ###
    # QICD algorithm for solving quantile regression coefficients with different
    #       penalties.
    #
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
    #       intercept, logic,
    #           intercept term indicator;
    #       criteria, character,
    #           'Lambda' stands for using L1-penalty criteria in Chernozhukov et al. (2010) 
    #           'BIC' stands for using BIC criteria to choose lambda,
    #           'CV' stands for using Cross validation to choose lambda;
    #       nfold, integer,
    #           number of fold in CV criteria;
    #       bic_c, numeric,
    #           a constant in BIC criteria, c.f. Peng and Wang (2014);
    #       l1_c, numeric,
    #           a constant in L1-penalty criteria, c.f. Chernozhukov et al. (2010);
    #       ...,
    #           further arguments to be passed to qicdalg().
    #
    # Output:
    #       beta_hat, (p x 1 or (p + 1) x 1) vector,
    #           penalized quantile regression coefficients estimates;
    #       lambda_min, numeric,
    #           penalty coefficient used during estimation.
    #
    # Refer:
    #       An Iterative Coordinate Descent Algorithm for High-dimensional Nonconvex 
    #           Penalized Quantile Regression (2014).
    #
    ###
    
    if(is.data.frame(X)){
        X <- data.matrix(X)
    }
    n <- nrow(X)
    p <- ncol(X)
    
    l <- length(lambda)
    
    beta_collect <- NULL
    crit_val <- NULL
    
    if(criteria == "Lambda"){
        sigma_hat <- apply(X ^ 2, 2, mean)
        M <- 1000
        lams <- rep(0, M)
        for(i in 1 : M){
            u <- runif(n)
            psi <- tcrossprod(tau - (u < tau), 1 / (sqrt(sigma_hat * tau * (1 - tau))))
            lams[i] <- max(abs(apply(X * psi, 2, mean))) * n
        }
        lambda_min <- l1_c * quantile(lams, 0.9) * sqrt(tau * (1 - tau)) / n
        beta_hat <- qicdalg(Y, X, tau = tau, lambda = lambda_min, intercept = intercept, ...)
    }else if(criteria == "BIC"){
        for(k in 1 : l){
            beta_new <- qicdalg(Y, X, tau = tau, lambda = lambda[k], intercept = intercept, ...)
            beta_collect <- cbind(beta_collect, beta_new)
            BIC <- .BicVal(Y, X, beta_new, tau, intercept, bic_c)
            crit_val <- c(crit_val, BIC)
        }
    }else{
        n1 <- floor(n / nfold)
        n2 <- n - n1 * nfold
        for(k in 1 : l){
            CV <- 0
            for(i in 1 : ((n - n2) / n1)){
                idx <- ((i - 1) * n1 + 1) : (i * n1)
                beta_new <- qicdalg(Y[-idx], X[-idx, ], beta_in, tau, lambda[k], intercept, ...)
                if(intercept){
                    CV <- CV + sum(Y[idx] - cbind(rep(1, n1), X[idx, ]) %*% beta_new) / n1
                }else{
                    CV <- CV + sum(Y[idx] - X[idx, ] %*% beta_new) / n1
                }
            }
            if(n2 != 0){
                idx <- (n - n2 + 1) : n
                beta_new <- qicdalg(Y[-idx], X[-idx, ], beta_in, tau, lambda[k], intercept, ...)
                if(intercept){
                    CV <- CV + sum(Y[idx] - cbind(rep(1, n2), X[idx, ]) %*% beta_new) / n2
                }else{
                    CV <- CV + sum(Y[idx] - X[idx, ] %*% beta_new) / n2
                }
            }
            beta_new <- qicdalg(Y, X, beta_in, tau, lambda[k], intercept)
            beta_collect <- cbind(beta_collect, beta_new)
            crit_val <- c(crit_val, CV)
        }
    }
    
    # find minimum criterion values
    if(!is.null(crit_val)){
        idx <- which.min(crit_val)
        if(length(idx) == 1){
            lambda_min <- lambda[idx]
            beta_hat <- beta_collect[, idx] 
        }else{
            lambda_min <- min(lambda[idx])
            beta_hat <- beta_collect[, which(lambda == lambda_min)]
        }
    }
    
    return(list(beta_hat = beta_hat, lambda_min = lambda_min))
}