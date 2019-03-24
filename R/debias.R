debias <- function(Y, X, beta_hat, tau = 0.5, intercept = TRUE, c_h = 0.01, delta = 2.0, 
                   mat_thres = 1e-8, ...){
    
    # Partially debias procedure for estimated penalized quantile regression coefficients
    #
    # Input:
    #       Y, (n x 1) vector,
    #           response;
    #       X, (n x p) matrix,
    #           predictors;
    #       beta_hat, (p x 1 or (p + 1) x 1) vector,
    #           pre-estimated coefficients vector;
    #       tau, numeric,
    #           current quantile level, should lie in (0, 1);
    #       intercept, logic,
    #           intercept term indicator;
    #       c_h, numeric,
    #           bandwith in estimating \vartau in Bradic and Kolar (2017);
    #       delta, numeric,
    #           parameter in precision matrix estimation, c.f. Cai, Liu and Xia (2014);
    #       mat_thres, numeric,
    #           mat_thres value used in precision matrix estimation,
    #       ..., 
    #           further arguments to be passed to qrfit().
    #
    # Output:
    #       beta_amend, (p x 1 or (p + 1) x 1) vector,
    #           debiased quantile regression coefficients estimates.
    #   
    # Refer:
    #       Uniform inference for high-dimensional quantile regression: linear funtionals
    #           and regression rank scores (2017).
    #
    
    n <- length(Y)
    p <- length(beta_hat)
    
    if(intercept){
        Xx <- cbind(rep(1, n), X)
    }else{
        Xx <- X
    }
    
    Xmean <- apply(Xx, 2, mean)
    Xred <- Xx - tcrossprod(rep(1, n), Xmean)
    Sigma_hat <- t(Xred) %*% Xred / (n - 1)
    Theta_hat <- matrix(0, p, p)
    
    for(ii in 1 : n){
        Theta_hat <- Theta_hat + (tcrossprod(Xred[ii, ]) - Sigma_hat) ^ 2
    }
    
    Lambda_hat <- delta * sqrt(Theta_hat * log(p) / n / n)
    Sigma_star <- Sigma_hat
    Sigma_star[abs(Sigma_hat) < Lambda_hat] <- 0
    eigs <- eigen(Sigma_star, symmetric = TRUE)
    eigval <- pmax(eigs$values, log(p) / n)
    Sigma_dag <- eigs$vectors %*% diag(eigval) %*% t(eigs$vectors)
    D_hat <- solve(Sigma_dag)
    
    s <- length(seq_len(p)[beta_hat > mat_thres])
    h <- c_h * (s * sqrt(log(max(p, n)) / n)) ^ (1/4)
    
    # Compute \tau_n
    if((h + tau) > 1 | (tau - h) < 0){
        stop("Bandwidith h is too wide")
    }
    up_h <- qrfit(Y, X, tau = tau + h, intercept = intercept, ...)
    down_h <- qrfit(Y, X, tau = tau - h, intercept = intercept, ...)
    taun <- drop(matrix(1, 1, n) %*% (Xx %*% (up_h$beta_hat - down_h$beta_hat)) / (2 * h * n))
    
    phi_fit <- drop(tau - (Y < Xx %*% beta_hat))
    appterm <- t(matrix(phi_fit, nrow = 1, ncol = n) %*% Xx)  / n
    
    idx <- seq_len(p)[beta_hat != 0]
    amend_est <- rep(0, p)
    amend_est[idx] <- drop(taun * (D_hat %*% appterm))[idx] / n
    beta_amend <- beta_hat + amend_est
    return(beta_amend)
}