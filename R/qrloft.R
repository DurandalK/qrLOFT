qrloft <- function(Y, X, beta_in, tau = 0.5, intercept = TRUE, type = NULL, B = 500, 
                   Cret = FALSE, ...){
    
    ###
    # A lack of fit test for high dimensional quantile regression
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
    #       intercept, logic,
    #           intercept term indicator;
    #       type, char,
    #           test type indicator, can be "LowDim" or "HighDim".
    #       ...
    #           other parameter passed to QRlof, such as Bootstrap sample size, B,
    #           correction indicator, Cret.
    # Output:
    #       test.stat, numeric,
    #           proposed test statistic;
    #       p.value, numeric,
    #           p value for proposed test.
    #
    ###
    
    if(is.data.frame(X)){
        X <- data.matrix(X)
    }
    n <- nrow(X)
    p <- ncol(X)
    
    if(is.null(type)){
        if(n / p > 10)
            type <- "LowDim"
        else
            type <- "HighDim"
    }
    
    if(type == "LowDim"){
        output <- .QRlof(X, Y, tau, intercept, B, Cret)
        Mnhat <- output[1]
        pval <- output[2]
    }else{
        if(missing(beta_in)){
            est <- qrfit(Y, X, tau = tau, intercept = intercept, ...)
            beta_in <- est$beta_hat
        }
        # compute estimated eps
        if(intercept){
            eps_hat <- drop(Y - cbind(rep(1, n), X) %*% beta_in)
        }else{
            eps_hat <- drop(Y - X %*% beta_in)
        }
        I <- seq(1, n)
        indx_hat <- I[eps_hat < 0]
        
        # separate the X
        if(length(indx_hat) == 0 | length(indx_hat) == n){
            stop("All estimated residuals are of same sign, an intercept term is suggested before re-estimation.")
        }        
        if(length(indx_hat) == 1){
            Xs <- as.matrix(t(X[indx_hat, ]))
            Xsc <- X[-indx_hat, ]
        }else if(length(indx_hat) == (n - 1)){
            Xs <- X[indx_hat, ]
            Xsc <- as.matrix(t(X[-indx_hat, ]))
        }else{
            Xs <- X[indx_hat, ]
            Xsc <- X[-indx_hat, ]
        }
        
        # Precision Matrix Estimation
        Beta <- .PrecisionMat(Xs, Xsc)
        # M_mu test statistics
        Muhat <- .MMuhat(Xs, Xsc, Beta)
        # M_sigma test statistics
        Mshat <- .MSighat(Xs, Xsc)
        # Proposed test statistics
        Mnhat <- max(Muhat - 2 * log(p) + log(log(p)), 
                     Mshat - 4 * log(p) + log(log(p)))
        
        # p-value
        pval <- 1 - exp(-(pi ^ (-0.5) + (8 * pi) ^ (-0.5)) * exp(-Mnhat/2))
    }
    
    return(list(test.stat = Mnhat, p.value = pval))
}