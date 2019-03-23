#include "QRMM.h"
#include "sort.h"
#include "eigen_sym.h"
#include "ProTest.h"
#include "akj.h"
using namespace Rcpp;

double PenltPr(int model, double beta_in, double lambda, double a){
    /*
    * Input: 
    *      model, integer,
    *          penalty iterm idicator, 0 for "SCAD", 1 for "MCP, 2 for "LASSO", others for "Penalized QR";
    *      beta_in, numeric,
    *          the i-th running coefficient;
    *      lambda, numeric,
    *          the current lambda in use;
    *      a, numeric,
    *          a pre-set constant used in "SCAD" and "MCP or a calculated weight for "Penalized QR".
    * 
    *  Output: 
    *      pprime, numericpenalty functions' primal form.
    */
    double pprime;
    if(model == 0){     // SCAD
        if(fabs(beta_in) >= 0 && fabs(beta_in) < lambda)
            pprime = (beta_in >= 0) ? lambda : -lambda;
        else if(fabs(beta_in) >= lambda && fabs(beta_in) <= a * lambda)
            pprime = (beta_in >= 0) ? (a * lambda - beta_in) / (a - 1) : (-a * lambda - beta_in) / (a - 1);
        else
            pprime = 0;
    }else if(model == 1){      // MCP
        if(fabs(beta_in) >= 0 && fabs(beta_in) < a * lambda)
            pprime = (beta_in >= 0) ? (lambda - beta_in / a) : (-lambda - beta_in / a);
        else
            pprime = 0;
    }else if(model == 2){      // LASSO
        pprime = (beta_in >= 0) ? lambda : -lambda;
    }else{      // Penalized QR
        pprime = (beta_in >= 0) ? (lambda * a) : (-lambda * a);
    }
    return pprime;
}

// [[Rcpp::export(.Cqicd)]]
NumericVector Cqicd(NumericVector Y, NumericMatrix X, NumericVector beta_in, double tau, int intercept, double lambda, int model, double a, double threshod, int maxint){
    /*
    * Non-convex penalized quantile regression method with SCAD, MCP or LASSO (greedy coordinate descent algorithm)
    * 
    * Input: 
    *       Y, (n x 1) vector,
    *           response;
    *       X, (n x p) matrix,
    *           predictor design matrix;
    *       beta_in, (p x 1 or (p + 1) x 1) vector,
    *           pre-estimated coefficients vector;
    *       tau, numeric,
    *           current quantile level, should lie in (0, 1);
    *       intercept, integer,
    *           intercept term indicator;
    *       lambda, numeric,
    *           current penalty level, should be (>0);
    *       model, integer,
    *           penalty model indicator, 0 for "SCAD", 1 for "MCP, 2 for "LASSO", 3 for "Penalized QR";
    *       a, numeric,
    *           a pre-set constant used in "SCAD" and "MCP;
    *       threshod, numeric,
    *           a threshold value used for MM-algorithm iteration controlling;
    *       maxint, integer,
    *           maximum iteration number for greedy coordinate descent and MM-algorithm.
    *           
    * Output:
    *       beta_out, (p x 1 or (p + 1) x 1) vector,
    *           estimated coefficient vector.
    * 
    * Refer: 
    *      An Iterative Coordinate Descent Algorithm for High-dimensional Nonconvex Penalized Quantile Regression (2014)
    *      Coordinate Descent Algorithm For LASSO Penalized Regression (2008)
    */
    double diffout, diffin, w_sum, tmp, sigma;
    int p = X.ncol(), n = X.nrow(), l = beta_in.size(), i, j, r, k = 0, idx;
    VecDoub z(n + 1), weight(n + 1), resd(n), beta_u(p + 1), beta_r(p + 1), beta_k(p + 1), aow(p + 1);
    MatDoub Xx(n, p + 1);
    
    NumericVector beta_out(l);
    
    double PenltPr(int model, double beta_in, double lambda, double a);
    
    if(model == 3){ // Calculate penalty weights sigmahat for penalized QR
        aow[0] = 0.0;
        for(j = 0; j < p; j++){
            sigma = 0.0;
            for(i = 0; i < n; i++){
                sigma += pow(X(i, j), 2);
            }
            aow[j + 1] = sqrt(sigma / n);
        }
    }else{
        for(j = 0; j < p + 1; j++){
            aow[j] = a;
        }
    }
    
    // Initialize output and in-run term Xx ([0, X] or [1, X]) Yy and beta_r ([0; beta_in] or [1; beta_in]).
    for(i = 0; i < n; i++){
        Xx[i][0] = intercept;
        for(j = 0; j < p; j++){
            Xx[i][j + 1] = X(i, j);
        }
    }
    beta_u[0] = 0.0;
    for(j = 0; j < l; j++){
        beta_u[j + 1 - intercept] = beta_in[j];
    }
    
    do{
        for(j = 0; j < p + 1; j++){
            beta_k[j] = beta_u[j];
        }
        diffout = 0.0;
        r = 0;
        do{
            for(j = 0; j < p + 1; j++){
                beta_r[j] = beta_u[j];
            }
            diffin = 0.0;
            for(idx = 0; idx < p + 1; idx++){
                for(i = 0; i < n; i++){
                    resd[i] = Y[i];
                    for(j = 0; j < p + 1; j++){
                        resd[i] -= Xx[i][j] * beta_u[j];
                    }
                }
                w_sum = 0.0;        // Weighted median method for updating beta_hat
                for(i = 0; i < n; i++){
                    weight[i] = (resd[i] > 0.0) ? (fabs(Xx[i][idx]) * tau / n) : (fabs(Xx[i][idx]) * (1 - tau) / n);
                    z[i] = (idx > 0) ? (resd[i] / Xx[i][idx] + beta_u[idx]) : (resd[i] + beta_u[0]);
                    w_sum += weight[i];
                }
                z[n] = 0.0;
                weight[n] = (idx > 0) ? PenltPr(model, fabs(beta_k[idx]), lambda, aow[idx]) : 0.0;
                w_sum += weight[n];
                sort2(z, weight);
                i = 0;
                tmp = w_sum - weight[0];
                while(tmp > (w_sum / 2)){
                    ++i;
                    tmp -= weight[i];
                }
                beta_u[idx] = z[i];
            }
            for(j = 0; j < p + 1; j++)
                diffin += pow((beta_u[j] - beta_r[j]), 2); 
            diffin = sqrt(diffin);
            r++;
        }while(r < maxint && diffin > threshod);
        for(j = 0; j < p + 1; j++)
            diffout += pow((beta_u[j] - beta_k[j]), 2); 
        diffout = sqrt(diffout);
        k++;
    }while(k < maxint && diffout > threshod);
    
    for(i = 0, j = 1 - intercept; j < p + 1; i++, j++){
        beta_out[i] = beta_u[j];
    }
    
    return beta_out;
}

NumericVector checkloss(NumericVector res, double tau){
    /*
    * Check-loss function
    * 
    * Input: 
    *      res, (n x 1) vector,
    *          residuals;
    *      tau, numeric,
    *          quantile level, should be in (0, 1).
    *  
    *  Output:
    *      out, (n x 1) vector,
    *          check-loss function value.
    */
    double pos, neg;
    int l = res.size(), i;
    NumericVector out(l);
    
    for(i = 0; i < l; i++){
        pos = (fabs(res[i]) + res[i]) / 2;
        neg = (fabs(res[i]) - res[i]) / 2;
        out[i] = (tau * pos + (1 - tau) * neg);
    }
    return out;
}

// [[Rcpp::export(.BicVal)]]
double BicVal(NumericVector Y, NumericMatrix X, NumericVector beta_in, double tau, int intercept, double C){
    /*
    * BIC value evaluation
    * 
    * Input: 
    *      Y, (n x 1) vector,
    *          response;
    *      X, (n x p) matrix,
    *          predictor design matrix;
    *      beta_in, (p x 1 or (p + 1) x 1) vector,
    *          pre-estimated coefficients vector;
    *      tau, numeric,
    *          current quantile level, should lie in (0, 1);
    *      intercept, integer,
    *          intercept term indicator;
    *      C, numeric,
    *          a controlling constant used in BIC evaluation, more details can be found in the literature.
    *          
    *  Output:
    *      BIC, numeric
    *          BIC value.
    * 
    * Refer: 
    *      An Iterative Coordinate Descent Algorithm for High-dimensional Nonconvex Penalized Quantile Regression(2014)
    */
    int n = X.nrow(), p = X.ncol(), l = beta_in.size(), df = 0, i, j;
    double BIC = 0.0;
    NumericVector res(n), fitval(n);
    
    NumericVector checkloss(NumericVector res, double tau);
    
    for(j = 0; j < l; j++){
        df += (fabs(beta_in[j]) > 1e-6) ? 1 : 0;
    }
    
    for(i = 0; i < n; i++){
        res[i] = intercept * beta_in[0];
        for(j = 0; j < p; j++){
            res[i] += X(i, j) * beta_in[j + intercept];                
        }
        res[i] = Y[i] - res[i]; 
    }
    fitval = checkloss(res, tau);
    
    for(i = 0; i < n; i++){
        BIC += fitval[i];
    }
    
    BIC = log(BIC) + df * log(log(n)) / n * (C * log(p));
    
    return BIC;
}

// [[Rcpp::export(.PrecisionMat)]]
NumericMatrix PrecisionMat(NumericMatrix Xs, NumericMatrix Xsc){
    /* 
    * Precision matrix estimation
    * 
    * Input: 
    *       Xs, n1 x p matrix,
    *           predictors design matrix whose fitted residuals are above(below) zero;
    *       Xsc, n2 x p matrix,
    *           predictors design matrix whose fitted residuals are below(above) zero;
    * Output: 
    *       Beta, p x p, matrix,
    *           estimated precision matrix;
    * 
    * Refer: 
    *       "Two-sample test of high dimensional means under dependence" (2014)
    *       "Adaptive thresholding for sparse covariance matrix estimation." (2011)
    */
    int p = Xs.ncol(), n1 = Xs.nrow(), n2 = Xsc.nrow(), i1, i2, j1, j2, k;
    double C = 2.0 * sqrt(log(p) / (n1 + n2)), threshold;
    MatDoub S_igma(p, p), T_heta(p, p), tmp1(p, p), tmp2(p, p), evtor(p, p);
    VecDoub Xsbar(p), Xscbar(p), diags(p), eval(p);
    
    NumericMatrix O_mega(p, p);
    
    // Xs_bar, Xsc_bar;
    for(j1 = 0; j1 < p; j1++){
        Xsbar[j1] = 0.0;
        Xscbar[j1] = 0.0;
        for(i1 = 0; i1 < n1; i1++){
            Xsbar[j1] += Xs(i1, j1);
        }
        Xsbar[j1] /= n1;
        for(i2 = 0; i2 < n2; i2++){
            Xscbar[j1] += Xsc(i2, j1);
        }
        Xscbar[j1] /= n2;
    }
    
    // Sigma_hat matrix estimation (up-half);
    for(j1 = 0; j1 < p; j1++){
        for(j2 = j1; j2 < p; j2++){
            S_igma[j1][j2] = 0.0;
            for(i1 = 0; i1 < n1; i1++){
                S_igma[j1][j2] += (Xs(i1, j1) - Xsbar[j1]) * (Xs(i1, j2) - Xsbar[j2]);
            }
            for(i2 = 0; i2 < n2; i2++){
                S_igma[j1][j2] += (Xsc(i2, j1) - Xscbar[j1]) * (Xsc(i2, j2) - Xscbar[j2]);   
            }
            S_igma[j1][j2] /= (n1 + n2);
        }
    }
    
    // store the original diagonal values in Sigma_hat (un-shrink version)
    for(j1 = 0; j1 < p; j1++){
        diags[j1] = S_igma[j1][j1];
    }
    
    // Theta_hat matrix estimation (up-half), update Sigma_hat to Sigma_hat^star (whole);
    for(j1 = 0; j1 < p; j1++){
        for(j2 = j1; j2 < p; j2++){
            T_heta[j1][j2] = 0.0;
            for(i1 = 0; i1 < n1; i1++)
                T_heta[j1][j2] += pow((Xs(i1, j1) - Xsbar[j1]) * (Xs(i1, j2) - Xsbar[j2]) - S_igma[j1][j2], 2);
            for(i2 = 0; i2 < n2; i2++)
                T_heta[j1][j2] += pow((Xsc(i2, j1) - Xscbar[j1]) * (Xsc(i2, j2) - Xscbar[j2]) - S_igma[j1][j2], 2);
            T_heta[j1][j2] /= (n1 + n2);
            if(fabs(S_igma[j1][j2]) < sqrt(T_heta[j1][j2]) * C){
                S_igma[j1][j2] = 0.0;
            }
            if(j2 > j1){
                S_igma[j2][j1] = S_igma[j1][j2];
            }
        }
    }
    
    // Obatin Omega_hat by solve Sigma_hat^star;
    
    Symmeig J(S_igma);
    eval = J.d;
    evtor = J.z;
    threshold = log(p) * (n1 + n2) / (n1 * n2);
    
    for(i1 = 0; i1 < p; i1 ++){
        eval[i1] = (eval[i1] > threshold) ? eval[i1] : threshold;
    }
    
    for(j1 = 0; j1 < p; j1++){
        for(j2 = j1; j2 < p; j2++){
            tmp1[j1][j2] = 0.0;
            for(k = 0; k < p; k++){
                tmp1[j1][j2] += evtor[j1][k] * eval[k] * evtor[j2][k];
            }
            tmp1[j2][j1] = tmp1[j1][j2];
        }
    }
    
    LUdcmp M = LUdcmp(tmp1);
    M.inverse(tmp2);
    
    for(j1 = 0; j1 < p; j1++)
        for(j2 = 0; j2 < p; j2++)
            O_mega(j1, j2) = tmp2[j1][j2];
    
    return O_mega;
}

// [[Rcpp::export(.MMuhat)]]
double MMuhat(NumericMatrix Xs, NumericMatrix Xsc, NumericMatrix Beta){
    /* 
    * M_mu test statistics computation
    * 
    * Input: 
    *       Xs, n1 x p matrix,
    *           predictors design matrix whose fitted residuals are above(below) zero;
    *       Xsc, n2 x p matrix,
    *           predictors design matrix whose fitted residuals are below(above) zero;
    * Output: 
    *       Mmu_hat, numeric,
    *           numeric, test statistics for high dimensional means.
    * 
    * Source: "Two-sample test of high dimensional means under dependence" (2014)
    */
    int n1 = Xs.nrow(), n2 = Xsc.nrow(), p = Xs.ncol(), i1, i2, j1, j2;
    double Mmu_hat;
    MatDoub Xshat(n1, p), Xschat(n2, p), bijs(p, p), bijsc(p, p);
    VecDoub Xshatbar(p), Xschatbar(p), bii(p), Dhat(p), tmp(p);
    
    // Xshatbar, Xschatbar;
    for(j1 = 0; j1 < p; j1++){
        Xshatbar[j1] = 0.0;
        for(i1 = 0; i1 < n1; i1++){
            Xshat[i1][j1] = 0.0;
            for(j2 = 0; j2 < p; j2++){
                Xshat[i1][j1] += Xs(i1, j2) * Beta(j2, j1);
            }
            Xshatbar[j1] += Xshat[i1][j1];
        }
        Xshatbar[j1] /= n1;
        Xschatbar[j1] = 0.0;
        for(i2 = 0; i2 < n2; i2++){
            Xschat[i2][j1] = 0.0;
            for(j2 = 0; j2 < p; j2++){
                Xschat[i2][j1] += Xsc(i2, j2) * Beta(j2, j1);
            }
            Xschatbar[j1] += Xschat[i2][j1];
        }
        Xschatbar[j1] /= n2;
    }
    
    // bijs and bijsc;
    for(j1 = 0; j1 < p; j1++){
        for(j2 = j1; j2 < p; j2++){
            bijs[j1][j2] = 0.0;
            bijsc[j1][j2] = 0.0;
            for(i1 = 0; i1 < n1; i1++){
                bijs[j1][j2] += (Xshat[i1][j1] - Xshatbar[j1]) * (Xshat[i1][j2] - Xshatbar[j2]);
            }
            for(i2 = 0; i2 < n2; i2++){
                bijsc[j1][j2] += (Xschat[i2][j1] - Xschatbar[j1]) * (Xschat[i2][j2] - Xschatbar[j2]);
            }
            if(j2 > j1){
                bijs[j2][j1] = bijs[j1][j2];
                bijsc[j2][j1] = bijsc[j1][j2];
            }
        }
    }
    
    // bii, D_hat;
    for(j1 = 0; j1 < p; j1++){
        bii[j1] = bijs[j1][j1] + bijsc[j1][j1];
        Dhat[j1] = Xshatbar[j1] - Xschatbar[j1];
        tmp[j1] = pow(Dhat[j1], 2) / bii[j1];
    }
    Mmu_hat = tmp[0];
    for(j1 = 1; j1 < p; j1++){
        Mmu_hat = (Mmu_hat > tmp[j1]) ? Mmu_hat : tmp[j1];
    }
    Mmu_hat *= (n1 * n2);
    return Mmu_hat;	
}

// [[Rcpp::export(.MSighat)]]
double MSighat(NumericMatrix Xs, NumericMatrix Xsc){
    /* 
    * M_sigma test statistics computation
    * 
    * Input: 
    *       Xs, n1 x p matrix,
    *           predictors design matrix whose fitted residuals are above(below) zero;
    *       Xsc, n2 x p matrix,
    *           predictors design matrix whose fitted residuals are below(above) zero;
    * Output: 
    *       Msig_hat,
    *           numeric, test statistics for high dimensional covariance matrix.
    * 
    * Source: "Two-Sample Covariance Matrix Testing and Support Recovery in High-Dimensional and Sparse Settings" (2013)
    */
    double Msig_hat;
    int n1 = Xs.nrow(), n2 = Xsc.nrow(), p = Xs.ncol(), i1, i2, j1, j2;
    VecDoub Xsbar(p), Xscbar(p), arry(p), maxvec(p);
    MatDoub Sigs(p, p), Sigsc(p, p), Gams(p, p), Gamsc(p, p);
    
    // Xs_bar, Xsc_bar;
    for(j1 = 0; j1 < p; j1++){
        Xsbar[j1] = 0;
        Xscbar[j1] = 0;
        for(i1 = 0; i1 < n1; i1++){
            Xsbar[j1] += Xs(i1, j1);
        }
        Xsbar[j1] /= n1;
        for(i2 = 0; i2 < n2; i2++){
            Xscbar[j1] += Xsc(i2, j1);
        }
        Xscbar[j1] /= n2;
    }
    
    // Sigma_s_hat, Sigma_sc_hat matrix estimation;
    for(j1 = 0; j1 < p; j1++){
        for(j2 = j1; j2 < p; j2++){
            Sigs[j1][j2] = 0;
            Sigsc[j1][j2] = 0;
            for(i1 = 0; i1 < n1; i1++)
                Sigs[j1][j2] += (Xs(i1, j1) - Xsbar[j1]) * (Xs(i1, j2) - Xsbar[j2]);
            Sigs[j1][j2] /= n1;
            for(i2 = 0; i2 < n2; i2++)
                Sigsc[j1][j2] += (Xsc(i2, j1) - Xscbar[j1]) * (Xsc(i2, j2) - Xscbar[j2]);
            Sigsc[j1][j2] /= n2;
            if(j2 > j1){
                Sigs[j2][j1] = Sigs[j1][j2];
                Sigsc[j2][j1] = Sigsc[j1][j2];
            }
        }
    }
    
    // Gamma_s_hat, Gamma_sc_hat matrix estimation;
    for(j1 = 0; j1 < p; j1++){
        for(j2 = j1; j2 < p; j2++){
            Gams[j1][j2] = 0;
            Gamsc[j1][j2] = 0;
            for(i1 = 0; i1 < n1; i1++)
                Gams[j1][j2] += pow((Xs(i1, j1) - Xsbar[j1]) * (Xs(i1, j2) - Xsbar[j2]) - Sigs[j1][j2], 2);
            Gams[j1][j2] /= n1;
            for(i2 = 0; i2 < n2; i2++)
                Gamsc[j1][j2] += pow((Xsc(i2, j1) - Xscbar[j1]) * (Xsc(i2, j2) - Xscbar[j2]) - Sigsc[j1][j2], 2);
            Gamsc[j1][j2] /= n2;
            if(j2 > j1){
                Gams[j2][j1] = Gams[j1][j2];
                Gamsc[j2][j1] = Gamsc[j1][j2];
            }
        }
    }
    
    // fractions
    for(j1 = 0; j1 < p; j1++){
        arry[0] = pow(Sigs[j1][0] - Sigsc[j1][0], 2) / (Gams[j1][0] / n1 + Gamsc[j1][0] / n2);
        maxvec[j1] = arry[0];
        for(j2 = 1; j2 < p; j2++){
            arry[j2] = pow(Sigs[j1][j2] - Sigsc[j1][j2], 2) / (Gams[j1][j2] / n1 + Gamsc[j1][j2] / n2);
            maxvec[j1] = (maxvec[j1] > arry[j2]) ? maxvec[j1] : arry[j2];
        }
    }
    Msig_hat = maxvec[0];
    for(j1 = 1; j1 < p; j1++)
        Msig_hat = (Msig_hat > maxvec[j1]) ? Msig_hat : maxvec[j1];
    
    return Msig_hat;
}

// [[Rcpp::export(.QRlof)]]
NumericVector QRlof(NumericMatrix X, NumericVector Y, double tau, int intercept, int B, int Cret){
    /*
     * A lack of fit test for multivariate quantile regression
     * 
     * Input: 
     *      Y, (n x 1) vector,
     *          response;
     *      X, (n x p) matrix,
     *          predictor design matrix;
     *      tau, numeric,
     *          current quantile level, should lie in (0, 1);
     *      intercept, integer,
     *          intercept term indicator; 
     *      B, integer,
     *          Bootstrap resample size;
     *      Cret, integer,
     *          indicator for using residual crrection, c.f. Feng et al. (2011)
     *  Output: 
     *      test_stat, 2 x 1 vector.
     *          first element is the test statistic,
     *          second element is the bootstrap p-value.
     *      
     */
    int n = X.nrow(), p = X.ncol(), i, j, l = p + intercept, m;
    double stat_bstrp = 0.0, Tstat;
    
    NumericMatrix Xx(n, l);
    NumericVector test_stat(2), beta_est(l);
    
    for(i = 0; i < n; i++){
        Xx(i, 0) = intercept;
        for(j = 0; j < p; j++){
            Xx(i, j + intercept) = X(i, j);
        }
    }
    
    beta_est = qrmm(Xx, Y, tau);
    Tstat = Tlof(Xx, Y, beta_est, tau);
    test_stat[0] = Tstat;
    
    if(B > 1){
        if(Cret == TRUE){
            NumericVector beta_bstrp(l), fitval(n), res(n), rnd(n), wgt(n), Yfit(n), vec1(l), vec2(n), ress(n), px(n), xlam(n);
            MatDoub sXX(l, l), tmp(l, l);
            double z = 0.0, f0 = 1.0, psi = 1.0, score = 1.0, h = -1.0, alpha = 0.5, kappa = 0.9;
            int iker = 0, nx = n, nz = 1;
            
            RNGScope scope;
            for(i = 0; i < n; i++){
                fitval[i] = 0.0;
                for(j = 0; j < l; j++)
                    fitval[i] += Xx(i, j) * beta_est[j];
                res[i] = Y[i] - fitval[i];
                ress[i] = res[i];
            }
            
            for(i = 0; i < l; i++){
                for(j = i; j < l; j++){
                    sXX[i][j] = 0.0;
                    sXX[j][i] = 0.0;
                }
            }
            for(m = 0; m < n; m++){
                for(i = 0; i < l; i++){
                    for(j = i; j < l; j++){
                        sXX[i][j] += Xx(m, i) * Xx(m, j);
                        sXX[j][i] = sXX[i][j];
                    }
                }
                px[m] = 1.0 / n;
                xlam[m] = 1.0;
            }
            
            SVD S = sXX;
            if(S.rank() == l){
                LUdcmp M = LUdcmp(sXX);
                M.inverse(tmp);
                std::sort(ress.begin(), ress.end());
                akj_(ress.begin(), &z, px.begin(), &iker, &f0, &psi, &score, &nx, &nz, &h, &alpha, &kappa, xlam.begin());
                for(m = 0; m < n; m++){
                    vec2[m] = 0.0;
                    for(j = 0; j < l; j++){
                        vec1[j] = 0.0;
                        for(i = 0; i < l; i++)
                            vec1[j] += Xx(m, i) * tmp[i][j];
                        vec2[m] += vec1[j] * Xx(m, j);
                    }
                    res[m] += vec2[m] / f0 * (tau - ((res[m] < 0) ? 1 : 0));
                }
            }
            
            for(m = 0; m < B; m++){
                rnd = runif(n, 0, 1);
                for(i = 0; i < n; i++){
                    wgt[i] = (rnd[i] >= tau) ? (2.0 * (1 - tau)) : (-2.0 * tau);
                    Yfit[i] = fitval[i] + fabs(res[i]) * wgt[i];
                }
                beta_bstrp = qrmm(Xx, Yfit, tau);
                Tstat = Tlof(Xx, Yfit, beta_bstrp, tau);
                stat_bstrp += (Tstat >= test_stat[0]) ? 1.0 : 0.0;
            }
            test_stat[1] = stat_bstrp / B;
        }else{
            NumericVector beta_bstrp(l), fitval(n), res(n), rnd(n), wgt(n), Yfit(n);
            
            RNGScope scope;
            for(i = 0; i < n; i++){
                fitval[i] = 0.0;
                for(j = 0; j < l; j++)
                    fitval[i] += Xx(i, j) * beta_est[j];
                res[i] = Y[i] - fitval[i];
            }
            
            for(m = 0; m < B; m++){
                rnd = runif(n, 0, 1);
                for(i = 0; i < n; i++){
                    wgt[i] = (rnd[i] >= tau) ? (2.0 * (1 - tau)) : (-2.0 * tau);
                    Yfit[i] = fitval[i] + fabs(res[i]) * wgt[i];
                }
                beta_bstrp = qrmm(Xx, Yfit, tau);
                Tstat = Tlof(Xx, Yfit, beta_bstrp, tau);
                stat_bstrp += (Tstat >= test_stat[0]) ? 1.0 : 0.0;
            }
            test_stat[1] = stat_bstrp / B;
        }
    }else{
        test_stat[1] = 1;
    }
    
    return test_stat;
}
