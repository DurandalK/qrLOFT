#include <Rcpp.h>
#include "nr3.h"
#include "ludcmp.h"
#include "svd.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector qrmm(NumericMatrix X, NumericVector Y, double tau){
    /*
     * Input: 
     *      X, n x p matrix,
     *          design matrix for quantile regression, intercept term is required to be included as the first column if needed;
     *      Y, n x 1 vector,
     *          response vector;
     *      tau, numeric,
     *          quantile level, should be in (0, 1).
     * 
     *  Output: 
     *      beta_out, p x 1 vector,
     *          estimated coefficients for quantile regression, first element is the estimated intercept term.
     *          
     *  Refer: 
     *      Quantile Regression via an MM Algorithm (1999)
     */
    int iter_in, iter_out = 0, maxint = 500, n = X.nrow(), p = X.ncol(), i, j, k;
    double threshod = 1e-6, eps = 0.5 - sqrt(0.25 - threshod / n), surr, lastsurr = 0.0, upd = 1.0;
    
    VecDoub res(n), wgt(n), v(n), d_elta(p);
    MatDoub xWx(p, p), tmp(p, p);
    NumericVector beta_out(p);
    
    // obtain the initial estimates for beta_out (Least-square-estimator);
    for(j = 0; j < p; j++){
        for(k = j; k < p; k++){
            xWx[j][k] = 0.0;
            for(i = 0; i < n; i++)
                xWx[j][k] += X(i, j) * X(i, k);
            if(k > j)
                xWx[k][j] = xWx[j][k];
        }
    }
    LUdcmp M = LUdcmp(xWx);
    if(fabs(M.det()) > 1e-6){
        M.inverse(tmp);
    }else{
        for(k = 0; k < p; k++)
            M.lu[k][k] += 1e-4;
        M.inverse(tmp);
    }
    
    SVD S = xWx;
    if(S.rank() == p){
        LUdcmp M = LUdcmp(xWx);
        M.inverse(tmp);
    }else{
        Rprintf("A singular design matrix encounted QR regression. \n");
        return 0;
    }
    
    for(k = 0; k < p; k++){
        beta_out[k] = 0.0;
        for(i = 0; i < n; i++){
            for(j = 0; j < p; j++)
                beta_out[k] += tmp[k][j] * X(i, j) * Y[i];
        }
    }
    // calculate lastsurr (Q_eps(theta^k | theta^k));
    lastsurr = 0.0;
    for(i = 0; i < n; i++){
        res[i] = 0.0;
        for(j = 0; j < p; j++)
            res[i] += X(i, j) * beta_out[j];
        res[i] = Y[i] - res[i];
        wgt[i] = 1 / (eps + fabs(res[i]));
        v[i] = 1 - 2 * tau - res[i] * wgt[i];
        lastsurr += pow(res[i], 2) * wgt[i] + (4 * tau - 2) * res[i];
    }
    
    // main iteration
    while(upd > threshod && iter_out < maxint){
        for(j = 0; j < p; j++){
            for(k = j; k < p; k++){
                xWx[j][k] = 0.0;
                for(i = 0; i < n; i++)
                    xWx[j][k] += X(i, j) * X(i, k) * wgt[i];
                if(k > j)
                    xWx[k][j] = xWx[j][k];
            }
        }
        
        SVD S = xWx;
        if(S.rank() == p){
            LUdcmp M = LUdcmp(xWx);
            M.inverse(tmp);
        }else{
            Rprintf("A singular design matrix encounted QR regression. \n");
            return 0;
        }
        
        for(k = 0; k < p; k++){
            d_elta[k] = 0.0;
            for(i = 0; i < n; i++){
                for(j = 0; j < p; j++){
                    d_elta[k] += tmp[k][j] * X(i, j) * v[i];
                }
            }
            d_elta[k] = -d_elta[k];
        }
        
        // obtain the appropriate fractional step
        iter_in = 0;
        do{
            surr = 0.0;
            for(i = 0; i < n; i++){
                res[i] = 0.0;
                for(j = 0; j < p; j++)
                    res[i] += X(i, j) * (beta_out[j] + d_elta[j] * pow(0.5, iter_in));
                res[i] = Y[i] - res[i];
                surr += pow(res[i], 2) * wgt[i] + (4 * tau - 2) * res[i];
            }
            iter_in++;
        }while(surr > lastsurr && iter_in < maxint);

        // update the convegence condition and lastsurr;
        upd = lastsurr - surr;
        lastsurr = 0.0;
        for(k = 0; k < p; k++)
            beta_out[k] += d_elta[k] * pow(0.5, iter_in - 1);
        for(i = 0; i < n; i++){
            res[i] = 0.0;
            for(j = 0; j < p; j++)
                res[i] += X(i, j) * beta_out[j];
            res[i] = Y[i] - res[i];
            wgt[i] = 1 / (eps + fabs(res[i]));
            v[i] = 1 - 2 * tau - res[i] * wgt[i];
            lastsurr += pow(res[i], 2) * wgt[i] + (4 * tau - 2) * res[i];
        }
        
        iter_out++;
    }

    return beta_out;
}
