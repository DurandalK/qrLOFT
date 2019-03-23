using namespace Rcpp;

double ENorm(NumericVector X, NumericVector Y){
    int p = X.size();
    double norm = 0.0;
    
    for(int i = 0; i < p; i++){
        norm += pow(X[i] - Y[i], 2);
    }
    norm = sqrt(norm);
    return norm;
}

double Tlof(NumericMatrix X, NumericVector Y, NumericVector beta_in, double tau){
    
    int n = X.nrow(), p = X.ncol(), i, j, k;
    NumericVector fitval(n), res(n), tmpi(p), tmpj(p);
    double tstat, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
    
    for(i = 0; i < n; i++){
        fitval[i] = 0.0;
        for(j = 0; j < p; j++){
            fitval[i] += X(i, j) * beta_in[j];         
        }
        res[i] = Y[i] - fitval[i];
    }
    
    for(i = 0; i < n; i++){
        for(j = i; j < n; j++){
            for(k = 0; k < p; k++){
                tmpi[k] = X(i, k);
                tmpj[k] = X(j, k);
            }
            if(res[i] * res[j] < 0){
                sum1 += ENorm(tmpi, tmpj);
            }else if((res[i] < 0) & (res[j] < 0)){
                sum2 += ENorm(tmpi, tmpj);
            }else{
                sum3 += ENorm(tmpi, tmpj);
            }
        }
    }
    tstat = sum1 / (n * n * tau * (1 - tau)) - sum2 / (pow(n * tau, 2)) - sum3 / (pow(n * (1 - tau), 2));
    tstat = tstat * (n * tau * (1 - tau));
    
    return tstat;
}