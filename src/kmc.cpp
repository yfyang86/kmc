#include "common.h"
#include <vector>
#include <Rcpp.h>
#include <iostream>

// Check the sign of a matrix
List RevCHECK(SEXP xx) {
    NumericMatrix x(xx);
    List re;
    NumericVector rr(1);
    
    if (x.nrow() != 2) {
        rr(0) = signcheck(x);
    } else {
        std::vector<int> pos_in;
        std::vector<int> neg_in;
        
        // Separate positive and negative indices
        for (int j = 0; j < x.ncol(); j++) {
            if (x(0, j) > 0) pos_in.push_back(j);
            if (x(0, j) < 0) neg_in.push_back(j);
        }
        
        NumericMatrix Amat(1, pos_in.size() * neg_in.size());
        int ind = 0;
        
        // Compute Amat values
        for (int i = 0; i < pos_in.size(); i++) {
            for (int j = 0; j < neg_in.size(); j++) {
                Amat(0, ind) = x(0, pos_in[i]) * x(1, neg_in[j]) - x(0, neg_in[j]) * x(1, pos_in[i]);
                ind++;
            }
        }
        
        rr(0) = signcheck(Amat);
    }
    
    re("flg") = rr;
    return re;
}

/*
* Compute omega and lambda values
*
* @param kmctime: T
* @param delta:   status
* @param lambda:  root finding
* @param gtmat:   g_1(X),...,g_p(X)
*/
List omegalambda(SEXP kmctime, SEXP delta, SEXP lambda, SEXP gtmat) {
    Environment stats("package:stats");
    RNGScope scope;
    
    NumericMatrix Gtmat(gtmat);
    NumericVector Kmctime(kmctime), Delta(delta);
    vector<double> Lambda = Rcpp::as<std::vector<double>>(lambda);
    
    int n = Gtmat.ncol();
    NumericVector S(n);
    int cenlocL = n - floor(sum(Delta));
    NumericVector cenloc(cenlocL);
    int intflg = 0;
    
    // Find censored locations
    for (int i = 0; i < n; i++) {
        if (Delta(i) < 0.5) {
            cenloc(intflg) = i;
            intflg++;
        }
    }
    
    Delta(n - 1) = 1; // Setting the last to be observed
    
    // Start iteration
    NumericVector uomega(n);
    uomega(0) = 1.0 / ((double)n - sum(Lambda, Gtmat, 0));
    double Scen = 0.0;
    
    for (int k = 1; k < n; k++) {
        if (Delta(k) > 0.5) { // Observed
            double tmp = 0.0;
            for (int i = 0; i < n; i++) {
                tmp += uomega(i);
                S(i) = 1 - tmp;
            }
            
            Scen = 0.0;
            if (cenloc(0) < (k - 1)) {
                intflg = 0;
                while (intflg < cenlocL && cenloc(intflg) <= (k - 1)) {
                    Scen += 1 / S(floor(cenloc(intflg)));
                    intflg++;
                }
            }
            
            uomega(k) = 1.0 / ((double)n - sum(Lambda, Gtmat, k) - Scen);
        }
    }
    
    List re;
    re("S") = S;
    re("omega") = uomega;
    re("gt") = Gtmat;
    
    return re;
}

/*
* Compute KMC data
*
* @param kmctime: T
* @param delta:   status
* @param lambda:  root finding
* @param gtmat:   g_1(X),...,g_p(X)
*/
List RCPP_KMCDATA(SEXP kmctime, SEXP delta, SEXP lambda, SEXP gtmat) {
    Environment stats("package:stats");
    RNGScope scope;
    
    NumericMatrix Gtmat(gtmat);
    NumericVector Kmctime(kmctime), Delta(delta);
    vector<double> Lambda = Rcpp::as<std::vector<double>>(lambda);
    
    int n = Gtmat.ncol();
    int p__ = Gtmat.nrow();
    NumericVector S(n);
    int cenlocL = n - floor(sum(Delta));
    NumericVector cenloc(cenlocL);
    int intflg = 0;
    
    // Find censored locations
    for (int i = 0; i < n; i++) {
        if (Delta(i) < 0.5) {
            cenloc(intflg) = i;
            intflg++;
        }
    }
    
    Delta(n - 1) = 1; // Setting the last to be observed
    
    // Start iteration
    NumericVector uomega(n);
    uomega(0) = 1.0 / ((double)n - sum(Lambda, Gtmat, 0));
    double Scen = 0.0;
    
    for (int k = 1; k < n; k++) {
        if (Delta(k) > 0.5) { // Observed
            double tmp = 0.0;
            for (int i = 0; i < n; i++) {
                tmp += uomega(i);
                S(i) = 1 - tmp;
            }
            
            Scen = 0.0;
            if (cenloc(0) < (k - 1)) {
                intflg = 0;
                while (intflg < cenlocL && cenloc(intflg) <= (k - 1)) {
                    Scen += 1.0 / S(floor(cenloc(intflg)));
                    intflg++;
                }
            }
            
            uomega(k) = 1.0 / ((double)n - sum(Lambda, Gtmat, k) - Scen);
        }
    }
    
    NumericVector chk(p__, 0.0);
    for (int i = 0; i < p__; i++) {
        for (int j = 0; j < n; j++) {
            chk(i) = chk(i) + Delta(j) * uomega(j) * Gtmat(i, j);
        }
    }
    
    NumericVector Gamma(n);
    for (int i = 0; i < n; i++) {
        Gamma(i) = 1 / S(i);
    }
    
    List re;
    re("omega") = uomega;
    re("gamma") = Gamma;
    re("S") = S;
    re("chk") = chk;
    
    return re;
}

extern "C" {
    /*
    * Compute KMC data without copying
    *
    * @param delta:  status
    * @param gtmat:  g_1(X),...,g_p(X)
    * @param uomega: omega values (output)
    * @param np:     dimensions (p, n)
    * @param chk:    check values (output)
    */
    void nocopy_kmc_data(int* delta, double* gtmat, double* uomega, int* np, double* chk) {
        int nn = np[1];  // col: sample size
        
        vector<double> S(nn);
        int cenlocL = nn;
        for (int i = 0; i < nn; ++i) {
            cenlocL -= delta[i];
        }
        
        vector<int> cenloc(cenlocL); // Index of Censoring
        int intflg = 0;
        for (int i = 0; i < nn; i++) {
            if (delta[i] == 0) {
                cenloc[intflg] = i;
                intflg++;
            }
        }
        
        uomega[0] = 1.0 / ((double)nn - gtmat[0]); // Iteration
        S[0] = 1.0 - uomega[0];
        double Scen = 0.0;
        int k = 0;
        
        while (k < nn - 1) {
            if (delta[k] == 0) {
                Scen += 1 / S[k];
            }
            
            uomega[k + 1] = (delta[k + 1] == 1 ? 1.0 / ((double)nn - gtmat[k + 1] - Scen) : 0.0);
            S[k + 1] = S[k] - uomega[k + 1];
            k++;
        }
        
        uomega[nn - 1] = (uomega[nn - 1] < 0 ? 0 : uomega[nn - 1]);
    }
}

extern "C" {
    /**
    * Compute KMC natively
    *
    * @param delta:     status (double array of length n)
    * @param lambda_gt: <lambda, g(T_i)> values (double array of length n)
    * @param w:         jumps (output, length = sum(delta))
    * @param np:        dimensions (n, p)
    */
    void kmc_native(double* delta, double* lambda_gt, double* w, int* np) {
        int n = np[1];
        std::vector<double> S(n);
        size_t i = 0;
        double n_double = (double)n;
        
        for (size_t j = 0; j < n; j++) {
            w[j] = 0.0;
        }
        
        w[i] = 1 / (n_double - lambda_gt[i]);
        double sumSjDelta0 = 0;
        S[i] = 1.0 - w[i];
        
        while (i < n - 1) {
            // Iterative Computing: KMC
            if (delta[i] < 0.5) {
                sumSjDelta0 += 1 / S[i];
            }
            
            w[i + 1] = (delta[i + 1] > 0.0 ? 1.0 / (n_double - lambda_gt[i + 1] - sumSjDelta0) : 0.0);
            S[i + 1] = S[i] - w[i + 1];
            i++;
        }
        
        #ifdef S_ROUTINE_
        // TODO: Implement S_ROUTINE_
        w[n - 1] = 1.0;
        for (size_t j = 0; j < n - 1; j++) {
            w[n - 1] -= w[j];
        }
        #endif
        
        w[n - 1] = (w[n - 1] < 0 ? 0 : w[n - 1]);
    }
}