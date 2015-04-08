#ifndef KMC_COMM__
#define KMC_COMM__

#include <R.h>
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
#include <vector>
using namespace Rcpp;
using namespace std;

double sum(NumericVector y){ // \sum y_i
double re=0.;
for (int i=0; i<y.size();i++) re+=y(i);
return re;
}

double sum(NumericVector x,NumericVector y){//<x,y>
double re=0.;
for (int i=0; i<y.size();i++) re+=y(i)*x(i);
return re;
}


double sum(vector<double> x,NumericMatrix mat,int col){// <x,mat[,col]>
double re=0.;
for (int i=0; i<mat.nrow();i++) re+=x[i]*mat(i,col);
return re;
}

double sum(double *x,double * mat,int col,int p){// <x,mat[,col]>
    double re=0.;
    for (int i=0; i<p;i++) re+=x[i]*mat[p*col+i];
    return re;
}



#endif
