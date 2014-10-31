#ifndef KMC_COMM__
#define KMC_COMM__

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


double sum(double x,NumericMatrix mat,int col){// <x,mat[,col]>
double re=0.;
for (int i=0; i<mat.nrow();i++) re+=x*mat(i,col);
return re;
}

double sum(double x,double * mat,int col,int p){// <x,mat[,col]>
    double re=0.;
    for (int i=0; i<p;i++) re+=x*mat[p*col+i];
    return re;
}



#endif
