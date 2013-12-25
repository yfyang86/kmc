#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;


double sum(NumericVector y){
double re=0.;
for (int i=0; i<y.size();i++) re+=y(i);
return re;
}

double sum(NumericVector x,NumericVector y){
double re=0.;
for (int i=0; i<y.size();i++) re+=y(i)*x(i);
return re;
}


double sum(double x,NumericMatrix mat,int col){
double re=0.;
for (int i=0; i<mat.nrow();i++) re+=x*mat(i,col);
return re;
}


List omegalambda(SEXP kmctime,SEXP delta,SEXP lambda,SEXP gtmat){
Environment stats("package:stats");
RNGScope scope;
NumericMatrix Gtmat(gtmat);
NumericVector Kmctime(kmctime),Delta(delta);
double Lambda=Rcpp::as<double>(lambda);

//int p=Gtmat.nrow();
int n=Gtmat.ncol();
double tmp=0.;
NumericVector S(n);
//NumericVector uncenloc(floor(sum(Delta)));
    int cenlocL=n-floor(sum(Delta));
NumericVector cenloc(cenlocL);
    int intflg=0;
for (int i=0;i<n;i++){
	if (Delta(i)<.5) {
        cenloc(intflg)=i;//==0
        intflg++;
    }
}
Delta(n-1)=1;//setting the last to be observed!

//start iteration
NumericVector uomega(n);
uomega(0)=1./((double)n-sum(Lambda,Gtmat,0));
double Scen=0.;
    
for (int k=1;k<n;k++){
	if (Delta(k)>.5){//==1
		tmp=0.;
		for (int i=0;i<n;i++){
			tmp+=uomega(i);
			S(i)=1-tmp;
		}
		Scen=0;
		if (cenloc(0)<(k-1)) {
			intflg=0;
			while (intflg<cenlocL  && cenloc(intflg)<=(k-1) ){
				//if (k>47990) printf("Check point: %d iteration %dINTFLG\n",k,intflg);
                Scen+=1/S(floor(cenloc(intflg)));
				intflg++;
                
			}
		}
	 uomega(k)=1./((double)n-sum(Lambda,Gtmat,k) -Scen);
	}
}
List re; re("S")=S;re("omega")=uomega;re("gt")=Gtmat;
return(re);
}
