library(Rcpp)
lambdaoo<-function(kmctime,delta,lambda,gtmat){
.Call('kmcomegalambda', PACKAGE = 'kmc',kmctime,delta,lambda,gtmat)
}
