library(Rcpp)
lambdaoo<-function(kmctime,delta,lambda,gtmat){
.Call('kmcomegalambda', PACKAGE = 'kmc',kmctime,delta,lambda,gtmat)
}

kmcdata_rcpp<-function(kmctime,delta,lambda,gtmat){
    .Call('kmcRCPP_KMCDATA', PACKAGE = 'kmc',kmctime,delta,lambda,gtmat)
}