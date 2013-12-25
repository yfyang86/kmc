#include <Rcpp.h>

using namespace Rcpp;

// rcpp_hello_world
List	omegalambda(SEXP kmctime,SEXP delta,SEXP lambda,SEXP gtmat);

RcppExport SEXP kmcomegalambda(SEXP kmctime,SEXP delta,SEXP lambda,SEXP gtmat) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        List __result = omegalambda( kmctime,delta,lambda,gtmat);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

