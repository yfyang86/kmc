
#### KMC Nov ####
library(Rcpp)
library(rootSolve)

src = "NumericVector kmc_native(NumericVector delta, NumericVector lambda_gt){
  int n = delta.size();
  NumericVector S(n);
  NumericVector w(n);
  size_t i = 0;
  double n_double = (double) n;
  for (size_t j =0;j < n; j++) w[j] = 0.;
  w[i] = 1/(n_double - lambda_gt[i]);
  double sumSjDelta0 = 0;
  S[i] = 1. -  w[i];
  while (i < n - 1){

    if (delta[i] < 0.5) {
      sumSjDelta0 += 1/S[i];
    }

    w[i+1] =  (delta[i+1]>0.0? 1./(n_double - lambda_gt[i+1] - sumSjDelta0):0.);
    S[i+1] = S[i] - w[i+1];
    i++;
  }
  //w[n-1] = 1.;
  //for (size_t j=0; j< n-1; j++) w[n-1] -=  w[j];
  w[n-1] = ( w[n-1] < 0 ? 0:w[n-1]);
  return w;
}
"

cppFunction(src)

kmc.clean <- function(kmc.time, delta){
  n <- length(kmc.time)
  dataOrder <- order(kmc.time, -delta)
  kmc.time <- kmc.time[dataOrder]
  delta <- delta[dataOrder]             ### changed 10/2018

  FirstUnCenLocation<-which(delta==1)[1];
  if (FirstUnCenLocation==n) {stop('Only one uncensored point.');}
  if (FirstUnCenLocation!=1){
    delta=delta[FirstUnCenLocation:n];
    kmc.time=kmc.time[FirstUnCenLocation:n];
  }
  delta[length(kmc.time)]=1;
  return (list(kmc.time=kmc.time,delta=delta));
}

kmc.solver_rec <- function(x, d, glist, init = 0){
  inputData = kmc.clean(x, d)
  p = length(glist)
  n = length(inputData[[1]])

  Gmat = matrix(0, p, n)
  for(i in 1:p) Gmat[i, ] = glist[[i]](inputData$kmc.time)

  fun.C <- function(lam){
      w = kmc_native(inputData$delta, lam%*%Gmat)
      return((Gmat%*%w));
    }

if (p==1){
 fun.C2 <- function(lam){
      w = kmc_native(inputData$delta, lam%*%Gmat)
      return(1-sum((w)));
    }
    Cmean = multiroot(start = -2+init, f = fun.C2 )
}else{
  Cmean = NULL
}

Cmu = multiroot(start = init, f = fun.C )
return(list(Cmu, Cmean));
}


#### Example ####

x <- c( 1, 1.5, 2, 3, 4.2, 5.0, 6.1, 5.3, 4.5, 0.9, 2.1, 4.3) # positive time
d <- c( 1,   1, 0, 1, 0, 1, 1, 1, 1, 0, 0,   1)               # status censored/uncensored

     
f<-function(x) { x-3.7}                                       # \sum f(ti) wi ~ 0 

myfun5 <- function( x)  { 
      x^2-16.5
     } 
     
g = list( f1=f,f2=myfun5) ;                                    # define constraint as a list

# dim == 2
kmc.solver_rec(x, d, g, c(1,1))

# dim == 1 test 
kmc.solver_rec(x, d, list(f), 1)

kmc.solver_rec(x, d, list(f), -1)


emplik::el.cen.EM2( x  = x, d = d, fun = f, mu = 0, maxit = 50)