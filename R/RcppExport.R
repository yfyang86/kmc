check_G_mat <- function(gmat){
    .Call('kmcRCPP_RevCHECK',PACKAGE='kmc', gmat)
}


lambdaoo<-function(kmctime,delta,lambda,gtmat){
.Call('kmcomegalambda', PACKAGE = 'kmc', kmctime,delta,lambda,gtmat)
}

kmcdata_rcpp<-function(kmctime,delta,lambda,gtmat){
    .Call('kmcRCPP_KMCDATA', PACKAGE = 'kmc', kmctime,delta,lambda,gtmat)
}

kmc_routine4<-function( 
		delta, # status 
		lambda,# root finding 
		gtmat # g(t) Matrix p X n 
){ 
    np=as.integer(dim(gtmat)) 
    chk=numeric(np[1]) 
    delta=as.integer(delta) 
    lambda=as.vector(lambda); 
    uomega=rep(0, np[2])
    # gtmat=as.numeric(gtmat) 
    # .C extension will ignore variable name! 
    re=.C( 
		"nocopy_kmc_data", delta, as.vector(lambda%*%gtmat),uomega,np,chk 
    ) 
    return(as.vector(gtmat %*% re[[3]]));
}

kmc_routine5 <- function(
    		delta, # status 
		    lambda,# root finding 
		    gtmat # g(t) Matrix p X n 
            ){
                np=as.integer(dim(gtmat))             
                w = rep(0., np[2])
                delta=as.integer(delta) 
                lambda=as.vector(lambda); 
                lambagt=as.vector(lambda%*%gtmat)
                re=.C("kmc_native", delta ,lambagt ,w , np)
                return(as.vector(re[[3]]));
}

kmc_routine5_1d <- Vectorize(function(lambda, delta, time, flist){
    sum(kmc_routine5(
      delta  = delta, 
      lambda = lambda, 
      gtmat  = matrix(
                unlist(lapply(1:length(flist), function(i) flist[[i]] (time)))
                , nrow=length(flist) 
                , byrow = T) )) - 1.}, vectorize.args = "lambda")

kmc_routine5_nd <- function(lambda, delta, time, flist){
  sum(kmc_routine5(
    delta  = delta, 
    lambda = lambda, 
    gtmat  = matrix(
      unlist(lapply(1:length(flist), function(i) flist[[i]] (time)))
      , nrow=length(flist) 
      , byrow = T) )) - 1.}

kmc_find0loc <- function(d){
 re=1L;
 d=as.integer(d);
 n=as.integer(length(d));
 return(.C("locLastZero",d,n,re)[[3]]);
}
