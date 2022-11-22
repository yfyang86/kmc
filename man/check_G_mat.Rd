\name{check_G_mat}
\alias{check_G_mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check the contraints of KMC
}
\description{
To derive the empirical likelihood with constraints, we need to make sure there are solutions.
Dines' method is used here to check whether the linear constraintsare proper or not.
}
\usage{
check_G_mat(gmat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gmat}{
      A p by n. Here p is the number of constraints, n is the number of observations. The matrix is defined in <doi: 10.1201/b18598>.
}
}

\value{
\item{flg}{
A flag:
- 0: not proper
- 1: proper
}
}

\references{
Dines, L. L.  (1926). On positive solutions of a system of linear equations 
\emph{Annals of Mathematics} 
\bold{pages}  386–392

Zhou, M. and Yang, Y. (2015).
A recursive formula for the Kaplan-Meier estimator with mean constraints and its application to empirical likelihood
\emph{Computational Statistics}. \bold{Online ISSN} 1613-9658.
}

\author{
Yifan Yang(yfyang.86@hotmail.com)
}

\examples{
#### A Proper Example ####
x <- c( 1, 1.5, 2, 3, 4.2, 5.0, 6.1, 5.3, 4.5, 0.9, 2.1, 4.3) 
d <- c( 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1) 
f1 <-function(x) { x - 3.7} 
f2 <- function(x) {x^2 - 16.5 } 
g <- list(f1, f2)
re = kmc.clean(x, d)
p = length(g)
n = length(re$kmc.time)
gmat<-matrix(0, p, n);
for(i in 1:p){
  gmat[i,] = g[[i]](re$kmc.time)
  }
# You may want to require(Rcpp) on some platforms (such Mac OSX-ARM)
# library(Rcpp)
# check_G_mat(gmat)
}


