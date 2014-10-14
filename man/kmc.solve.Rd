\name{kmc.solve}
\alias{kmc.solve}
\title{
Calculate NPMLE with constriants for right censored data
    %%  ~~function to do ... ~~
}
\description{
This function calculate the Kaplan-Meier estimator with mean constraints recursively. 
\deqn{El(F)=\prod_{i=1}{n}(\Delta F(T_i))^{\delta_i}(1-F(T_i))^{1-\delta_i}}
with constraints
\deqn{\sum_i g(T_i)\Delta F(T_i)=0,\quad,i=1,2,...}
It uses Lagrange multiplier directly.
}
\usage{
kmc.solve(x, d, g, em.boost = T, using.num = T, using.Fortran =
                 T, using.C = F, tmp.tag = T, rtol = 1e-09, control =
                 list(nr.it = 20, nr.c = 1, em.it = 3),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
      Possitive time
}
  \item{d}{
      Status, 0: right censored; 1 uncensored
}
  \item{g}{
list of contraint functions. It should be a list of functions list(f1,f2,...) 
}
  \item{em.boost}{
  	 logical asking whether to use EM to get the initial value, default=TRUE. See 'Details' for EM control.
}
\item{using.num}{
  	 logical asking whether to use numeric derivative in iterations, default=TRUE.
}
\item{using.Fortran}{
  	 logical asking whether to use Fortran in root solving, default=F.
}
\item{using.C}{
	logical asking whether to use Rcpp in iteraruib, default=T. This option will promote the performance of KMC algorithm. Development version works on one constraint only, otherwise it will generate a Error information. It won't work on using.num=F.
}
\item{tmp.tag}{
	Development version needs it, keep it as TRUE.
}
\item{rtol}{
	Tolerance used in rootSolve(multiroot) package, see 'rootSolve::multiroot'.
}
\item{control}{
  nr.it controls max iterations allowed in N-R algorithm default=20, nr.c is the scaler used in N-R algorithm default=1,em.it is max iteration if use EM algorithm (em.boost) to get the initial value of lambda, default=3.
}
\item{...}{
Unspecified yet.
}
}
\details{
The development version won't check the constraint is proper or not. Besides, a extrordinary large number will be generated when numeric method/constraint(s) fails.
}
\value{
 a list with the following components:
 \item{loglik.ha}{
 	log empirical  likelihood without constraints
 }
 \item{loglik.h0}{
 	log empirical likelihood with constraints
 }
 \item{"-2llr"}{
 	-2 Log empirical likelihood
 }
 \item{phat}{
 	\deqn{\Delta F(T_i)}
 }
 \item{pvalue}{
  p-value of the test
 }
 \item{df}{
  Degree(s) of freedom. It equals the number of constraints.
 }
 \item{lambda}{
  The lambda is the Lagrangian multiplier described in reference.
 }
}
\references{
    A Recursive Formula for the Kaplan-Meier Estimator with Constraint
%% ~put references to the literature/web site here ~
}
\author{
    Mai Zhou(mai@ms.uky.edu), Yifan Yang(yifan.yang@uky.edu)
%%  ~~who you are~~
}


\seealso{
	  	\code{\link{plotkmc2D}}, \code{\link{print.kmcS3}}.
}
\examples{
x <- c( 1, 1.5, 2, 3, 4.2, 5.0, 6.1, 5.3, 4.5, 0.9, 2.1, 4.3) # positive time
d <- c( 1,   1, 0, 1, 0, 1, 1, 1, 1, 0, 0,   1)               # status censored/uncensored

#################
# dim =1
#################

f<-function(x) { x-3.7}                                       # \sum f(ti) wi ~ 0 
g=list( f=f) ;                                                #define constraint as a list

kmc.solve( x,d,g) ;                                           #using default
kmc.solve( x,d,g,using.C=TRUE) ;                              #using Rcpp

#################
# dim =2
#################

myfun5 <- function( x)  { 
 x^2-16.5
} 
g=list( f1=f,f2=myfun5) ;                                     #define constraint as a list
kmc.solve( x,d,g) ->re0;
print(re0);                                                   #print method, plain ASCII
print(re0,5);                                                 #print method, digits=5,plain ASCII
print(re0,type='md');                                         #print method, output GIT Markdown
#summary(re0) is under developing.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
