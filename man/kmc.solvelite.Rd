\name{kmc.solve}
\alias{kmc.solve}
\title{
Calculate NPMLE with constriants for right censored data
    %%  ~~function to do ... ~~
}
\description{
This function calculate the Kaplan-Meier estimator with mean constraints recursively. 
\deqn{El(F)=\prod_{i=1}^{n}(\Delta F(T_i))^{\delta_i}(1-F(T_i))^{1-\delta_i}}
with constraints
\deqn{\sum_i g(T_i)\Delta F(T_i)=0,\quad,i=1,2,...}
It uses Lagrange multiplier directly. This function is a lite version of kmc.solve.
}
\usage{
kmc.solve(x, d, g, rtol = 1e-09, 
          control = list(nr.it = 20, nr.c = 1, em.it = 3),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
     Non-negative real vector. The observed time.
}
  \item{d}{
      0/1 vector. Censoring status indictator, 0: right censored; 1 uncensored
}
  \item{g}{
list of contraint functions. It should be a list of functions list(f1,f2,...) 
}

\item{rtol}{
	Tolerance used in rootSolve(multiroot) package, see 'rootSolve::multiroot'.
}
\item{control}{
  A list. The entry nr.it controls max iterations allowed in N-R algorithm default=20; nr.c is the scaler used in N-R algorithm default=1; em.it is max iteration if use EM algorithm (em.boost) to get the initial value of lambda, default=3.
}
\item{...}{
Unspecified yet.
}
}
\details{
The function check_G_mat checks whether the solution space is null or not under the constraint. But due to the computational complexity, it will detect at most two conditions.
}
\value{
 A list with the following components:
 \item{loglik.ha}{
 	The log empirical  likelihood without constraints
 }
 \item{loglik.h0}{
 	The log empirical likelihood with constraints
 }
 \item{"-2LLR"}{
 	The -2 Log empirical likelihood ratio
 }
 \item{phat}{
 	\deqn{\Delta F(T_i)}
 }
 \item{pvalue}{
  The p-value of the test
 }
 \item{df}{
  Degree(s) of freedom. It equals the number of constraints.
 }
 \item{lambda}{
  The lambda is the Lagrangian multiplier described in reference.
 }
}
\references{
Zhou, M. and Yang, Y. (2015).
A recursive formula for the Kaplan-Meier estimator with mean constraints and its application to empirical likelihood
\emph{Computational Statistics}
\bold{Online ISSN} 1613-9658.
}

\author{
    Mai Zhou(mai@ms.uky.edu), Yifan Yang(yfyang.86@hotmail.com)
}


\seealso{
	  	\code{\link{plotkmc2D}}.
}
\examples{
# positive time
x <- c( 1, 1.5, 2, 3, 4.2, 5.0, 6.1, 5.3, 4.5, 0.9, 2.1, 4.3) 
# status censored/uncensored
d <- c( 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1)               

#################
# dim =1
#################

f <- function(x) {x-3.7}                  # \sum f(ti) wi ~ 0 
g <- list(f=f)                            #define constraint as a list

kmc.solve(x, d, g)                        #using default
kmc.solvelite(x, d, g)                    #using default

#################
# dim =2
#################

myfun5 <- function( x)  { 
 x^2-16.5
} 
g <- list(f1=f, f2=myfun5) ;              #define constraint as a list
re0 <- kmc.solvelite(x, d, g);

###################################################
# Real Data: Stanford Heart Transplant
###################################################
# library(survival)
# stanford5 <- stanford2[!is.na(stanford2$t5), ]
# y=log10(stanford5$time)
# d <- stanford5$status
# g <- list(f = function(x) {(x-2.4 < 0.) - 0.5} )  # \sum I(x_i<2.4) wi = 0.5)
# kmc.solvelite(y, d, g)
##################################
## Example (Cont'd): -2LLR Curve
##################################
# iters = 100
# scale = 2.
# starter = 2.5
# result = rep(0., 20)
# observe_range = starter + (1:iters)/(iters * scale)
# for (i in 1:iters){
#         g <- list(f = function(x) {(x- (starter + i/(iters * scale)) < 0.) - 0.5})
#         result[i] = kmc.solvelite(y, d, g)$`-2LLR`
# }
# plot(x = observe_range, result, xlab = 'time', ylab = '-2LLR', type = 'b')
}