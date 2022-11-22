\name{kmc.clean}
\alias{kmc.clean}

\title{
Perform Data Clean for the kmc Algorithm
}
\description{
The kmc.clean function clean the (kmc.time, delta) for the randomized censored data:

- Reorder the data according to the observed time and status;

- Clean the (right) censored data point(s) if they happen before the first uncesored data.

- If there are ties in the data. For the time points contain ties, e.g. \deqn{(T_{i_s}, d_{i_s}), i_s \in S \forall j \in S, T_{j} \equiv T}, we re-arranged the data in a manner that those with d=1 are ordered ahead of those with d=0. As d=0 indicates the data point is right censored, such procedure is trivial. }

\usage{
kmc.clean(kmc.time, delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kmc.time}{
     Non-negative real vector. The observed time.
}
  \item{delta}{
      0/1 vector. Censoring status indictator, 0: right censored; 1 uncensored
}
}

\value{
 A list with the following components:
 \item{kmc.time}{
 	The cleaned observed time.
 }
 \item{delta}{
 	The cleaned censoring status indictator, 0: right censored; 1 uncensored
 }
}

\references{
Zhou, M. and Yang, Y. (2015).
A recursive formula for the Kaplan-Meier estimator with mean constraints and its application to empirical likelihood
\emph{Computational Statistics}
\bold{Online ISSN} 1613-9658.
}

\author{
Yifan Yang(yfyang.86@hotmail.com)
}

\examples{
x <- c( 1, 1.5, 2, 3, 4.2, 5.0, 6.1, 5.3, 4.5, 0.9, 2.1, 4.3) 
d <- c( 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1) 
kmc.clean(x, d)
}


