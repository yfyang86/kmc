\name{kmc.bcm.test}
\alias{kmc.bcm.test}
\title{
Empirical-likelihood ratio test of beta in the binary-choice /
current-status (Case-1 interval-censored) model.
}
\description{
Tests the regression coefficient \eqn{\beta} in the model
\deqn{y_i = \mathbf{1}\{\beta^\top x_i + \varepsilon_i > 0\}, \quad
      \delta_i = 1 - y_i,}
with i.i.d.\ residual law \eqn{F_\varepsilon} treated
nonparametrically. The test uses Owen's empirical-likelihood ratio
applied to a Buckley-James-style estimating function built on top of
the PAVA NPMLE of \eqn{F_\varepsilon}. With the scale-identification
convention that one coordinate of \eqn{\beta} is fixed by dropping
the corresponding component of the estimating function,
\deqn{
  -2\log\mathrm{ELR}(\beta_0)\;\xrightarrow{d}\;\chi^{2}_{p-1}
}
under regularity (see references). Two variants of the estimating
function are available:
\itemize{
\item \code{centered = FALSE} (default): raw Buckley-James
estimating function \eqn{g_i = x_i\,\hat\varepsilon_i}. Finite-sample
\eqn{\chi^2} calibration relies on a PAVA-induced cancellation
discussed in the references.
\item \code{centered = TRUE}: centered estimating function
\eqn{\tilde g_i = (x_i - \hat m(e_i))\,\hat\varepsilon_i} where
\eqn{\hat m} is a Nadaraya-Watson kernel regression of \eqn{x} on
the residual axis \eqn{e_i = -\beta^\top x_i}. The centering
discharges the no-bias condition automatically via the tower
property; this is the rigorous variant.
}
}
\usage{
kmc.bcm.test(X, delta, beta, centered = FALSE, h = NULL, drop = 1L)
}
\arguments{
  \item{X}{
  \eqn{n\times p} matrix of covariates. The identification convention
  is to include no intercept and to fix one coordinate of \eqn{\beta}
  (the \code{drop}-th by default) to a known value.
  }
  \item{delta}{
  Integer 0/1 vector of length \eqn{n}. In BCM notation,
  \eqn{\delta_i = 1} when the latent event has occurred at or before
  the observation time (equivalently, \eqn{\varepsilon_i \le -\beta^\top x_i}).
  }
  \item{beta}{
  Numeric vector of length \eqn{p} giving the regression coefficient
  to be tested. Component \code{beta[drop]} should equal its
  scale-fixing value.
  }
  \item{centered}{
  Logical. If \code{TRUE}, use the centered estimating function with
  Nadaraya-Watson kernel correction; if \code{FALSE} (default), use
  the raw Buckley-James estimating function.
  }
  \item{h}{
  Bandwidth for the Nadaraya-Watson regression (centered case only).
  If \code{NULL} (default), uses
  \eqn{h_n = 1.5\,\hat\sigma_e\,n^{-1/3}} (under-smoothed
  relative to Horowitz's MSE-optimal \eqn{n^{-1/5}} to ensure
  asymptotic bias control).
  }
  \item{drop}{
  Integer in \code{1:p}. Index of the coordinate of \eqn{\beta} that
  is fixed by the scale identification. The corresponding column of
  the estimating function matrix is removed before Owen-EL.
  Default \code{1L}.
  }
}
\value{
A list of class \code{kmcS3} with components:
\item{-2LLR}{
The \eqn{-2\log} empirical-likelihood ratio; asymptotically
\eqn{\chi^2_{p-1}} under \eqn{H_0:\beta=\beta_0}.
}
\item{df}{
Degrees of freedom, equal to \eqn{p-1}.
}
\item{pvalue}{
The asymptotic p-value computed against \eqn{\chi^2_{p-1}}.
}
\item{centered}{Whether the centered variant was used.}
\item{drop}{Which coordinate was dropped for scale identification.}
\item{h}{Bandwidth used (NA if \code{centered=FALSE}).}
\item{lambda}{Dual variable from \code{emplik::el.test}.}
\item{convergence}{1 if the EL inner solver succeeded, 0 otherwise.}
}
\references{
Owen, A. (2001). \emph{Empirical Likelihood}. Chapman and Hall/CRC.

Cosslett, S. R. (1983). Distribution-free maximum likelihood
estimator of the binary choice model. \emph{Econometrica}, 51(3),
765--782.

Murphy, S. A., and van der Vaart, A. W. (1997). Semiparametric
likelihood ratio inference. \emph{Annals of Statistics}, 25(4),
1471--1509.

Yang, Y., and Zhou, M. (2026). \emph{Empirical Likelihood under
Censoring: Computational Methods for Survival and Discrete-Choice
Data}. (Chapter 3, Theorems 3.1 and 3.10).
}
\author{
Yifan Yang \email{yfyang.86@hotmail.com}, Mai Zhou
}
\seealso{
\code{\link{kmc.bjtest}} for the analogous test in the right-censored
AFT model; \code{\link{kmc.solve}} for the constrained KM estimator.
}
\examples{
## Simulated binary-choice data with two continuous covariates,
## logistic residuals, and true beta = c(1, 0.5):
set.seed(42)
n     <- 800
beta0 <- c(1.0, 0.5)
X     <- matrix(runif(n * 2, -1, 1), n, 2)
eps   <- rlogis(n)
delta <- as.integer(as.numeric(X \%*\% beta0) + eps <= 0)

## Raw BJ-residual EL test:
out_raw <- kmc.bcm.test(X, delta, beta0, centered = FALSE)
out_raw$`-2LLR`
out_raw$pvalue

## Centered EL test (rigorous variant):
out_cn  <- kmc.bcm.test(X, delta, beta0, centered = TRUE)
out_cn$`-2LLR`
out_cn$pvalue

## Compare with a deliberately wrong beta -- the test should reject:
out_wrong <- kmc.bcm.test(X, delta, beta = c(1, -1.5), centered = FALSE)
out_wrong$`-2LLR`     # large
out_wrong$pvalue      # small
}
\keyword{survival}
\keyword{nonparametric}
