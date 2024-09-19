# kmc.solve2 <- function(x, d, g, em.boost = T, using.num = T, using.Fortran = T, using.C = F, tmp.tag = T, rtol = 1E-9, control = list(nr.it = 20, nr.c = 1, em.it = 3), ...) {}

x <- c(1, 1.5, 2, 3, 4.2, 5.0, 6.1, 5.3, 4.5, 0.9, 2.1, 4.3)
d <- c(1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1)
f <- function(x) {
  x - 3.7
} # \sum f(ti) wi ~ 0
g <- list(f = f)
# define constraint as a list
kmc.solve(x, d, g, using.C = TRUE)
# using Rcpp
#################
# dim =2
#################
myfun5 <- function(x) {
  x^2 - 16.5
}
g <- list(f1 = f, f2 = myfun5)
# define constraint as a list
re0 <- kmc.solve(x, d, g)
control <- list(nr.it = 20, nr.c = 1, em.it = 3)


#' The `kmc.solve2` function calculates the omega and lambda
#' for the kmc data.
#' @param x: time
#' @param d {0,1} indicator of observer/censored
#' @param g the constraints
#' @param rtol the tolerance
#' @param control the control parameters
#' @param ... other parameters
#' @return a list of the results
kmc.solve2 <- function(
    x,
    d,
    g,
    rtol = 1E-9,
    control = list(nr.it = 20, nr.c = 1, em.it = 3), ...) {
  ###### checking PHASE 1         ######
  if (length(unique(d)) != 1) {
    if (!setequal(unique(d), c(0, 1))) stop("Status must be 0/1")
  } else {
    if (d[1] != 1) stop("Status must be 0/1")
  } # check d = 0/1 or all 1's
  if (sum(d) < length(g)) {
    stop("Number of observation MUST be greater than numbers of constraints")
  }
  if ("nr.it" %in% names(control)) {
    nr.it <- control$nr.it
    if (nr.it < 10) nr.it <- 10
  } else {
    nr.it <- 20
  } # NR iteration
  if ("nr.c" %in% names(control)) {
    nr.c <- control$nr.c
    if (nr.c > 1) {
      warning("In N-R iteration, C should be between 0 and 1")
      nr.c <- 1
    }
  } else {
    nr.c <- 1
  } # NR scaler
  if ("em.it" %in% names(control)) {
    em.it <- control$em.it
    if (em.it > 10) em.it <- 10
  } else {
    em.it <- 3
  } # EM iteration
  experimental <- F
  if ("experimental" %in% names(control)) {
    experimental <- T
  }
  default.init <- rep(0., length(g))
  if ("default.init" %in% names(control)) {
    default.init <- control[["default.init"]]
  }
  ###### end of checking PHASE 1  ######


  #' The `omega.lambda_naive` function calculates the omega and lambda
  #' for the kmc data.
  #' @param kmc.time: time
  #' @param delta {0,1} indicator of observer/censored
  #' @param lambda the weight for the constraints
  #' @param g the constraints
  #' @param gt.mat the constraints matrix
  omega.lambda_naive <- function(kmc.time, delta, lambda, g, gt.mat) {
    # iter
    p <- length(g) # the number of constraint
    n <- length(kmc.time)
    if (p == 1) {
      ### 1 constraint: offered by Dr Zhou.
      n <- length(delta)
      delta[n] <- 1
      u.omega <- rep(0, n) ##### all be zero to begin
      u.omega[1] <- 1 / (n - lambda * gt.mat[1]) ## first entry
      S <- rep(1.0, n)
      S.cen <- 0
      for (k in 2:n) {
        if (delta[k] == 0) {
          S[k] <- S[k - 1] - u.omega[k - 1]
          S.cen <- S.cen + 1 / S[k]
        } else {
          u.omega[k] <- 1 / (n - lambda * gt.mat[k] - S.cen)
          S[k] <- S[k - 1] - u.omega[k - 1]
        }
      }
      # return(list(S=S,omega=u.omega, mea= sum(gt*u.omega)))
      return(list(S = S, omega = u.omega, gt = gt.mat))
    } else {
      uncen.loc <- which(delta == 1)
      cen.loc <- which(delta == 0)
      delta[n] <- 1
      #########################################
      u.omega <- numeric(n)
      u.omega[1] <- 1 / (n - sum(lambda * gt.mat[, 1]))
      for (k in 2:n) {
        if (delta[k] == 1) {
          S <- 1 - cumsum(u.omega) # need to update every kmc.time(add in one entry each kmc.time)
          SCenLoc <- cen.loc[cen.loc %in% (1:(k - 1))]
          S.cen <- 0
          if (length(SCenLoc) != 0) {
            S.cen <- sum(1 / S[SCenLoc])
          }
          u.omega[k] <- 1 / (n - sum(lambda * gt.mat[, k]) - S.cen)
          # cat(':::',sum(omega),'\n')
        }
      }
      return(list(S = S, omega = u.omega, gt = gt.mat))
    }
  }


  re <- kmc.clean(kmc.time = x, delta = d)
  kmc.time <- re$kmc.time
  delta <- re$delta

  p <- length(g)
  ## by default, we assume the leading p observations are uncensored
  if (sum(delta == 1) < p) {
    warning("Number of uncensored observations must be greater than the number of constraints")
  }
  delta[1:p] <- 1
  n <- length(delta)
  gt.mat <- matrix(0, p, n)
  for (i in 1:p) gt.mat[i, ] <- g[[i]](kmc.time)

  init.lam <- rep(0, length(g))

  ## compute the NULL hypothesis
  re0 <- omega.lambda_naive(
    kmc.time = kmc.time,
    delta = delta, lambda = init.lam,
    g = g, gt.mat = gt.mat
  )
  loglik.null <- kmc.el(delta, re0$omega, re0$S)

  ## compute the alternative hypothesis

  kmc.comb123 <- function(x) {
    return(kmc_routine4(lambda = x, delta = delta, gtmat = gt.mat))
  }

  mini_iter <- 10
  iter <- 0
  lambda <- init.lam

  while (iter < mini_iter) {
    result_h1 <- multiroot(kmc.comb123, start = lambda, ctol = rtol, maxiter = 512)
    lambda <- result_h1$root
    iter <- iter + result_h1$iter
    rtol <- rtol / 2
  }

  re1 <- omega.lambda_naive(
    kmc.time = kmc.time,
    delta = delta,
    lambda = lambda,
    g = g, gt.mat = gt.mat
  )
  loglik.ha <- kmc.el(delta, re1$omega, re1$S)

  re.tmp <- list(
    loglik.null = loglik.null,
    loglik.h0 = loglik.ha,
    "-2LLR" = 2. * (loglik.null - loglik.ha),
    g = g,
    time = x,
    status = d,
    phat = re1$omega,
    pvalue = 1 - pchisq(-2 * (loglik.ha - loglik.null), df = length(g)),
    lambda = lambda
  )

  class(re.tmp) <- "kmcS3"

  return(re.tmp)
}
