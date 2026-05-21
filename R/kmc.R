#' The `kmc.el` function calculates the empirical likelihood
#' for the kmc data.
#' @param delta {0,1} indicator of observer/censored
#' @param omega the weight for the uncensored data
#' @param S the weight for the censored data
kmc.el <- function(delta, omega, S) {
  llog <- function(z, eps) {
    ans <- z
    avoidNA <- !is.na(z)
    lo <- (z < eps) & avoidNA
    ans[lo] <- log(eps) - 1.5 + 2 * z[lo] / eps - 0.5 * (z[lo] / eps)^2
    ans[!lo] <- log(z[!lo])
    ans
  }
  n <- length(S)
  val1 <- omega[delta == 1]
  val2 <- S[delta == 0]
  eps <- 1e-7
  sum(llog(val1, .001 / n^2)[val1 > eps]) + sum(llog(val2, .001 / n^2)[val2 > eps])
}

#' The `kmc.clean` function clean the (kmc.time, delta) for
#' the randomized censored data.
#' 1. No tie: the kmc.time and delta are re-arranged in an increasing order.
#' 2. Tie: for the time points contain ties,
#'            e.g. (T_{i_s}, d_{i_s}), i_s \in S \forall j \in S, T_{j} \equiv T
#'            we re-arranged the data in a manner that those with d=1 are ordered
#'            ahead of those with d=0. As d=0 indicates the data point is right
#'            censored, such procedure is trivial.
#' @param kmc.time: time
#' @param delta {0,1} indicator of observer/censored
#' @param tie method to break the tie, default is "fake"
#' @example
#'
kmc.clean <- function(kmc.time, delta, tie = 'fake') {
  n <- length(kmc.time)
  dataOrder <- order(kmc.time, -delta)
  kmc.time <- kmc.time[dataOrder]
  tie_loc <- which(diff(kmc.time) == 0)
  if (length(tie_loc) > 0 ){
    if (tie == 'fake') {
      kmc.time[tie_loc] <- kmc.time[tie_loc] + seq_along(tie_loc) * (1E-12)
    }
  }
  # make sure the way to break the tie is consistent with the delta's
  delta <- delta[dataOrder]
  FirstUnCenLocation <- which(delta == 1)[1]
  if (FirstUnCenLocation == n) {
    stop("Only one uncensored point.")
  }
  if (FirstUnCenLocation != 1) {
    delta <- delta[FirstUnCenLocation:n]
    kmc.time <- kmc.time[FirstUnCenLocation:n]
  }
  delta[length(kmc.time)] <- 1
  return(list(kmc.time = kmc.time, delta = delta))
}



omega.lambda <- cmpfun(function(kmc.time, delta, lambda, g, gt.mat) {
  n <- length(kmc.time)
  delta[n] <- 1

  # Compute lambda' * g(t_k) for each observation k
  # Works for both scalar lambda (p=1) and vector lambda (p>1)
  lg <- colSums(lambda * gt.mat)

  u.omega <- numeric(n)
  S <- rep(1.0, n)
  S.cen <- 0

  u.omega[1] <- 1 / (n - lg[1])
  for (k in 2:n) {
    S[k] <- S[k - 1] - u.omega[k - 1]
    if (delta[k] == 0) {
      S.cen <- S.cen + 1 / S[k]
    } else {
      u.omega[k] <- 1 / (n - lg[k] - S.cen)
    }
  }

  return(list(S = S, omega = u.omega, gt = gt.mat))
})

kmc.data <- cmpfun(function(kmc.time, delta, lambda, g, gt.mat, using.C = F) {
  # my omega contains 0, it is length of n
  # need S omega
  if (using.C) {
    # tmp<-lambdaoo(kmc.time,delta,lambda,gt.mat);
    return(kmcdata_rcpp(kmctime = kmc.time, delta = delta, lambda = lambda, gtmat = gt.mat))
  } else {
    tmp <- omega.lambda(kmc.time, delta, lambda, g, gt.mat)
  }
  b <- delta * tmp$omega
  # check.constriant=apply(t(b*t(tmp$gt)),1,sum);
  check.constriant <- rowSums(t(b * t(tmp$gt)))
  gama <- 1 / (tmp$S)
  return(list(omega = tmp$omega, gamma = gama, S = tmp$S, chk = check.constriant))
})

omega.lambda12 <- cmpfun(function(kmc.time, delta, lambda, g, gt.mat) {
  # iter
  cumsumx <- function(x) apply(x, 1, cumsum)
  p <- length(g) # the number of constraint
  n <- length(kmc.time)
  uncen.loc <- which(delta == 1)
  cen.loc <- which(delta == 0)
  delta[n] <- 1
  #########################################
  u.omega <- numeric(n)
  udev.omega <- matrix(0, p, n)
  u.omega[1] <- 1 / (n - sum(lambda * gt.mat[, 1]))
  udev.omega[, 1] <- u.omega[1]^2 * gt.mat[, 1]
  for (k in 2:n) {
    if (delta[k] == 1) {
      S <- 1 - cumsum(u.omega) # need to update every kmc.time(add in one entry each kmc.time)
      SCenLoc <- cen.loc[cen.loc %in% (1:(k - 1))]
      S.cen <- 0
      S.cen2 <- 0
      if (length(SCenLoc) != 0) {
        S.cen <- sum(1 / S[SCenLoc])
        S.cen2 <- sum(1 / (S[SCenLoc]^2) * cumsumx(matrix(udev.omega[, 1:(k - 1)], ncol = k - 1))[SCenLoc])
      }
      u.omega[k] <- 1 / (n - sum(lambda * gt.mat[, k]) - S.cen)
      udev.omega[k] <- u.omega[k]^2 * (gt.mat[, k] + S.cen2)
      # cat(':::',sum(omega),'\n')
    }
  }
  return(list(S = S, omega = u.omega, gt = gt.mat, omega.dev = sum(delta * gt.mat * udev.omega)))
})

kmc.data12 <- function(kmc.time, delta, lambda, g, gt.mat) {
  # my omega contains 0, it is length of n
  # need S omega
  tmp <- omega.lambda12(kmc.time, delta, lambda, g, gt.mat)
  b <- delta * tmp$omega
  # check.constriant=apply(t(b*t(tmp$gt)),1,sum);
  check.constriant <- rowSums(t(b * t(tmp$gt)))
  gama <- 1 / (tmp$S)
  return(list(omega = tmp$omega, gamma = tmp$gamma, S = tmp$S, chk = check.constriant, domega = tmp$omega.dev))
}



# Internal: validate status vector and constraint count
.validate_kmc_inputs <- function(d, g) {
  if (length(unique(d)) != 1) {
    if (!setequal(unique(d), c(0, 1))) stop("Status must be 0/1")
  } else {
    if (d[1] != 1) stop("Status must be 0/1")
  }
  if (sum(d) < length(g)) {
    stop("Number of observation MUST be greater than numbers of constraints")
  }
}

# Internal: parse and validate control parameters
.parse_kmc_control <- function(control, p) {
  nr.it <- if ("nr.it" %in% names(control)) max(control$nr.it, 10) else 20
  nr.c <- if ("nr.c" %in% names(control)) {
    if (control$nr.c > 1) {
      warning("In N-R iteration, C should be between 0 and 1")
      1
    } else {
      control$nr.c
    }
  } else {
    1
  }
  em.it <- if ("em.it" %in% names(control)) min(control$em.it, 10) else 3
  experimental <- "experimental" %in% names(control)
  default.init <- if ("default.init" %in% names(control)) {
    control[["default.init"]]
  } else {
    rep(0., p)
  }
  list(nr.it = nr.it, nr.c = nr.c, em.it = em.it,
       experimental = experimental, default.init = default.init)
}

kmc.solve <- function(x, d, g, em.boost = T, using.num = T, using.Fortran = T, using.C = F, tmp.tag = T, rtol = 1E-9, control = list(nr.it = 20, nr.c = 1, em.it = 3), ...) {
  .validate_kmc_inputs(d, g)
  ctrl <- .parse_kmc_control(control, length(g))
  nr.it <- ctrl$nr.it
  nr.c <- ctrl$nr.c
  em.it <- ctrl$em.it
  experimental <- ctrl$experimental
  default.init <- ctrl$default.init

  # Data preprocessing: sort, break ties, ensure proper censoring structure
  re <- kmc.clean(kmc.time = x, delta = d)
  kmc.time <- re$kmc.time
  delta <- re$delta

  p <- length(g)
  if (tmp.tag) delta[1:p] <- 1
  n <- length(delta)
  gt.mat <- matrix(0, p, n)
  for (i in 1:p) gt.mat[i, ] <- g[[i]](kmc.time)

  # Objective function for root finding (C implementation)
  kmc.comb123 <- function(x) {
    kmc_routine4(lambda = x, delta = delta, gtmat = gt.mat)
  }

  # Experimental mode (early return)
  if (experimental) {
    fun.C <- function(lam) {
      w <- kmc_routine5(delta, lam, gt.mat)
      return(gt.mat %*% w)
    }
    Cmean <- NULL
    if (p == 1) {
      fun.C2 <- function(lam) {
        w <- kmc_routine5(delta, lam, gt.mat)
        return(1 - sum(w))
      }
      Cmean <- multiroot(start = -2 + default.init, f = fun.C2)
    }
    Cmu <- multiroot(start = default.init, f = fun.C)
    return(list(Cmu, Cmean))
  }

  # Compute initial lambda via EM boost or zero initialization
  if (em.boost) {
    if (p == 1) {
      em_re <- el.cen.EM.kmc(x = kmc.time, d = delta, fun = g[[1]],
                              mu = 0, maxit = em.it, debug.kmc = F)
      init.lam <- (n - 1 / em_re$prob[1]) / g[[1]](em_re$times[1])
    } else if (p == 2) {
      em_re <- el.cen.EM2.kmc(x = kmc.time, d = delta,
                               fun = function(x) cbind(g[[1]](x), g[[2]](x)),
                               mu = c(0, 0), maxit = 5, debug.kmc = F)
      del.loc <- which(delta == 1)[1:2]
      tmp <- c(0, 0)
      if (del.loc[2] != 2) {
        tmp[2] <- sum(as.numeric(delta[1:(del.loc[2] - 1)] == 0) /
                        rep(1 - em_re$prob[1], 2))
      }
      UD <- cbind(g[[1]](em_re$times[1:2]), g[[2]](em_re$times[1:2]))
      init.lam <- as.vector(solve(UD) %*% (n - 1 / em_re$prob[del.loc] - tmp))
    } else {
      init.lam <- rep(0, p)
    }
  } else {
    init.lam <- rep(0, p)
  }

  # Root finding for lambda
  if (using.num || (p != 1)) {
    lambda <- multiroot(kmc.comb123, start = init.lam,
                        ctol = rtol, useFortran = using.Fortran)$root
  } else {
    # Custom Newton-Raphson for single-constraint analytic derivative
    kmc.comb12 <- Vectorize(function(x) {
      re <- kmc.data12(kmc.time, delta, lambda = x, g, gt.mat = gt.mat)
      list(x = re$chk, dev = re$domega)
    })
    multiroot.nr <- function(f_, xinit, it = nr.it, C = nr.c, tol = 1E-9) {
      if (C * tol > 1) C <- ceiling(1 / tol / 10)
      re <- xinit
      for (i in 1:it) {
        tmp <- f_(re)
        if (abs(tmp[[1]]) < tol) break
        re <- re - tmp[[1]] / tmp[[2]] * C
      }
      if (i == it) message("May not converge.")
      re
    }
    lambda <- multiroot.nr(f_ = kmc.comb12, xinit = init.lam, it = 15, C = 1, tol = rtol)
  }

  # Null hypothesis log-likelihood (unconstrained KM)
  if (em.boost & (p == 1)) {
    loglik.null <- WKM(kmc.time, delta)$logel
  } else {
    re0 <- omega.lambda(kmc.time = kmc.time, delta = delta,
                        lambda = 0, g = g, gt.mat = gt.mat)
    loglik.null <- kmc.el(delta, re0$omega, re0$S)
  }

  # Alternative hypothesis log-likelihood
  result <- tryCatch(
    omega.lambda(kmc.time, delta, lambda, g, gt.mat = gt.mat),
    error = function(cond) {
      message(cond)
      list(S = NA, omega = NA, gt = NA)
    }
  )

  if (!is.na(result$S[1])) {
    loglik.ha <- kmc.el(delta, result$omega, result$S)
    llr <- -2 * (loglik.ha - loglik.null)
    re.tmp <- list(
      loglik.null = loglik.null,
      loglik.h0 = loglik.ha,
      "-2LLR" = llr,
      g = g,
      time = x,
      status = d,
      phat = result$omega,
      pvalue = 1 - pchisq(llr, df = p),
      lambda = lambda
    )
    if (llr > 100) warning("\nThe results may be not feasible!\n")
  } else {
    re.tmp <- list(
      loglik.null = loglik.null,
      loglik.h0 = NA,
      "-2LLR" = NA,
      g = g,
      time = x,
      status = d,
      phat = NA,
      pvalue = NA,
      df = NA,
      lambda = NA
    )
  }
  class(re.tmp) <- "kmcS3"
  return(re.tmp)
}

#' The `kmc.solve2` function calculates the omega and lambda
#' for the kmc data.
#' @param x: time
#' @param d {0,1} indicator of observer/censored
#' @param g the constraints
#' @param rtol the tolerance
#' @param control the control parameters
#' @param ... other parameters
#' @return a list of the results
kmc.solvelite <- function(
    x,
    d,
    g,
    rtol = 1E-9,
    control = list(nr.it = 20, nr.c = 1, em.it = 3), ...) {
  .validate_kmc_inputs(d, g)
  .parse_kmc_control(control, length(g))  # validates control params

  re <- kmc.clean(kmc.time = x, delta = d)
  kmc.time <- re$kmc.time
  delta <- re$delta

  p <- length(g)
  if (sum(delta == 1) < p) {
    warning("Number of uncensored observations must be greater than the number of constraints")
  }
  delta[1:p] <- 1
  n <- length(delta)
  gt.mat <- matrix(0, p, n)
  for (i in 1:p) gt.mat[i, ] <- g[[i]](kmc.time)

  init.lam <- rep(0, p)

  # Null hypothesis (lambda = 0, unconstrained KM)
  re0 <- omega.lambda(kmc.time = kmc.time, delta = delta,
                      lambda = init.lam, g = g, gt.mat = gt.mat)
  loglik.null <- kmc.el(delta, re0$omega, re0$S)

  # Root finding with iterative refinement
  kmc.comb_inner <- function(x) {
    kmc_routine4(lambda = x, delta = delta, gtmat = gt.mat)
  }

  mini_iter <- 10
  iter <- 0
  lambda <- init.lam
  current_rtol <- rtol

  while (iter < mini_iter) {
    result_h1 <- multiroot(kmc.comb_inner, start = lambda,
                           ctol = current_rtol, maxiter = 512)
    lambda <- result_h1$root
    iter <- iter + result_h1$iter
    current_rtol <- current_rtol / 2
  }

  # Alternative hypothesis
  re1 <- omega.lambda(kmc.time = kmc.time, delta = delta,
                      lambda = lambda, g = g, gt.mat = gt.mat)
  loglik.ha <- kmc.el(delta, re1$omega, re1$S)

  convergence <- 1 - (sum(re1$omega < 0.) > 0.)
  llr <- 2. * (loglik.null - loglik.ha)

  re.tmp <- list(
    loglik.null = loglik.null,
    loglik.h0 = loglik.ha,
    "-2LLR" = llr,
    g = g,
    time = x,
    status = d,
    phat = re1$omega,
    pvalue = 1 - pchisq(llr, df = p),
    lambda = lambda,
    convergence = convergence
  )
  class(re.tmp) <- "kmcS3"
  return(re.tmp)
}


kmc.bjtest <- function(
    y, d, x, beta, init.st = "naive") {
  n <- length(y)
  x <- as.matrix(x)
  xdim <- dim(x)
  if (xdim[1] != n) {
    stop("check dim of x")
  }
  if (length(beta) != xdim[2]) {
    stop("check dim of x and beta")
  }
  e <- y - as.vector(x %*% beta)
  e_eps_ <- (1:length(e)) * 1e-12
  e <- e + e_eps_
  ordere <- order(e, -d)
  esort <- e[ordere]
  dsort <- d[ordere]
  xsort <- as.matrix(x[ordere, ])
  dsort[length(dsort)] <- 1
  temp0 <- WKM(esort, dsort, zc = 1:n)
  pKM <- temp0$jump
  temp <- redistF(y = esort, d = dsort, Fdist = pKM)
  weight <- temp$weight / n
  A <- matrix(0, ncol = xdim[2], nrow = n)
  for (i in 1:n) {
    if (dsort[i] == 1) {
      A[i, ] <- t(as.matrix(weight[1:i, i])) %*% xsort[1:i, ]
      A[i, ] <- A[i, ] / pKM[i]
    }
  }

  gt.matrix <- t(A * esort)
  delta <- dsort
  kmc.time <- esort

  kmc.comb123 <- function(x) {
    kmc_routine4(lambda = x, delta = delta, gtmat = gt.matrix) -> re
    return(re)
  }


  u.lambda2 <- function(re = el.cen.EM2.kmc(x = kmc.time, d = delta, fun = function(t, q) {
                          t * q
                        }, mu = c(0, 0), maxit = 5, debug.kmc = F, q = A)) {
    del.loc <- which(delta == 1)[1:2]
    tmp <- c(0, 0)
    if (del.loc[2] != 2) tmp[2] <- sum(as.numeric(delta[1:(del.loc[2] - 1)] == 0) / (rep(1 - re$prob[1], 2)))
    UD <- cbind(gt.matrix[1:2, 1], gt.matrix[1:2, 2])
    uu.lambda <- as.vector(
      solve(UD) %*% (n - 1 / re$prob[del.loc] - tmp)
    )
    # debug oupur lambda: print(uu.lambda)
    uu.lambda
  }

  if (init.st == "naive") {
    init.lam <- c(0, 0)
  } else {
    init.lam <- u.lambda2()
  }
  # cat("init.lam:\t",init.lam,'\tLAM')
  multiroot(kmc.comb123, start = init.lam, useFortran = T, rtol = 1e-9, atol = 1e-9, ctol = 1e-9)$root -> lambda
  # cat(lambda,'\n')
  omega.lambda(kmc.time = esort, delta = delta, lambda = lambda, g = NULL, gt.mat = gt.matrix) -> result ## set lambda=0, it compute KM-est
  temp2 <- kmc.el(delta, result$omega, result$S)
  pnew <- result$omega
  logel1 <- temp0$logel
  logel2 <- temp2
  list(prob = pnew, logel = logel1, logel2 = logel2, `-2LLR` = 2 *
    (logel1 - logel2), convergence = c(0, 1)[(abs(sum(pnew) - 1) < 0.001) + 0])
}


plotkmc2D <- function(resultkmc, flist = list(f1 = function(x) {
                        x
                      }, f2 = function(x) {
                        x^2
                      }), range0 = c(0.2, 3, 20)) {
  tmp.df <- length(resultkmc$g)
  xx <- resultkmc[["-2LLR"]]
  xl <- seq(0, max(6, xx + 2), 0.01)
  plot(xl, dchisq(xl, df = tmp.df), type = "l", main = "Kaplan-Meier Estimator with Constraint", xlab = "X", ylab = "Probabilty")
  points(xx, dchisq(xx, df = tmp.df), col = "red", lty = 2, type = "h")

  if (tmp.df == 2) {
    # X11();
    theta0 <- -do.call(c, lapply(resultkmc$g, function(x) (x(0))))
    x.grid <- seq((theta0[1] - range0[1]), (theta0[1] + range0[1]), length.out = range0[3])
    y.grid <- seq((theta0[2] - range0[2]), (theta0[2] + range0[2]), length.out = range0[3])
    tmp.z <- matrix(0, range0[3], range0[3])
    for (ii in 1:range0[3]) {
      for (jj in 1:range0[3]) {
        tmpg <- list(f1 = function(xuu) {
          flist[[1]](xuu) - tmp1
        }, f2 = function(xuu) {
          flist[[2]](xuu) - tmp2
        })
        tmp1 <- x.grid[ii]
        tmp2 <- y.grid[jj]
        tmp.z[ii, jj] <- kmc.solve(resultkmc$time, resultkmc$status, tmpg)[[2]]
      }
    }
    contour(x.grid, y.grid, tmp.z)
    points(theta0[1], theta0[2], main = "CI", col = "red")
  }

  par(mfrow = c(1, 1))
  return(list(X = x.grid, Y = y.grid, Z = tmp.z))
}


# ============================================================
# kmc.bcm.test
#   Empirical-likelihood ratio test of the regression coefficient
#   beta in the binary-choice / current-status (Case-1 interval
#   censoring) model
#
#       y_i = 1{ beta^T x_i + eps_i > 0 },  delta_i = 1 - y_i.
#
#   Identification convention: beta is identified only up to
#   positive scale.  One coordinate (default: the first) is
#   dropped from the Owen-EL constraint system; df = p - 1.
#
#   Two test variants are provided:
#     centered = FALSE  (default): raw Buckley-James estimating
#         function  g_i = x_i * hat_eps_i.  Finite-sample chi^2
#         calibration relies on a fortuitous PAVA-induced
#         cancellation; see Section 3.6 of the book.
#     centered = TRUE: centered estimating function tilde g_i =
#         (x_i - hat m(e_i)) * hat_eps_i where hat m is a
#         Nadaraya-Watson kernel regression of x on the residual
#         axis e.  Discharges the no-bias condition (A7) by the
#         tower property (Theorem 3.10).
#
#   Returns a kmcS3 list with the same shape as kmc.bjtest output.
# ============================================================

.kmc_bcm_pava_eps <- function(X, delta, beta) {
  ## Sort by e_i = -beta^T x_i (BCM convention: delta=1 means
  ## chose the BASE, i.e. eps <= -beta^T x); compute PAVA NPMLE
  ## hatF on the sorted delta and return BJ-imputed residuals.
  e   <- -as.numeric(X %*% beta)
  ord <- order(e)
  es  <- e[ord]
  ds  <- delta[ord]
  Xs  <- X[ord, , drop = FALSE]
  n   <- length(es)

  Fhat <- isoreg(ds)$yf
  ## jumps of Fhat in sorted order:  omega_j = Fhat_j - Fhat_{j-1}
  omega    <- diff(c(0, Fhat))
  cum_w_e  <- cumsum(omega * es)
  total_we <- cum_w_e[n]

  EPS <- 1e-12
  hat_eps <- numeric(n)
  for (i in seq_len(n)) {
    if (ds[i] == 1L && Fhat[i] > EPS) {
      hat_eps[i] <- cum_w_e[i] / Fhat[i]
    } else if (ds[i] == 0L && Fhat[i] < 1 - EPS) {
      hat_eps[i] <- (total_we - cum_w_e[i]) / (1 - Fhat[i])
    }
    ## else boundary; hat_eps[i] stays at 0
  }
  list(Xs = Xs, hat_eps = hat_eps, e_sorted = es, delta_sorted = ds, order = ord)
}

.kmc_bcm_nw_mean <- function(Xs, es, h = NULL) {
  ## Nadaraya-Watson kernel regression of Xs (n x p) on es (n)
  ## with Gaussian kernel.  Default bandwidth: Silverman's rule
  ## n^(-1/3) scaling (under-smoothing to discharge bias).
  n <- length(es)
  if (is.null(h)) h <- 1.5 * sd(es) * n^(-1/3)
  if (h <= 0) h <- 1e-6
  D  <- outer(es, es, "-") / h
  W  <- exp(-0.5 * D * D)
  Wn <- W / pmax(rowSums(W), 1e-12)
  Wn %*% Xs
}

kmc.bcm.test <- function(X, delta, beta, centered = FALSE, h = NULL, drop = 1L) {
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (length(delta) != n) stop("length(delta) must equal nrow(X)")
  if (length(beta)  != p) stop("length(beta) must equal ncol(X)")
  if (!all(delta %in% c(0L, 1L)) && !all(delta %in% c(0, 1)))
    stop("delta must be 0/1")
  if (drop < 1L || drop > p) stop("'drop' must be in 1..ncol(X)")

  bj <- .kmc_bcm_pava_eps(X, as.integer(delta), beta)

  if (centered) {
    mhat <- .kmc_bcm_nw_mean(bj$Xs, bj$e_sorted, h = h)
    g    <- (bj$Xs - mhat) * bj$hat_eps
  } else {
    g <- bj$Xs * bj$hat_eps
  }

  ## scale identification: drop one coordinate
  keep <- setdiff(seq_len(p), as.integer(drop))
  g <- g[, keep, drop = FALSE]

  df <- ncol(g)
  res <- tryCatch(
    emplik::el.test(g, mu = rep(0, df)),
    error = function(e) list(`-2LLR` = NA_real_, lambda = rep(NA_real_, df))
  )
  llr <- res[["-2LLR"]]
  pval <- if (is.na(llr) || llr < 0) NA_real_ else 1 - pchisq(llr, df = df)

  out <- list(
    `-2LLR` = llr,
    df = df,
    pvalue = pval,
    centered = centered,
    drop = drop,
    h = if (centered) (if (is.null(h)) 1.5 * sd(bj$e_sorted) * n^(-1/3) else h) else NA_real_,
    lambda = if (!is.null(res$lambda)) res$lambda else rep(NA_real_, df),
    n = n,
    p = p,
    beta = beta,
    convergence = if (is.na(llr) || llr < 0) 0L else 1L
  )
  class(out) <- "kmcS3"
  out
}
