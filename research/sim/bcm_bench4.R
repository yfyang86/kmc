## sim/bcm_bench4.R
##
## Rigorous-v3 BCM-EL benchmark with a third reference solver.
##
## Adds, on top of the bcm_bench3.R design:
##   * A direct damped-Newton solver on the Owen-EL dual,
##     written from scratch (independent of emplik::el.test),
##     as a sanity check that kmc.bcm.test's underlying solver
##     gives the correct -2 log ELR.
##   * Reports max |stat_kmc - stat_newton| across replicates
##     to verify the two agree to ~1e-9.
##   * Two-sided comparison: kmc.bcm.test ~ Newton (same problem,
##     different solver -> should agree) vs EM (different problem,
##     different solver -> systematic disagreement).

suppressPackageStartupMessages({
  library(compiler); library(emplik); library(rootSolve); library(survival)
})

setwd("/Users/yifanyang/Git/Dissertation")
source("kmc/R/kmcmaskel.R", chdir = TRUE)
source("kmc/R/kmc.R",       chdir = TRUE)
source("sim/bcm_sim.R",     chdir = TRUE)

## --------------------------------------------------------------
## Direct damped-Newton on the Owen-EL dual.
##
## Same setup as kmc.bcm.test: PAVA NPMLE of F at fixed beta,
## BJ-imputed residuals, estimating function g_i = X_i hat_eps_i,
## one coordinate dropped for scale identification.
##
## We solve the dual lambda directly:
##     0 = sum_{i=1}^{n} g_i / (1 + lambda^T g_i)
## by damped Newton's method.  Statistic:
##     -2 log ELR  =  2 sum log(1 + lambda^T g_i).
## --------------------------------------------------------------
newton_bcm_lr <- function(X, delta, beta, max_iter = 50L, tol = 1e-10) {
  X <- as.matrix(X)
  n <- length(delta)

  e   <- -as.numeric(X %*% beta)
  ord <- order(e); es <- e[ord]; ds <- delta[ord]
  Xs  <- X[ord, , drop = FALSE]

  Fhat   <- isoreg(ds)$yf
  Fhat   <- pmin(pmax(Fhat, 1e-12), 1 - 1e-12)
  omega  <- diff(c(0, Fhat))
  cum_we <- cumsum(omega * es); total_we <- cum_we[n]

  hat_eps <- numeric(n)
  for (i in seq_len(n)) {
    if (ds[i] == 1L && Fhat[i] > 1e-12) {
      hat_eps[i] <- cum_we[i] / Fhat[i]
    } else if (ds[i] == 0L && Fhat[i] < 1 - 1e-12) {
      hat_eps[i] <- (total_we - cum_we[i]) / (1 - Fhat[i])
    }
  }

  g  <- (Xs * hat_eps)[, -1L, drop = FALSE]
  p  <- ncol(g)

  lambda <- rep(0, p)
  converged <- FALSE
  conv_iter <- max_iter

  for (it in seq_len(max_iter)) {
    z    <- 1 + as.numeric(g %*% lambda)
    if (any(z <= 1e-12)) {
      lambda <- lambda * 0.5
      next
    }
    grad <- colSums(g / z)              # length p; we want this = 0
    ## Jacobian of F(lambda) = sum g/z w.r.t. lambda is
    ##   J = -sum g_i g_i^T / z_i^2.
    ## Newton update: lambda_new = lambda - J^{-1} F = lambda + H^{-1} grad
    ## with H = sum g_i g_i^T / z_i^2 = crossprod(g/z).
    H    <- crossprod(g / z)            # p x p, positive definite
    step_dir <- tryCatch(solve(H, grad), error = function(e) NULL)
    if (is.null(step_dir)) break

    alpha <- 1
    new_lambda <- lambda + alpha * step_dir
    while (any(1 + as.numeric(g %*% new_lambda) <= 1e-12) && alpha > 1e-10) {
      alpha      <- alpha * 0.5
      new_lambda <- lambda + alpha * step_dir
    }
    if (alpha <= 1e-10) break
    lambda <- new_lambda

    if (max(abs(grad)) < tol) {
      converged <- TRUE
      conv_iter <- it
      break
    }
    conv_iter <- it
  }

  z <- 1 + as.numeric(g %*% lambda)
  stat <- 2 * sum(log(z[z > 0]))
  list(stat = stat, lambda = lambda, n_iter = conv_iter, converged = converged)
}

## --------------------------------------------------------------
## EM iteration with diagnostics (same as bcm_bench3.R)
## --------------------------------------------------------------
em_bcm_lr_diag <- function(X, delta, beta, max_iter = 50L, tol = 1e-6) {
  n <- length(delta); p <- ncol(X)
  e   <- -as.numeric(X %*% beta)
  ord <- order(e); es <- e[ord]; ds <- delta[ord]
  Xs  <- X[ord, , drop = FALSE]
  w   <- rep(1 / n, n)
  prev_stat <- Inf
  stat <- NA_real_
  conv_iter <- max_iter
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    cum_w   <- cumsum(w)
    cum_w_d <- cumsum(w * ds)
    cs      <- cum_w_d / pmax(cum_w, 1e-14)
    Fhat    <- isoreg(cs)$yf
    Fhat    <- pmin(pmax(Fhat, 1e-12), 1 - 1e-12)

    omega   <- diff(c(0, Fhat))
    cum_we  <- cumsum(omega * es); total_we <- cum_we[n]
    hat_eps <- numeric(n)
    for (i in seq_len(n)) {
      if (ds[i] == 1L && Fhat[i] > 1e-12) {
        hat_eps[i] <- cum_we[i] / Fhat[i]
      } else if (ds[i] == 0L && Fhat[i] < 1 - 1e-12) {
        hat_eps[i] <- (total_we - cum_we[i]) / (1 - Fhat[i])
      }
    }
    g <- (Xs * hat_eps)[, -1L, drop = FALSE]
    res <- tryCatch(emplik::el.test(g, mu = rep(0, ncol(g))),
                    error = function(e) NULL)
    if (is.null(res)) { stat <- NA_real_; conv_iter <- it; break }
    stat   <- res[["-2LLR"]]
    lambda <- res$lambda
    denom <- 1 + as.numeric(g %*% lambda)
    w_new <- 1 / (n * denom)
    w_new[!is.finite(w_new) | w_new < 0] <- 0
    if (sum(w_new) > 0) w_new <- w_new / sum(w_new) else w_new <- w
    if (max(abs(w_new - w)) < tol && abs(stat - prev_stat) < tol) {
      conv_iter <- it; converged <- TRUE; w <- w_new; break
    }
    w <- w_new; prev_stat <- stat
    conv_iter <- it
  }
  list(stat = stat, n_iter = conv_iter, converged = converged)
}

gen_gamma <- function(n, beta_true) {
  p     <- length(beta_true)
  X     <- matrix(runif(n * p, -1, 1), n, p)
  eps   <- rgamma(n, shape = 2, rate = 1) - 2
  delta <- as.integer(as.numeric(X %*% beta_true) + eps <= 0)
  list(X = X, delta = delta)
}

boot_ratio_ci <- function(num, den, B = 1000L, conf = 0.95) {
  num <- num[num > 0]; den <- den[den > 0]
  if (!length(num) || !length(den)) return(c(NA, NA))
  rats <- numeric(B)
  for (b in seq_len(B)) {
    n_b <- median(sample(num, length(num), replace = TRUE))
    d_b <- median(sample(den, length(den), replace = TRUE))
    rats[b] <- d_b / n_b
  }
  alpha <- (1 - conf) / 2
  unname(quantile(rats, c(alpha, 1 - alpha)))
}

bench_cell <- function(n, dgp, beta_true, reps) {
  t_kmc <- s_kmc <- numeric(reps)
  t_new <- s_new <- numeric(reps)
  t_em  <- s_em  <- numeric(reps)
  for (r in seq_len(reps)) {
    d  <- dgp(n, beta_true)

    t0 <- proc.time()
    s_kmc[r] <- kmc.bcm.test(d$X, d$delta, beta_true,
                             centered = FALSE)[["-2LLR"]]
    t_kmc[r] <- (proc.time() - t0)[3L]

    t0 <- proc.time()
    nn <- newton_bcm_lr(d$X, d$delta, beta_true)
    t_new[r] <- (proc.time() - t0)[3L]
    s_new[r] <- nn$stat

    t0 <- proc.time()
    em <- em_bcm_lr_diag(d$X, d$delta, beta_true, max_iter = 50L)
    t_em[r]  <- (proc.time() - t0)[3L]
    s_em[r]  <- em$stat
  }
  list(t_kmc = t_kmc, s_kmc = s_kmc,
       t_new = t_new, s_new = s_new,
       t_em  = t_em,  s_em  = s_em)
}

summarise_cell <- function(b, df = 1L) {
  ok_kmc <- !is.na(b$s_kmc) & b$s_kmc >= 0
  ok_new <- !is.na(b$s_new) & b$s_new >= 0
  ok_em  <- !is.na(b$s_em)  & b$s_em  >= 0
  both_solved <- ok_kmc & ok_new
  list(
    med_kmc = median(b$t_kmc[ok_kmc]),
    med_new = median(b$t_new[ok_new]),
    med_em  = median(b$t_em[ok_em]),
    speedup_kmc_vs_em  = median(b$t_em[ok_em]) / median(b$t_kmc[ok_kmc]),
    speedup_kmc_vs_em_ci = boot_ratio_ci(b$t_kmc[ok_kmc], b$t_em[ok_em]),
    max_abs_diff_kmc_new = max(abs(b$s_kmc[both_solved] - b$s_new[both_solved])),
    max_abs_diff_kmc_em  = max(abs(b$s_kmc[ok_kmc & ok_em] - b$s_em[ok_kmc & ok_em])),
    t1_kmc = mean(b$s_kmc[ok_kmc] > qchisq(.95, df)),
    t1_new = mean(b$s_new[ok_new] > qchisq(.95, df)),
    t1_em  = mean(b$s_em[ok_em]   > qchisq(.95, df))
  )
}

main_bench4 <- function() {
  beta_true <- c(1.0, 0.5)
  dgps <- list(
    logistic = gen_logistic,
    normal   = gen_normal,
    gamma    = gen_gamma
  )
  ns <- c(200L, 500L, 1000L, 2000L, 5000L)
  reps_for <- function(n) if (n <= 500L) 800L else
                          if (n <= 1500L) 400L else
                          if (n <= 3000L) 200L else 100L

  cat("\nRigorous-v3 BCM-EL benchmark (with direct-Newton third solver)\n")
  cat(strrep("=", 100), "\n", sep = "")
  cat(sprintf("%-8s %5s %5s | %8s %8s %8s | %16s | %12s | %5s %5s %5s\n",
              "dgp", "n", "reps",
              "kmc (s)", "newton(s)", "em (s)",
              "|kmc-new| max",
              "|kmc-em| max",
              "T-I k", "T-I n", "T-I e"))
  cat(strrep("-", 100), "\n", sep = "")
  cells <- list()
  for (dgp_name in names(dgps)) {
    for (n in ns) {
      set.seed(2028L + 19L * n + nchar(dgp_name))
      reps <- reps_for(n)
      cell <- bench_cell(n, dgps[[dgp_name]], beta_true, reps)
      s    <- summarise_cell(cell, df = 1L)
      cat(sprintf("%-8s %5d %5d | %8.4f %8.4f %8.4f | %16.2e | %12.2f | %.3f %.3f %.3f\n",
                  dgp_name, n, reps,
                  s$med_kmc, s$med_new, s$med_em,
                  s$max_abs_diff_kmc_new,
                  s$max_abs_diff_kmc_em,
                  s$t1_kmc, s$t1_new, s$t1_em))
      cells[[paste(dgp_name, n, sep = "_")]] <- list(
        meta = list(dgp = dgp_name, n = n, reps = reps),
        raw = cell, summary = s
      )
    }
  }
  saveRDS(cells, "sim/results/bcm_bench4.rds")
  invisible(cells)
}

if (sys.nframe() == 0L) main_bench4()
