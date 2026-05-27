## sim/bcm_sim.R
##
## Simulation study for the empirical-likelihood test of beta in the
## binary-choice / current-status model:
##
##   y_i = 1{ beta^T x_i + eps_i > 0 },   delta_i = 1 - y_i = I{eps_i <= -beta^T x_i}.
##
## Identification: beta is identified only up to positive scale; the
## intercept is absorbed into the residual location and is not identifiable
## when F_eps is unspecified. Following Manski/Cosslett/Klein-Spady we drop
## the intercept and fix beta_1 = 1.
##
## TEST.  We use the Owen empirical likelihood with Buckley-James-style
## estimating equations (Qin & Lawless 1994, Owen 2001, applied to the
## current-status problem as described in this dissertation chapter):
##
##   g_i(beta_0) = X_i * hat_eps_i,
##
## where  hat_eps_i = E[ eps | delta_i, X_i ]  is computed by plugging the
## constrained PAVA NPMLE  hat F^{(beta_0)}  into the Buckley-James formulae:
##
##   delta_i = 1 :  hat_eps_i = sum_{j: e_j <= e_i} omega_j e_j(beta_0) / hat F(e_i)
##   delta_i = 0 :  hat_eps_i = sum_{j: e_j >  e_i} omega_j e_j(beta_0) / (1 - hat F(e_i))
##
## Under H0: beta = beta_0,  E[ g_i ] = E[ X_i eps_i ] = 0 since eps is
## independent of X with zero mean.  The EL ratio statistic
##
##   gamma_n(beta_0) = 2 sum log(1 + lambda^T g_i)
##
## with  lambda  the dual variable, is asymptotically chi^2_p by
## Qin & Lawless (1994), where p = dim(X) - 1 is the number of free
## slope coefficients after scale identification.
##
## This is the actual test that the chapter's `el.test.wt2` solver
## implements, and the one we calibrate here.

suppressPackageStartupMessages({
  library(emplik)
})

EPS_TRIM <- 1e-12

## ---------------------------------------------------------
## 1. PAVA NPMLE of F at fixed beta, returned at sorted positions
## ---------------------------------------------------------
pava_F <- function(e, delta) {
  o <- order(e)
  list(F_sorted = isoreg(delta[o])$yf,
       e_sorted = e[o],
       d_sorted = delta[o],
       order    = o)
}

## ---------------------------------------------------------
## 2. Buckley-James imputed residuals at beta_0
## ---------------------------------------------------------
bj_residuals <- function(X, delta, beta0) {
  e   <- -as.numeric(X %*% beta0)
  fit <- pava_F(e, delta)
  Fh  <- fit$F_sorted
  es  <- fit$e_sorted
  ds  <- fit$d_sorted
  ord <- fit$order
  n   <- length(es)
  Xs  <- X[ord, , drop = FALSE]

  # Jumps of F at sorted positions.  omega_j = F_j - F_{j-1}, with F_0 := 0.
  omega <- diff(c(0, Fh))
  cum_w_e <- cumsum(omega * es)             # sum_{j <= i} omega_j e_j
  total_w_e <- cum_w_e[n]

  hat_eps <- numeric(n)
  for (i in seq_len(n)) {
    if (ds[i] == 1L && Fh[i] > EPS_TRIM) {
      hat_eps[i] <- cum_w_e[i] / Fh[i]
    } else if (ds[i] == 0L && Fh[i] < 1 - EPS_TRIM) {
      hat_eps[i] <- (total_w_e - cum_w_e[i]) / (1 - Fh[i])
    } # boundary case left at 0
  }

  list(g = Xs * hat_eps, Xs = Xs, hat_eps = hat_eps, ord = ord)
}

## ---------------------------------------------------------
## 3. Owen EL test statistic for H0: beta = beta_0
##    After scale identification, only p-1 components of g are linearly
##    independent.  We test the p-1 free components.
## ---------------------------------------------------------
bcm_lr_owen <- function(X, delta, beta0) {
  bj <- bj_residuals(X, delta, beta0)
  g  <- bj$g
  # Drop the first component (the scale-fixing identification removes one
  # degree of freedom in the system X * eps = 0).  Equivalently the linear
  # system has rank p-1 generically.
  if (ncol(g) >= 2L) g <- g[, -1L, drop = FALSE]
  res <- tryCatch(
    emplik::el.test(g, mu = rep(0, ncol(g))),
    error = function(e) list(`-2LLR` = NA_real_, conv = 1L)
  )
  list(stat = res[["-2LLR"]], df = ncol(g))
}

## ---------------------------------------------------------
## 4. Data-generating processes
## ---------------------------------------------------------
gen_logistic <- function(n, beta_true) {
  p     <- length(beta_true)
  X     <- matrix(runif(n * p, -1, 1), n, p)
  eps   <- rlogis(n)
  delta <- as.integer(as.numeric(X %*% beta_true) + eps <= 0)
  list(X = X, delta = delta)
}

gen_normal <- function(n, beta_true) {
  p     <- length(beta_true)
  X     <- matrix(runif(n * p, -1, 1), n, p)
  eps   <- rnorm(n)
  delta <- as.integer(as.numeric(X %*% beta_true) + eps <= 0)
  list(X = X, delta = delta)
}

gen_gamma_centered <- function(n, beta_true) {
  p     <- length(beta_true)
  X     <- matrix(runif(n * p, -1, 1), n, p)
  eps   <- rgamma(n, shape = 2, rate = 1) - 2   # mean-zero, skewed
  delta <- as.integer(as.numeric(X %*% beta_true) + eps <= 0)
  list(X = X, delta = delta)
}

## ---------------------------------------------------------
## 5. One Monte-Carlo cell
## ---------------------------------------------------------
run_cell <- function(n, beta_true, dgp, reps, seed) {
  set.seed(seed)
  stats <- numeric(reps)
  times <- numeric(reps)
  df_used <- NA_integer_
  for (r in seq_len(reps)) {
    d  <- dgp(n, beta_true)
    t0 <- proc.time()
    out <- bcm_lr_owen(d$X, d$delta, beta_true)
    times[r] <- (proc.time() - t0)[3L]
    stats[r] <- out$stat
    df_used  <- out$df
  }
  list(stats = stats, times = times, df = df_used)
}

## ---------------------------------------------------------
## 6. Summarise one cell
## ---------------------------------------------------------
summarise_cell <- function(cell) {
  s    <- cell$stats
  df   <- cell$df
  good <- !is.na(s) & is.finite(s) & s >= 0
  s    <- s[good]
  list(
    df          = df,
    n_ok        = sum(good),
    censor_rate = NA_real_,
    type1_10    = mean(s > qchisq(0.90, df)),
    type1_05    = mean(s > qchisq(0.95, df)),
    type1_01    = mean(s > qchisq(0.99, df)),
    mean_stat   = mean(s),
    median_stat = median(s),
    sd_stat     = sd(s),
    ks_p        = suppressWarnings(ks.test(s, "pchisq", df = df))$p.value,
    median_time = median(cell$times[good])
  )
}

## ---------------------------------------------------------
## 7. Main driver
## ---------------------------------------------------------
main <- function(reps     = 1000L,
                 out_path = "sim/results/sim_results.rds") {
  cat(sprintf("\nBCM Owen-EL simulation: reps=%d per cell\n", reps))
  cat(strrep("=", 70), "\n", sep = "")

  beta_p2 <- c(1.0,  0.5)            # df = 1
  beta_p3 <- c(1.0,  0.5, -0.3)      # df = 2
  beta_p4 <- c(1.0,  0.5, -0.3, 0.4) # df = 3

  cfgs <- list()
  ns   <- c(200L, 500L, 1000L, 2000L)

  for (n in ns) {
    cfgs[[length(cfgs) + 1L]] <- list(
      name = sprintf("logistic_p2_n%d", n),
      beta = beta_p2, dgp = gen_logistic, n = n
    )
    cfgs[[length(cfgs) + 1L]] <- list(
      name = sprintf("normal_p2_n%d", n),
      beta = beta_p2, dgp = gen_normal, n = n
    )
    cfgs[[length(cfgs) + 1L]] <- list(
      name = sprintf("logistic_p3_n%d", n),
      beta = beta_p3, dgp = gen_logistic, n = n
    )
    cfgs[[length(cfgs) + 1L]] <- list(
      name = sprintf("gamma_p3_n%d", n),
      beta = beta_p3, dgp = gen_gamma_centered, n = n
    )
  }

  cells <- list()
  for (cfg in cfgs) {
    cat(sprintf("  cell %-22s ... ", cfg$name)); flush.console()
    t0 <- proc.time()
    cell  <- run_cell(cfg$n, cfg$beta, cfg$dgp, reps,
                      seed = 2026L + 17L * cfg$n + nchar(cfg$name))
    el    <- (proc.time() - t0)[3L]
    summ  <- summarise_cell(cell)
    cat(sprintf("%5.1fs | df=%d t1(.05)=%.3f mean=%.2f med=%.2f KS p=%.3f\n",
                el, summ$df, summ$type1_05, summ$mean_stat, summ$median_stat,
                summ$ks_p))
    cells[[cfg$name]] <- list(meta = cfg, raw = cell, summary = summ)
  }

  saveRDS(cells, file = out_path)
  cat("\nSaved to ", out_path, "\n", sep = "")
  invisible(cells)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  reps <- if (length(args) >= 1L) as.integer(args[1L]) else 1000L
  outp <- if (length(args) >= 2L) args[2L] else "sim/results/sim_results.rds"
  main(reps = reps, out_path = outp)
}
