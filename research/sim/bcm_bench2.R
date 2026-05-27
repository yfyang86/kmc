## sim/bcm_bench2.R
##
## Rigorous head-to-head benchmark:
##
##   kmc.bcm.test (raw + centered, v0.4-4) vs naive iterative-EM BCM-EL
##
## across 3 residual distributions (logistic, normal, Gamma(2,1)-2)
## and n in {200, 500, 1000, 2000, 5000}.
##
## Reports median timing with bootstrap 95% CI on the speedup, and a
## numerical-agreement check |-2LLR_kmc + -2LLR_em| max so that the
## two procedures are demonstrably solving the same null.

suppressPackageStartupMessages({
  library(compiler); library(emplik); library(rootSolve); library(survival)
})

setwd("/Users/yifanyang/Git/Dissertation")
source("kmc/R/kmcmaskel.R", chdir = TRUE)
source("kmc/R/kmc.R",       chdir = TRUE)
source("sim/bcm_sim.R",     chdir = TRUE)

## --------------------------------------------------------------
## Naive iterative-EM BCM-EL solver
##   E-step: impute residuals under current Fhat
##   M-step: Owen-EL inner solve via el.test for lambda, then
##           update weights and re-PAVA Fhat under EL weighting.
## --------------------------------------------------------------
em_bcm_lr <- function(X, delta, beta, max_iter = 50L, tol = 1e-6) {
  n <- length(delta); p <- ncol(X)
  e <- -as.numeric(X %*% beta)
  ord <- order(e); es <- e[ord]; ds <- delta[ord]
  Xs <- X[ord, , drop = FALSE]
  w  <- rep(1 / n, n)
  prev_stat <- Inf
  stat <- NA_real_; iter <- 0L
  for (it in seq_len(max_iter)) {
    iter <- it
    ## (M-step for F): weighted-PAVA on ds with weights w
    cum_w   <- cumsum(w)
    cum_w_d <- cumsum(w * ds)
    cs      <- cum_w_d / pmax(cum_w, 1e-14)
    Fhat    <- isoreg(cs)$yf
    Fhat    <- pmin(pmax(Fhat, 1e-12), 1 - 1e-12)

    ## (E-step): impute residuals
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

    ## (M-step for lambda): Owen-EL inner
    res <- tryCatch(emplik::el.test(g, mu = rep(0, ncol(g))),
                    error = function(e) NULL)
    if (is.null(res)) { stat <- NA_real_; break }
    stat   <- res[["-2LLR"]]
    lambda <- res$lambda
    ## update weights from lambda
    denom <- 1 + as.numeric(g %*% lambda)
    w_new <- 1 / (n * denom)
    w_new[!is.finite(w_new) | w_new < 0] <- 0
    if (sum(w_new) > 0) w_new <- w_new / sum(w_new) else w_new <- w

    if (max(abs(w_new - w)) < tol &&
        abs(stat - prev_stat) < tol) {
      w <- w_new; break
    }
    w <- w_new; prev_stat <- stat
  }
  list(stat = stat, n_iter = iter)
}

## --------------------------------------------------------------
## DGPs --- reuse sim/bcm_sim.R for logistic / normal; add gamma here.
## --------------------------------------------------------------
gen_gamma <- function(n, beta_true) {
  p     <- length(beta_true)
  X     <- matrix(runif(n * p, -1, 1), n, p)
  eps   <- rgamma(n, shape = 2, rate = 1) - 2
  delta <- as.integer(as.numeric(X %*% beta_true) + eps <= 0)
  list(X = X, delta = delta)
}

## --------------------------------------------------------------
## Bootstrap CI on the speedup (ratio of medians)
## --------------------------------------------------------------
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

## --------------------------------------------------------------
## One cell
## --------------------------------------------------------------
bench_cell <- function(n, dgp, beta_true, reps) {
  t_kmc_raw <- t_kmc_cn <- t_em <- numeric(reps)
  s_kmc_raw <- s_kmc_cn <- s_em <- numeric(reps)
  for (r in seq_len(reps)) {
    d  <- dgp(n, beta_true)

    t0 <- proc.time()
    s_kmc_raw[r] <- kmc.bcm.test(d$X, d$delta, beta_true,
                                 centered = FALSE)[["-2LLR"]]
    t_kmc_raw[r] <- (proc.time() - t0)[3L]

    t0 <- proc.time()
    s_kmc_cn[r]  <- kmc.bcm.test(d$X, d$delta, beta_true,
                                 centered = TRUE)[["-2LLR"]]
    t_kmc_cn[r]  <- (proc.time() - t0)[3L]

    t0 <- proc.time()
    em <- em_bcm_lr(d$X, d$delta, beta_true)
    t_em[r]  <- (proc.time() - t0)[3L]
    s_em[r]  <- em$stat
  }
  list(
    t_kmc_raw = t_kmc_raw, t_kmc_cn = t_kmc_cn, t_em = t_em,
    s_kmc_raw = s_kmc_raw, s_kmc_cn = s_kmc_cn, s_em = s_em
  )
}

summarise_cell <- function(b, df = 1L) {
  ok_raw <- !is.na(b$s_kmc_raw) & b$s_kmc_raw >= 0
  ok_cn  <- !is.na(b$s_kmc_cn)  & b$s_kmc_cn  >= 0
  ok_em  <- !is.na(b$s_em)      & b$s_em      >= 0
  list(
    med_raw = median(b$t_kmc_raw[ok_raw]),
    med_cn  = median(b$t_kmc_cn[ok_cn]),
    med_em  = median(b$t_em[ok_em]),
    speedup_raw = median(b$t_em[ok_em]) / median(b$t_kmc_raw[ok_raw]),
    speedup_cn  = median(b$t_em[ok_em]) / median(b$t_kmc_cn[ok_cn]),
    ci_speedup_raw = boot_ratio_ci(b$t_kmc_raw[ok_raw], b$t_em[ok_em]),
    ci_speedup_cn  = boot_ratio_ci(b$t_kmc_cn[ok_cn],  b$t_em[ok_em]),
    t1_kmc_raw = mean(b$s_kmc_raw[ok_raw] > qchisq(.95, df)),
    t1_kmc_cn  = mean(b$s_kmc_cn[ok_cn]   > qchisq(.95, df)),
    t1_em      = mean(b$s_em[ok_em]       > qchisq(.95, df)),
    max_abs_diff = max(abs(b$s_kmc_raw[ok_raw & ok_em] - b$s_em[ok_raw & ok_em]))
  )
}

## --------------------------------------------------------------
## Driver
## --------------------------------------------------------------
main_bench2 <- function() {
  beta_true <- c(1.0, 0.5)
  dgps <- list(
    logistic = gen_logistic,
    normal   = gen_normal,
    gamma    = gen_gamma
  )
  ns       <- c(200L, 500L, 1000L, 2000L, 5000L)
  reps_for <- function(n) if (n <= 500L) 1000L else
                          if (n <= 1500L) 500L else
                          if (n <= 3000L) 250L else 80L

  cat("Rigorous BCM-EL benchmark\n")
  cat(strrep("=", 86), "\n", sep = "")
  cat(sprintf("%-9s %5s %5s | %8s %8s %8s | %14s | %s\n",
              "dgp", "n", "reps",
              "kmc raw", "kmc cn", "EM",
              "speedup raw [CI]", "stat-agree"))
  cat(strrep("-", 86), "\n", sep = "")

  cells <- list()
  for (dgp_name in names(dgps)) {
    dgp <- dgps[[dgp_name]]
    for (n in ns) {
      set.seed(2026L + 17L * n + nchar(dgp_name))
      reps <- reps_for(n)
      cell <- bench_cell(n, dgp, beta_true, reps)
      s    <- summarise_cell(cell, df = 1L)
      cat(sprintf("%-9s %5d %5d | %8.4f %8.4f %8.4f | %5.2fx [%.2f, %.2f] | %.2e\n",
                  dgp_name, n, reps,
                  s$med_raw, s$med_cn, s$med_em,
                  s$speedup_raw, s$ci_speedup_raw[1L], s$ci_speedup_raw[2L],
                  s$max_abs_diff))
      cells[[paste(dgp_name, n, sep = "_")]] <- list(meta = list(dgp = dgp_name, n = n),
                                                      raw = cell, summary = s)
    }
  }
  saveRDS(cells, "sim/results/bcm_bench2.rds")
  invisible(cells)
}

if (sys.nframe() == 0L) main_bench2()
