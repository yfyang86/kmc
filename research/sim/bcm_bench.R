## sim/bcm_bench.R
##
## Runtime + Type-I head-to-head: kmc.bcm.test (KMC-style Owen-EL via
## emplik::el.test) versus a naive EM iteration for the BCM-EL null,
## measured on the same DGP across n in {200, 500, 1000, 2000, 5000}.
##
## The naive EM iterates the missing-data step (impute residuals
## under the current F-hat) and the Lagrange-multiplier update,
## without exploiting PAVA's closed form. This is the empirical
## stand-in for the EM-based alternative that the BCM chapter used
## to describe as "open work" prior to the kmc v0.4-4 release.

suppressPackageStartupMessages({
  library(compiler); library(emplik); library(rootSolve); library(survival)
})

setwd("/Users/yifanyang/Git/Dissertation")
source("kmc/R/kmcmaskel.R", chdir = TRUE)
source("kmc/R/kmc.R",       chdir = TRUE)
source("sim/bcm_sim.R",     chdir = TRUE)  # gen_logistic etc.

## --------------------------------------------------------------
## Naive EM iteration for the BCM-EL null at fixed beta.  At each
## iteration we (a) impute residuals under the current F-hat, (b)
## solve a one-step Owen-EL for the dual lambda via a small Newton
## step, then (c) re-estimate F-hat by the closed form analogous to
## one E-M cycle.  This is intentionally slow to match what an EM
## treatment would look like before the kmc package added a closed
## BCM solver.
## --------------------------------------------------------------
em_bcm_lr <- function(X, delta, beta, max_iter = 50L, tol = 1e-6) {
  n <- length(delta); p <- ncol(X)
  e <- -as.numeric(X %*% beta)
  ord <- order(e); es <- e[ord]; ds <- delta[ord]; Xs <- X[ord, , drop = FALSE]

  ## start at uniform F-hat over the sorted positions
  Fhat <- pmin(seq_len(n) / n - 0.5 / n + 0.5 / n, 1 - 1e-12)

  prev_stat <- Inf
  for (it in seq_len(max_iter)) {
    ## E-step: impute residuals via current F-hat
    omega <- diff(c(0, Fhat))
    cum_w_e <- cumsum(omega * es); total_we <- cum_w_e[n]
    hat_eps <- numeric(n)
    for (i in seq_len(n)) {
      if (ds[i] == 1L && Fhat[i] > 1e-12) {
        hat_eps[i] <- cum_w_e[i] / Fhat[i]
      } else if (ds[i] == 0L && Fhat[i] < 1 - 1e-12) {
        hat_eps[i] <- (total_we - cum_w_e[i]) / (1 - Fhat[i])
      }
    }
    g <- (Xs * hat_eps)[, -1L, drop = FALSE]  # drop first coord for scale

    ## M-step: solve Owen EL dual via emplik
    res <- tryCatch(emplik::el.test(g, mu = rep(0, ncol(g))),
                    error = function(e) NULL)
    if (is.null(res)) break
    stat <- res[["-2LLR"]]

    ## update F-hat: pool delta into PAVA (this is the M-step closed form)
    Fhat_new <- isoreg(ds)$yf
    Fhat_new <- pmin(pmax(Fhat_new, 1e-12), 1 - 1e-12)
    if (max(abs(Fhat_new - Fhat)) < tol && abs(stat - prev_stat) < tol) {
      Fhat <- Fhat_new; break
    }
    Fhat <- Fhat_new; prev_stat <- stat
  }
  list(stat = stat, n_iter = it)
}

## --------------------------------------------------------------
## Benchmark runner
## --------------------------------------------------------------
benchmark_cell <- function(n, reps, beta_true = c(1.0, 0.5),
                           dgp = gen_logistic) {
  set.seed(2026 + n)
  t_kmc  <- numeric(reps)
  t_em   <- numeric(reps)
  s_kmc  <- numeric(reps)
  s_em   <- numeric(reps)
  for (r in seq_len(reps)) {
    d  <- dgp(n, beta_true)
    t0 <- proc.time()
    s_kmc[r] <- kmc.bcm.test(d$X, d$delta, beta_true)[["-2LLR"]]
    t_kmc[r] <- (proc.time() - t0)[3L]
    t0 <- proc.time()
    em <- em_bcm_lr(d$X, d$delta, beta_true)
    t_em[r]  <- (proc.time() - t0)[3L]
    s_em[r]  <- em$stat
  }
  list(n = n, reps = reps,
       med_t_kmc = median(t_kmc), med_t_em = median(t_em),
       t1_kmc = mean(s_kmc > qchisq(.95, 1), na.rm = TRUE),
       t1_em  = mean(s_em  > qchisq(.95, 1), na.rm = TRUE),
       mean_kmc = mean(s_kmc, na.rm = TRUE),
       mean_em  = mean(s_em,  na.rm = TRUE),
       speedup = median(t_em) / median(t_kmc))
}

main_bench <- function() {
  cat("Wall-clock head-to-head: kmc.bcm.test  vs.  naive-EM BCM-EL\n")
  cat(strrep("=", 70), "\n", sep = "")
  cat(sprintf("%6s %6s %12s %12s %10s %10s %10s\n",
              "n", "reps", "kmc (s/rep)", "EM (s/rep)",
              "speedup", "T-I (kmc)", "T-I (EM)"))
  cells <- list()
  for (n in c(200L, 500L, 1000L, 2000L, 5000L)) {
    reps <- if (n <= 1000L) 200L else (if (n <= 2000L) 100L else 30L)
    cell <- benchmark_cell(n, reps)
    cat(sprintf("%6d %6d %12.4f %12.4f %10.1fx %10.3f %10.3f\n",
                n, reps, cell$med_t_kmc, cell$med_t_em, cell$speedup,
                cell$t1_kmc, cell$t1_em))
    cells[[as.character(n)]] <- cell
  }
  saveRDS(cells, "sim/results/bcm_bench.rds")
  invisible(cells)
}

if (sys.nframe() == 0L) main_bench()
