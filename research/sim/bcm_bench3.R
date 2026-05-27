## sim/bcm_bench3.R
##
## Rigorous-v2 BCM-EL head-to-head benchmark.
##
## Improvements over bcm_bench2.R:
##   1. More replicates at large n (n=5000 -> 200 reps).
##   2. Track EM iteration counts (median, max, % hitting max_iter).
##   3. Two EM variants: short (max_iter=50) and long (max_iter=500),
##      to rule out under-convergence as the cause of miscalibration.
##   4. KS p-values for the full distribution of -2LLR against chi^2_1,
##      not just Type-I at 0.05.
##   5. Bootstrap 95% CIs on Type-I rate and on speedup.
##   6. Log-log timing regression to estimate the empirical complexity
##      exponent: time ~ n^alpha for kmc and EM separately.
##   7. Cross-tabulate kmc.bcm.test vs EM statistics on the same
##      Monte-Carlo replicates (paired comparison, not marginal).

suppressPackageStartupMessages({
  library(compiler); library(emplik); library(rootSolve); library(survival)
})

setwd("/Users/yifanyang/Git/Dissertation")
source("kmc/R/kmcmaskel.R", chdir = TRUE)
source("kmc/R/kmc.R",       chdir = TRUE)
source("sim/bcm_sim.R",     chdir = TRUE)

## --------------------------------------------------------------
## EM iteration with diagnostics returned
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

    if (max(abs(w_new - w)) < tol &&
        abs(stat - prev_stat) < tol) {
      conv_iter <- it
      converged <- TRUE
      w <- w_new
      break
    }
    w <- w_new; prev_stat <- stat
    conv_iter <- it
  }
  list(stat = stat, n_iter = conv_iter, converged = converged,
       hit_max = (conv_iter >= max_iter) && !converged)
}

## --------------------------------------------------------------
## Bootstrap helpers
## --------------------------------------------------------------
boot_ci <- function(x, stat_fn, B = 1000L, conf = 0.95) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA, NA))
  vals <- replicate(B, stat_fn(sample(x, length(x), replace = TRUE)))
  alpha <- (1 - conf) / 2
  unname(quantile(vals, c(alpha, 1 - alpha)))
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

gen_gamma <- function(n, beta_true) {
  p     <- length(beta_true)
  X     <- matrix(runif(n * p, -1, 1), n, p)
  eps   <- rgamma(n, shape = 2, rate = 1) - 2
  delta <- as.integer(as.numeric(X %*% beta_true) + eps <= 0)
  list(X = X, delta = delta)
}

## --------------------------------------------------------------
## One cell (with EM tracked twice: short and long max_iter)
## --------------------------------------------------------------
bench_cell <- function(n, dgp, beta_true, reps) {
  t_kmc <- s_kmc <- numeric(reps)
  t_em_short <- s_em_short <- iters_em_short <- numeric(reps)
  hit_max_short <- logical(reps); converged_short <- logical(reps)
  t_em_long <- s_em_long <- iters_em_long <- numeric(reps)
  hit_max_long <- logical(reps); converged_long <- logical(reps)

  for (r in seq_len(reps)) {
    d  <- dgp(n, beta_true)

    t0 <- proc.time()
    s_kmc[r] <- kmc.bcm.test(d$X, d$delta, beta_true,
                             centered = FALSE)[["-2LLR"]]
    t_kmc[r] <- (proc.time() - t0)[3L]

    t0 <- proc.time()
    em_s <- em_bcm_lr_diag(d$X, d$delta, beta_true, max_iter = 50L)
    t_em_short[r]      <- (proc.time() - t0)[3L]
    s_em_short[r]      <- em_s$stat
    iters_em_short[r]  <- em_s$n_iter
    hit_max_short[r]   <- em_s$hit_max
    converged_short[r] <- em_s$converged

    t0 <- proc.time()
    em_l <- em_bcm_lr_diag(d$X, d$delta, beta_true, max_iter = 500L)
    t_em_long[r]      <- (proc.time() - t0)[3L]
    s_em_long[r]      <- em_l$stat
    iters_em_long[r]  <- em_l$n_iter
    hit_max_long[r]   <- em_l$hit_max
    converged_long[r] <- em_l$converged
  }
  list(
    t_kmc = t_kmc, s_kmc = s_kmc,
    t_em_short = t_em_short, s_em_short = s_em_short,
    iters_em_short = iters_em_short, hit_max_short = hit_max_short,
    converged_short = converged_short,
    t_em_long = t_em_long, s_em_long = s_em_long,
    iters_em_long = iters_em_long, hit_max_long = hit_max_long,
    converged_long = converged_long
  )
}

summarise_cell <- function(b, df = 1L) {
  ok_kmc <- !is.na(b$s_kmc) & b$s_kmc >= 0
  ok_em  <- !is.na(b$s_em_short) & b$s_em_short >= 0
  ok_eml <- !is.na(b$s_em_long)  & b$s_em_long  >= 0
  q05 <- qchisq(.95, df)

  t1_kmc <- mean(b$s_kmc[ok_kmc] > q05)
  t1_em  <- mean(b$s_em_short[ok_em] > q05)
  t1_eml <- mean(b$s_em_long[ok_eml] > q05)
  t1_kmc_ci <- boot_ci(b$s_kmc[ok_kmc],     function(x) mean(x > q05))
  t1_em_ci  <- boot_ci(b$s_em_short[ok_em], function(x) mean(x > q05))
  t1_eml_ci <- boot_ci(b$s_em_long[ok_eml], function(x) mean(x > q05))

  ks_kmc <- suppressWarnings(ks.test(b$s_kmc[ok_kmc], "pchisq", df = df))$p.value
  ks_em  <- suppressWarnings(ks.test(b$s_em_short[ok_em], "pchisq", df = df))$p.value
  ks_eml <- suppressWarnings(ks.test(b$s_em_long[ok_eml], "pchisq", df = df))$p.value

  list(
    med_kmc = median(b$t_kmc[ok_kmc]),
    med_em  = median(b$t_em_short[ok_em]),
    med_eml = median(b$t_em_long[ok_eml]),
    speedup = median(b$t_em_short[ok_em]) / median(b$t_kmc[ok_kmc]),
    speedup_ci = boot_ratio_ci(b$t_kmc[ok_kmc], b$t_em_short[ok_em]),
    t1_kmc = t1_kmc, t1_em = t1_em, t1_eml = t1_eml,
    t1_kmc_ci = t1_kmc_ci, t1_em_ci = t1_em_ci, t1_eml_ci = t1_eml_ci,
    ks_kmc = ks_kmc, ks_em = ks_em, ks_eml = ks_eml,
    em_iters_med = median(b$iters_em_short[ok_em]),
    em_iters_max = max(b$iters_em_short[ok_em]),
    em_hit_max_frac = mean(b$hit_max_short),
    eml_iters_med = median(b$iters_em_long[ok_eml]),
    eml_iters_max = max(b$iters_em_long[ok_eml]),
    eml_hit_max_frac = mean(b$hit_max_long)
  )
}

## --------------------------------------------------------------
## Driver with empirical complexity exponent fit
## --------------------------------------------------------------
main_bench3 <- function() {
  beta_true <- c(1.0, 0.5)
  dgps <- list(
    logistic = gen_logistic,
    normal   = gen_normal,
    gamma    = gen_gamma
  )
  ns <- c(200L, 500L, 1000L, 2000L, 5000L)
  reps_for <- function(n) if (n <= 500L) 1500L else
                          if (n <= 1500L) 800L else
                          if (n <= 3000L) 400L else 200L

  cat("\nRigorous-v2 BCM-EL benchmark\n")
  cat(strrep("=", 110), "\n", sep = "")
  cat(sprintf("%-8s %5s %5s | %8s %8s %8s | %7s | %5s %5s | %5s %5s | %5s %5s\n",
              "dgp", "n", "reps",
              "kmc (s)", "em50 (s)", "em500 (s)",
              "speedup",
              "T-I", "KS-p", "T-I", "KS-p", "T-I", "KS-p"))
  cat(sprintf("%-8s %5s %5s | %8s %8s %8s | %7s | %11s | %11s | %11s\n",
              "", "", "", "", "", "", "",
              "kmc", "em50", "em500"))
  cat(strrep("-", 110), "\n", sep = "")

  cells <- list()
  for (dgp_name in names(dgps)) {
    dgp <- dgps[[dgp_name]]
    for (n in ns) {
      set.seed(2027L + 23L * n + nchar(dgp_name))
      reps <- reps_for(n)
      cell <- bench_cell(n, dgp, beta_true, reps)
      s    <- summarise_cell(cell, df = 1L)
      cat(sprintf("%-8s %5d %5d | %8.4f %8.4f %8.4f | %5.1fx | %.3f %.3f | %.3f %.3f | %.3f %.3f\n",
                  dgp_name, n, reps,
                  s$med_kmc, s$med_em, s$med_eml, s$speedup,
                  s$t1_kmc, s$ks_kmc,
                  s$t1_em,  s$ks_em,
                  s$t1_eml, s$ks_eml))
      cells[[paste(dgp_name, n, sep = "_")]] <- list(
        meta = list(dgp = dgp_name, n = n, reps = reps),
        raw = cell, summary = s
      )
    }
  }

  ## empirical complexity exponent: time ~ n^alpha
  cat("\nEmpirical complexity (median time ~ n^alpha):\n")
  for (dgp_name in names(dgps)) {
    times_kmc <- sapply(ns, function(n) cells[[paste(dgp_name, n, sep = "_")]]$summary$med_kmc)
    times_em  <- sapply(ns, function(n) cells[[paste(dgp_name, n, sep = "_")]]$summary$med_em)
    times_em_l<- sapply(ns, function(n) cells[[paste(dgp_name, n, sep = "_")]]$summary$med_eml)
    ## guard against zero medians at small n (timing resolution): drop any 0
    valid_kmc <- times_kmc > 0
    valid_em  <- times_em  > 0
    valid_eml <- times_em_l > 0
    alpha_kmc <- if (sum(valid_kmc) >= 3) coef(lm(log(times_kmc[valid_kmc]) ~ log(ns[valid_kmc])))[2L] else NA
    alpha_em  <- if (sum(valid_em)  >= 3) coef(lm(log(times_em[valid_em])   ~ log(ns[valid_em])))[2L] else NA
    alpha_eml <- if (sum(valid_eml) >= 3) coef(lm(log(times_em_l[valid_eml])~ log(ns[valid_eml])))[2L] else NA
    cat(sprintf("  %-8s: kmc alpha = %.2f, em50 alpha = %.2f, em500 alpha = %.2f\n",
                dgp_name, alpha_kmc, alpha_em, alpha_eml))
  }

  ## convergence diagnostics
  cat("\nEM-50 convergence diagnostics:\n")
  cat(sprintf("  %-9s %5s | %14s | %18s\n", "dgp", "n", "median iters", "frac hit max_iter"))
  for (dgp_name in names(dgps)) {
    for (n in ns) {
      s <- cells[[paste(dgp_name, n, sep = "_")]]$summary
      cat(sprintf("  %-9s %5d | %14.1f | %18.3f\n",
                  dgp_name, n, s$em_iters_med, s$em_hit_max_frac))
    }
  }

  cat("\nEM-500 convergence diagnostics:\n")
  cat(sprintf("  %-9s %5s | %14s | %18s\n", "dgp", "n", "median iters", "frac hit max_iter"))
  for (dgp_name in names(dgps)) {
    for (n in ns) {
      s <- cells[[paste(dgp_name, n, sep = "_")]]$summary
      cat(sprintf("  %-9s %5d | %14.1f | %18.3f\n",
                  dgp_name, n, s$eml_iters_med, s$eml_hit_max_frac))
    }
  }

  saveRDS(cells, "sim/results/bcm_bench3.rds")
  invisible(cells)
}

if (sys.nframe() == 0L) main_bench3()
