## sim/bcm_rates.R
##
## Empirical verification of the convergence rates claimed in the
## BCM-EL Wilks proof:
##
##   (i)   PAVA NPMLE rate:   ||hat F - F_0||_{2, F_0} = O_p(n^{-1/3})
##   (ii)  Score CLT:         sqrt(n) g_bar_n(hat F) is approx N(0, Sigma)
##   (iii) Bias smallness:    sqrt(n) E[g_bar_n(hat F)] -> 0   (no-bias)
##
## (iii) is the technical assumption (A7) of Theorem 6.1; we verify
## it empirically here because the BJ-residual estimating function is
## not the literal efficient score and (A7) is not a free lunch.

source("sim/bcm_sim.R", chdir = FALSE)

EPS_TRIM <- 1e-12

## ---------------------------------------------------------
## (i) L_2 rate of PAVA NPMLE at fixed beta_0
##     compute ||hat F^{(beta_0)} - F_0||_2 across n,
##     averaged over Monte-Carlo replicates.
## ---------------------------------------------------------
l2_rate <- function(n_grid, reps, beta_true, dgp_fn, F_true) {
  out <- matrix(NA_real_, length(n_grid), 2L,
                dimnames = list(NULL, c("n", "l2_norm")))
  for (k in seq_along(n_grid)) {
    n <- n_grid[k]
    norms <- numeric(reps)
    for (r in seq_len(reps)) {
      d  <- dgp_fn(n, beta_true)
      e  <- -as.numeric(d$X %*% beta_true)
      o  <- order(e)
      Fh <- isoreg(d$delta[o])$yf
      Ft <- F_true(e[o])
      norms[r] <- sqrt(mean((Fh - Ft)^2))
    }
    out[k, ] <- c(n, mean(norms))
  }
  as.data.frame(out)
}

## ---------------------------------------------------------
## (ii) Score CLT and (iii) bias under H0
##      compute the per-replicate empirical mean of g, then
##      examine its sample distribution.
## ---------------------------------------------------------
gbar_distribution <- function(n, reps, beta_true, dgp_fn, drop_first = TRUE) {
  p <- length(beta_true)
  q <- if (drop_first) p - 1L else p
  gbar <- matrix(NA_real_, reps, q)
  for (r in seq_len(reps)) {
    d  <- dgp_fn(n, beta_true)
    bj <- bj_residuals(d$X, d$delta, beta_true)
    g  <- bj$g
    if (drop_first) g <- g[, -1L, drop = FALSE]
    gbar[r, ] <- colMeans(g)
  }
  list(gbar = gbar,
       gbar_mean    = colMeans(gbar),
       gbar_sd      = apply(gbar, 2, sd),
       sqrtn_gbar_mean = sqrt(n) * colMeans(gbar))
}

## ---------------------------------------------------------
## Main
## ---------------------------------------------------------
main_rates <- function(reps_F  = 300L,
                       reps_g  = 5000L,
                       out_dir = "sim/results") {
  beta_true <- c(1.0, 0.5)
  F_true <- plogis  # logistic CDF

  cat("\n=== (i) L_2 rate of hat F at beta_0 (logistic, p=2) ===\n")
  ns_F <- c(100L, 200L, 400L, 800L, 1600L, 3200L)
  rate_F <- l2_rate(ns_F, reps_F, beta_true, gen_logistic, F_true)
  print(rate_F)

  ## fit log-log slope:
  fit <- lm(log(l2_norm) ~ log(n), data = rate_F)
  slope <- coef(fit)[[2L]]
  cat(sprintf("\n  empirical slope of log ||hat F - F_0||_2 vs log n:  %.3f\n",
              slope))
  cat(sprintf("  theoretical slope:                                  -0.333\n"))

  cat("\n=== (ii) score CLT and (iii) bias check at n=2000 ===\n")
  ns_g <- c(500L, 1000L, 2000L, 4000L)
  bias_tab <- data.frame(n = integer(0), comp = integer(0),
                         sqrtn_gbar_mean = numeric(0),
                         gbar_se = numeric(0),
                         shapiro_p = numeric(0))
  for (n in ns_g) {
    out <- gbar_distribution(n, reps_g, beta_true, gen_logistic)
    for (k in seq_along(out$gbar_mean)) {
      ## Shapiro-Wilk normality test on a subsample (5000-sample limit)
      sw_p <- suppressWarnings(
        shapiro.test(sample(out$gbar[, k], min(4500L, nrow(out$gbar))))
      )$p.value
      bias_tab <- rbind(bias_tab, data.frame(
        n = n,
        comp = k,
        sqrtn_gbar_mean = out$sqrtn_gbar_mean[k],
        gbar_se = out$gbar_sd[k] / sqrt(reps_g),
        shapiro_p = sw_p
      ))
    }
  }
  print(bias_tab, digits = 3)

  ## Save artefacts
  saveRDS(list(rate_F = rate_F, bias_tab = bias_tab,
               F_slope = slope),
          file = file.path(out_dir, "rates_results.rds"))

  ## ----- figure -----
  pdf("fig/bcm_rates.pdf", width = 8.0, height = 4.0)
  par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.0, 1.0))

  plot(log(rate_F$n), log(rate_F$l2_norm), pch = 19, cex = 1.0,
       xlab = expression(log * n),
       ylab = expression(log * "||" * hat(F) - F[0] * "||"[2]),
       main = "PAVA NPMLE L2-rate")
  abline(fit, lwd = 2, col = "steelblue")
  ref <- range(log(rate_F$n))
  lines(ref, fit$coef[[1L]] + (-1/3) * (ref - mean(ref)) +
        mean(log(rate_F$l2_norm)) - (-1/3) * (ref - mean(ref))[1],
        lty = 2, col = "grey50")
  legend("topright",
         legend = c(sprintf("fit: slope %.3f", slope),
                    "reference: slope -1/3"),
         col = c("steelblue", "grey50"), lty = c(1, 2), lwd = 2, bty = "n",
         cex = 0.85)

  ## Plot histogram of sqrt(n) * gbar at largest n, last coord
  last_n <- ns_g[length(ns_g)]
  out_last <- gbar_distribution(last_n, reps_g, beta_true, gen_logistic)
  sqrtn_g <- sqrt(last_n) * out_last$gbar[, 1L]
  hist(sqrtn_g, breaks = 50, freq = FALSE,
       main = bquote("CLT check: " * sqrt(n) ~ bar(g)[n] * ",  n=" * .(last_n)),
       xlab = expression(sqrt(n) ~ bar(g)[n]))
  xs <- seq(min(sqrtn_g), max(sqrtn_g), length.out = 200)
  lines(xs, dnorm(xs, mean = mean(sqrtn_g), sd = sd(sqrtn_g)),
        col = "firebrick", lwd = 2)
  mtext(sprintf("mean = %.3f,  sd = %.3f", mean(sqrtn_g), sd(sqrtn_g)),
        side = 3, cex = 0.8, line = 0.2, adj = 0.97)
  dev.off()
  cat("\nWrote fig/bcm_rates.pdf\n")

  invisible(list(rate_F = rate_F, bias_tab = bias_tab))
}

if (sys.nframe() == 0L) {
  set.seed(20260521)
  main_rates()
}
