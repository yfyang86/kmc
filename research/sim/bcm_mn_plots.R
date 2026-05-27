## sim/bcm_mn_plots.R
##
## Two figures for the multinomial chapter:
##   fig/bcm_mn_qq.pdf   -- 2x2 Q-Q panel of rank-EL statistic vs chi^2 for
##                          four (residual, p, J) combinations at n=1500
##   fig/bcm_mn_power.pdf -- power curves at local alternatives, J=3 logistic

source("sim/bcm_rank.R", chdir = FALSE)

set.seed(20260522)

## ---------------------------------------------------------
## QQ panel
## ---------------------------------------------------------
qq_cells <- list(
  list(label = "Gumbel,  p=2, J=3",
       n = 1500L, betas = cbind(c(1.0, 0.5), c(-0.5, 1.0)),
       eps_dist = "gumbel", df = 3L),
  list(label = "Normal,  p=2, J=3",
       n = 1500L, betas = cbind(c(1.0, 0.5), c(-0.5, 1.0)),
       eps_dist = "normal", df = 3L),
  list(label = "Logistic, p=2, J=3",
       n = 1500L, betas = cbind(c(1.0, 0.5), c(-0.5, 1.0)),
       eps_dist = "logistic", df = 3L),
  list(label = "Gumbel,  p=3, J=3",
       n = 1500L, betas = cbind(c(1.0, 0.5, -0.3), c(-0.5, 1.0, 0.2)),
       eps_dist = "gumbel", df = 5L)
)
reps_qq <- 1500L

run_qq_cell <- function(cfg) {
  stats <- numeric(reps_qq)
  for (r in seq_len(reps_qq)) {
    d <- gen_multinomial(cfg$n, cfg$betas, eps_dist = cfg$eps_dist)
    stats[r] <- rank_el_lr(d$X, d$y, cfg$betas)$stat
  }
  s <- stats[!is.na(stats) & is.finite(stats) & stats >= 0]
  ks_p <- suppressWarnings(ks.test(s, "pchisq", df = cfg$df))$p.value
  list(stats = s, ks_p = ks_p)
}

cat("Computing four QQ cells (n=1500, 1500 reps each)...\n")
qq_results <- lapply(qq_cells, function(cfg) {
  cat("  ", cfg$label, "...\n")
  c(list(meta = cfg), run_qq_cell(cfg))
})

pdf("fig/bcm_mn_qq.pdf", width = 8.0, height = 7.5)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.4, 1.0))
for (cell in qq_results) {
  s   <- cell$stats; m <- length(s)
  pp  <- (seq_len(m) - 0.5) / m
  qx  <- qchisq(pp, df = cell$meta$df)
  qy  <- sort(s)
  lim <- range(c(qx, qy))
  plot(qx, qy, pch = ".", cex = 0.8, xlim = lim, ylim = lim,
       xlab = bquote("Quantiles of " * chi[.(cell$meta$df)]^2),
       ylab = "Empirical quantiles", main = cell$meta$label,
       cex.main = 0.95)
  abline(0, 1, lty = 2, col = "grey50")
  mtext(sprintf("KS p = %.3f", cell$ks_p),
        side = 3, adj = 0.97, line = 0.2, cex = 0.75)
}
dev.off()
cat("Wrote fig/bcm_mn_qq.pdf\n")

## ---------------------------------------------------------
## Power curves (J=3 logistic, p=2, vary beta_2 component 1)
## ---------------------------------------------------------
power_cells <- list(
  list(label = "n=500",  n = 500L,  betas = cbind(c(1.0, 0.5), c(-0.5, 1.0)),
       eps_dist = "gumbel", df = 3L),
  list(label = "n=1000", n = 1000L, betas = cbind(c(1.0, 0.5), c(-0.5, 1.0)),
       eps_dist = "gumbel", df = 3L),
  list(label = "n=2000", n = 2000L, betas = cbind(c(1.0, 0.5), c(-0.5, 1.0)),
       eps_dist = "gumbel", df = 3L)
)
h_grid <- c(-0.5, -0.3, -0.15, -0.05, 0, 0.05, 0.15, 0.3, 0.5)
reps_pow <- 800L

run_power <- function(cfg) {
  crit <- qchisq(0.95, df = cfg$df)
  power <- numeric(length(h_grid))
  for (k in seq_along(h_grid)) {
    hh <- h_grid[k]
    betas_test <- cfg$betas
    betas_test[1L, 2L] <- betas_test[1L, 2L] + hh
    rej <- 0L
    for (r in seq_len(reps_pow)) {
      d   <- gen_multinomial(cfg$n, cfg$betas, eps_dist = cfg$eps_dist)
      out <- rank_el_lr(d$X, d$y, betas_test)
      if (!is.na(out$stat) && out$stat > crit) rej <- rej + 1L
    }
    power[k] <- rej / reps_pow
  }
  data.frame(h = h_grid, power = power)
}

cat("\nComputing power curves at 3 sample sizes (800 reps per h)...\n")
power_results <- lapply(power_cells, function(cfg) {
  cat("  ", cfg$label, "...\n")
  list(meta = cfg, power = run_power(cfg))
})

pdf("fig/bcm_mn_power.pdf", width = 7.0, height = 5.5)
par(mar = c(4.2, 4.4, 0.6, 1.0))
plot(NA, xlim = c(-0.55, 0.55), ylim = c(0, 1),
     xlab = expression("local alternative  " * h *
                       "  (test " * H[0] * ": " * beta[2*","*1] * " = " *
                       beta[2*","*1]^0 + h),
     ylab = "rejection probability at level 0.05")
abline(h = 0.05, lty = 3, col = "grey60")
abline(v = 0,    lty = 3, col = "grey60")
cols <- c("grey60", "black", "steelblue")
ltys <- c(2, 1, 1)
labs <- sapply(power_results, function(x) x$meta$label)
for (k in seq_along(power_results)) {
  pw <- power_results[[k]]$power
  lines(pw$h, pw$power, col = cols[k], lwd = 2, lty = ltys[k])
  points(pw$h, pw$power, col = cols[k], pch = 19, cex = 0.7)
}
legend("topleft", legend = labs, col = cols, lty = ltys, lwd = 2,
       bty = "n", cex = 0.95)
dev.off()
cat("Wrote fig/bcm_mn_power.pdf\n")

saveRDS(list(qq = qq_results, power = power_results),
        file = "sim/results/bcm_mn_figs.rds")
