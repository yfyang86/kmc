## sim/bcm_plots.R
##
## Produce QQ plot and power-curve plot from saved RDS files for the
## BCM-EL chapter.

cells <- readRDS("sim/results/sim_results.rds")
pow   <- readRDS("sim/results/power_results.rds")

dir.create("fig", showWarnings = FALSE)

## ---------------------------------------------------------
## 1.  QQ panels: stats vs chi^2_{df} for the n=1000 cells
##     (logistic & normal at df=1; logistic & gamma at df=2)
## ---------------------------------------------------------
qq_one <- function(s, df, ttl, ks_p) {
  s <- s[!is.na(s) & is.finite(s) & s >= 0]
  m <- length(s)
  pp  <- (seq_len(m) - 0.5) / m
  qx  <- qchisq(pp, df)
  qy  <- sort(s)
  lim <- range(c(qx, qy))
  plot(qx, qy, pch = ".", cex = 0.7, xlim = lim, ylim = lim,
       xlab = bquote("Quantiles of " * chi[.(df)]^2),
       ylab = "Empirical quantiles",
       main = ttl, cex.main = 0.95)
  abline(0, 1, lty = 2, col = "grey50")
  mtext(sprintf("KS p = %.3f", ks_p), side = 3, cex = 0.75, line = 0.2,
        adj = 0.97)
}

pdf("fig/bcm_qq.pdf", width = 8.0, height = 7.5)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.2, 1.0))
for (key in c("logistic_p2_n1000", "normal_p2_n1000",
              "logistic_p3_n1000", "gamma_p3_n1000")) {
  ce <- cells[[key]]
  ttl <- switch(key,
                logistic_p2_n1000 = "Logistic residuals, p=2, n=1000",
                normal_p2_n1000   = "Normal residuals,   p=2, n=1000",
                logistic_p3_n1000 = "Logistic residuals, p=3, n=1000",
                gamma_p3_n1000    = "Gamma residuals,    p=3, n=1000")
  qq_one(ce$raw$stats, ce$summary$df, ttl, ce$summary$ks_p)
}
dev.off()
cat("Wrote fig/bcm_qq.pdf\n")

## ---------------------------------------------------------
## 2.  Power curves (across n) at df=1, logistic
## ---------------------------------------------------------
pdf("fig/bcm_power.pdf", width = 7.0, height = 5.5)
par(mar = c(4.2, 4.4, 0.8, 1.0))
plot(NA, xlim = c(-0.55, 1.05), ylim = c(0, 1),
     xlab = expression("Local alternative " * h * "  (test H"[0] * ":  " *
                       beta[2] * " = " * beta[2]^0 * " + h)"),
     ylab = "Rejection probability at level 0.05")
abline(h = 0.05, lty = 3, col = "grey60")
abline(v = 0,    lty = 3, col = "grey60")

cols <- c("logistic_p2_n500"  = "grey60",
          "logistic_p2_n1000" = "black",
          "logistic_p2_n2000" = "steelblue",
          "logistic_p3_n1000" = "firebrick")
ltys <- c("logistic_p2_n500"  = 2L,
          "logistic_p2_n1000" = 1L,
          "logistic_p2_n2000" = 1L,
          "logistic_p3_n1000" = 5L)
labs <- c("logistic_p2_n500"  = "p=2, n=500   (df=1)",
          "logistic_p2_n1000" = "p=2, n=1000  (df=1)",
          "logistic_p2_n2000" = "p=2, n=2000  (df=1)",
          "logistic_p3_n1000" = "p=3, n=1000  (df=2)")

for (key in names(cols)) {
  pw <- pow[[key]]$power
  lines(pw$h, pw$power, col = cols[[key]], lwd = 2, lty = ltys[[key]])
  points(pw$h, pw$power, col = cols[[key]], pch = 19, cex = 0.7)
}
legend("topleft", legend = labs[names(cols)],
       col = cols, lty = ltys, lwd = 2, bty = "n", cex = 0.9)
dev.off()
cat("Wrote fig/bcm_power.pdf\n")
