## sim/kmc_runtime.R
##
## fig/kmc_runtime.pdf -- log-log runtime comparison from the
## KMC chapter's tabulated benchmarks (Table 2.2 / Tabexp1 in
## the published Zhou-Yang 2015 paper).  Visualises the EM-vs-KMC
## speedup that the chapter currently only reports in numerical form.

## Data extracted from Table 2.2 / tabexp1 in KMC.tex:
##   N \in {200, 1000, 2000, 5000}
##   Three censoring-rate regimes:
##     60% censored (beta=1.5)
##     41% censored (beta=0.7)
##     17% censored (beta=0.2)
##   Two solvers: EM, KMC (numeric derivative)

ns        <- c(200, 1000, 2000, 5000)

em_60     <- c(0.175,  3.503, 13.935, 73.562)
em_41     <- c(0.064,  1.058,  4.104, 22.878)
em_17     <- c(0.014,  0.117,  0.425,  2.702)

kmc_60    <- c(0.011,  0.106,  0.349,  1.801)
kmc_41    <- c(0.010,  0.115,  0.385,  2.367)
kmc_17    <- c(0.008,  0.071,  0.240,  1.220)

pdf("fig/kmc_runtime.pdf", width = 9.5, height = 4.5)
par(mfrow = c(1, 2), mar = c(4.3, 4.6, 2.7, 1.2))

## --- Panel A: log-log time vs n ---
plot(NA, log = "xy",
     xlim = range(ns),
     ylim = range(c(em_60, em_41, em_17, kmc_60, kmc_41, kmc_17,
                     0.005, 100)),
     xlab = expression("sample size  " * N),
     ylab = "average seconds per replicate (log scale)",
     main = "(A)  EM vs KMC: time scaling")

## EM curves
lines(ns,  em_60, col = "firebrick", lty = 1, lwd = 1.6)
points(ns, em_60, col = "firebrick", pch = 19)
lines(ns,  em_41, col = "firebrick", lty = 2, lwd = 1.6)
points(ns, em_41, col = "firebrick", pch = 17)
lines(ns,  em_17, col = "firebrick", lty = 5, lwd = 1.6)
points(ns, em_17, col = "firebrick", pch = 15)

## KMC curves
lines(ns,  kmc_60, col = "black", lty = 1, lwd = 1.6)
points(ns, kmc_60, col = "black", pch = 19)
lines(ns,  kmc_41, col = "black", lty = 2, lwd = 1.6)
points(ns, kmc_41, col = "black", pch = 17)
lines(ns,  kmc_17, col = "black", lty = 5, lwd = 1.6)
points(ns, kmc_17, col = "black", pch = 15)

legend("topleft", inset = c(0, 0),
       legend = c(
         "EM, 60% censored",
         "EM, 41% censored",
         "EM, 17% censored",
         "KMC, 60% censored",
         "KMC, 41% censored",
         "KMC, 17% censored"),
       col = c(rep("firebrick", 3), rep("black", 3)),
       lty = rep(c(1, 2, 5), 2),
       pch = rep(c(19, 17, 15), 2),
       lwd = 1.6, cex = 0.78, bty = "n")

## fit empirical complexity exponents on log-log
fit_slope <- function(x, y) coef(lm(log(y) ~ log(x)))[2L]
alpha_em_60  <- fit_slope(ns, em_60)
alpha_em_41  <- fit_slope(ns, em_41)
alpha_em_17  <- fit_slope(ns, em_17)
alpha_kmc_60 <- fit_slope(ns, kmc_60)
alpha_kmc_41 <- fit_slope(ns, kmc_41)
alpha_kmc_17 <- fit_slope(ns, kmc_17)
mtext(sprintf("EM slope ~ %.2f, KMC slope ~ %.2f",
              mean(c(alpha_em_60, alpha_em_41, alpha_em_17)),
              mean(c(alpha_kmc_60, alpha_kmc_41, alpha_kmc_17))),
      side = 3, line = -1.4, cex = 0.78, adj = 0.95)

## --- Panel B: speedup vs n ---
speedup_60 <- em_60 / kmc_60
speedup_41 <- em_41 / kmc_41
speedup_17 <- em_17 / kmc_17
plot(NA, log = "xy",
     xlim = range(ns),
     ylim = c(1, max(c(speedup_60, speedup_41, speedup_17)) * 1.2),
     xlab = expression("sample size  " * N),
     ylab = "EM time / KMC time (log scale)",
     main = "(B)  Speedup factor")
abline(h = 1, col = "grey60", lty = 3)
lines(ns,  speedup_60, col = "steelblue", lty = 1, lwd = 1.6)
points(ns, speedup_60, col = "steelblue", pch = 19)
lines(ns,  speedup_41, col = "steelblue", lty = 2, lwd = 1.6)
points(ns, speedup_41, col = "steelblue", pch = 17)
lines(ns,  speedup_17, col = "steelblue", lty = 5, lwd = 1.6)
points(ns, speedup_17, col = "steelblue", pch = 15)
legend("topleft",
       legend = c("60% censored", "41% censored", "17% censored"),
       col = "steelblue", lty = c(1, 2, 5), pch = c(19, 17, 15),
       lwd = 1.6, cex = 0.9, bty = "n")

dev.off()
cat("Wrote fig/kmc_runtime.pdf\n")
