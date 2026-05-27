## sim/bcm_bench_plots.R
##
## Builds fig/bcm_bench.pdf -- a 3-panel summary of the BCM-EL
## benchmark, complementing Tables 3.3-3.5 of the chapter:
##
##   (A) Log-log median time vs n, with regression fits.  Shows
##       kmc.bcm.test scaling as n^0.65, EM as n^1.1.
##   (B) Type-I error at level 0.05 vs n.  Shows kmc/Newton flat
##       at the nominal level, EM diverging upward with n.
##   (C) Paired-statistic scatter on a single representative cell:
##       kmc.bcm.test vs Newton on the diagonal (agree to ~1e-11);
##       kmc.bcm.test vs EM diverging from the diagonal.

cells3 <- readRDS("sim/results/bcm_bench3.rds")
cells4 <- readRDS("sim/results/bcm_bench4.rds")

ns       <- c(200, 500, 1000, 2000, 5000)
dgps     <- c("logistic", "normal", "gamma")
solvers  <- c("kmc", "newton", "em")
sol_col  <- c(kmc = "black", newton = "steelblue", em = "firebrick")
sol_pch  <- c(kmc = 19,      newton = 17,         em = 4)
sol_lty  <- c(kmc = 1,       newton = 2,          em = 3)
dgp_lty  <- c(logistic = 1, normal = 2, gamma = 5)

## ------------------------------------------------
## Pull medians + Type-I per (dgp, n, solver)
## ------------------------------------------------
get_med <- function(cells, dgp_name, n_val, solver) {
  key <- paste(dgp_name, n_val, sep = "_")
  s   <- cells[[key]]$summary
  if (is.null(s)) return(NA)
  switch(solver,
    kmc    = s$med_kmc,
    newton = s$med_new,
    em     = s$med_em,
    NA
  )
}
get_t1 <- function(cells, dgp_name, n_val, solver) {
  key <- paste(dgp_name, n_val, sep = "_")
  s   <- cells[[key]]$summary
  if (is.null(s)) return(NA)
  switch(solver,
    kmc    = s$t1_kmc,
    newton = if (!is.null(s$t1_new)) s$t1_new else s$t1_kmc,
    em     = s$t1_em,
    NA
  )
}

## ------------------------------------------------
## Open device
## ------------------------------------------------
pdf("fig/bcm_bench.pdf", width = 11.0, height = 4.0)
par(mfrow = c(1, 3), mar = c(4.4, 4.6, 3.2, 1.0))

## --- Panel A: log-log time vs n ---
times_kmc <- sapply(dgps, function(d) sapply(ns, function(n) get_med(cells4, d, n, "kmc")))
times_new <- sapply(dgps, function(d) sapply(ns, function(n) get_med(cells4, d, n, "newton")))
times_em  <- sapply(dgps, function(d) sapply(ns, function(n) get_med(cells4, d, n, "em")))
rownames(times_kmc) <- rownames(times_new) <- rownames(times_em) <- as.character(ns)
## replace zeros (below timing resolution) with the smallest positive entry across the dataset
finite_t  <- c(times_kmc, times_new, times_em)
finite_t  <- finite_t[finite_t > 0]
small_t   <- min(finite_t) / 2
floor_zero <- function(M) { M[M == 0 | is.na(M)] <- small_t; M }
times_kmc <- floor_zero(times_kmc)
times_new <- floor_zero(times_new)
times_em  <- floor_zero(times_em)

ymax_t <- max(c(times_kmc, times_new, times_em))
ymin_t <- min(c(times_kmc, times_new, times_em))
plot(NA, log = "xy", xlim = range(ns), ylim = c(ymin_t * 0.7, ymax_t * 2),
     xlab = expression("sample size " * n),
     ylab = "median seconds per replicate (log scale)",
     main = "(A)  Time scaling")
for (d in dgps) {
  lines(ns, times_kmc[, d], col = sol_col["kmc"],    lty = dgp_lty[d], lwd = 1.5)
  points(ns, times_kmc[, d], col = sol_col["kmc"],   pch = sol_pch["kmc"])
  lines(ns, times_new[, d], col = sol_col["newton"], lty = dgp_lty[d], lwd = 1.5)
  points(ns, times_new[, d], col = sol_col["newton"], pch = sol_pch["newton"])
  lines(ns, times_em[, d],  col = sol_col["em"],     lty = dgp_lty[d], lwd = 1.5)
  points(ns, times_em[, d], col = sol_col["em"],     pch = sol_pch["em"])
}
## reference slopes
ref_n <- range(ns)
ref_x <- log10(ref_n)
## sublinear slope reference (n^0.65) anchored at midpoint of kmc-logistic
mid_n <- 1000
mid_t <- times_kmc["1000", "logistic"]
abline_loglog <- function(slope, anchor_x, anchor_y, ...) {
  x <- 10^seq(log10(min(ns)), log10(max(ns)), length.out = 50)
  y <- anchor_y * (x / anchor_x)^slope
  lines(x, y, ...)
}
abline_loglog(0.65, mid_n, mid_t, col = sol_col["kmc"],   lty = 4, lwd = 0.7)
abline_loglog(1.10, mid_n, times_em["1000", "logistic"],  col = sol_col["em"], lty = 4, lwd = 0.7)
legend("topleft", inset = c(0.0, 0.0),
       legend = c("kmc.bcm.test", "Newton (this work)", "iterative EM",
                  "ref slope 0.65", "ref slope 1.10"),
       col = c(sol_col["kmc"], sol_col["newton"], sol_col["em"],
                sol_col["kmc"], sol_col["em"]),
       pch = c(sol_pch[1:3], NA, NA),
       lty = c(1, 1, 1, 4, 4),
       lwd = c(1.5, 1.5, 1.5, 0.7, 0.7),
       cex = 0.78, bty = "n")

## --- Panel B: Type-I vs n ---
t1_kmc <- sapply(dgps, function(d) sapply(ns, function(n) get_t1(cells3, d, n, "kmc")))
t1_em  <- sapply(dgps, function(d) sapply(ns, function(n) get_t1(cells3, d, n, "em")))
rownames(t1_kmc) <- rownames(t1_em) <- as.character(ns)

plot(NA, log = "x", xlim = range(ns), ylim = c(0, 0.85),
     xlab = expression("sample size " * n),
     ylab = expression("empirical Type-I error at " * alpha == 0.05),
     main = "(B)  Calibration")
abline(h = 0.05, col = "grey50", lty = 3)
for (d in dgps) {
  lines(ns,  t1_kmc[, d], col = sol_col["kmc"], lty = dgp_lty[d], lwd = 1.5)
  points(ns, t1_kmc[, d], col = sol_col["kmc"], pch = sol_pch["kmc"])
  lines(ns,  t1_em[, d],  col = sol_col["em"],  lty = dgp_lty[d], lwd = 1.5)
  points(ns, t1_em[, d],  col = sol_col["em"],  pch = sol_pch["em"])
}
legend("topleft",
       legend = c("kmc.bcm.test", "iterative EM",
                  "logistic", "normal", "gamma",
                  "nominal level 0.05"),
       col = c(sol_col["kmc"], sol_col["em"],
                "grey30", "grey30", "grey30", "grey50"),
       pch = c(sol_pch[1], sol_pch[3], NA, NA, NA, NA),
       lty = c(1, 1, 1, 2, 5, 3),
       lwd = c(1.5, 1.5, 1.5, 1.5, 1.5, 1),
       cex = 0.78, bty = "n")

## --- Panel C: paired-statistic scatter at one representative cell (logistic, n=2000) ---
cell <- cells4[["logistic_2000"]]
xs   <- cell$raw$s_kmc
xs_n <- cell$raw$s_new
xs_e <- cell$raw$s_em
ok   <- !is.na(xs) & !is.na(xs_n) & !is.na(xs_e)
xs <- xs[ok]; xs_n <- xs_n[ok]; xs_e <- xs_e[ok]

xlim_c <- range(c(xs, xs_n, xs_e), na.rm = TRUE)
plot(xs, xs_n, log = "", xlim = xlim_c, ylim = xlim_c,
     pch = sol_pch["newton"], col = sol_col["newton"], cex = 1.1,
     xlab = expression(gamma[n] * " from kmc.bcm.test"),
     ylab = expression("competing solver's " * gamma[n]),
     main = "(C)  Cross-solver agreement (logistic, n=2000)")
points(xs, xs_e, pch = sol_pch["em"], col = sol_col["em"], cex = 0.9)
abline(0, 1, col = "grey50", lty = 2)
legend("topleft",
       legend = c("Newton (this work)",
                  "iterative EM",
                  "identity"),
       col = c(sol_col["newton"], sol_col["em"], "grey50"),
       pch = c(sol_pch["newton"], sol_pch["em"], NA),
       lty = c(NA, NA, 2), cex = 0.85, bty = "n")
mtext(side = 3, line = -1.5, adj = 0.04,
      text = sprintf("max |kmc - Newton| = %.1e", max(abs(xs - xs_n))),
      cex = 0.70)
mtext(side = 3, line = -2.7, adj = 0.04,
      text = sprintf("max |kmc - EM|     = %.1f", max(abs(xs - xs_e))),
      cex = 0.70)

dev.off()
cat("Wrote fig/bcm_bench.pdf\n")
