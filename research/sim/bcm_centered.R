## sim/bcm_centered.R
##
## Implements the centered estimating function for which (A7) holds
## rigorously:
##
##   tilde g_F (delta, X) = ( X - hat m(e) ) * hat eps(delta, X; F)
##
## where hat m(e) is a Nadaraya-Watson kernel regression of X on
## e = -beta_0^T X.  The first-order plug-in bias vanishes
## identically for tilde g by the tower property,
##
##   E[ (X - E[X|e]) K0(e) (F - F0)(e) ]
##     = E[ E[X - E[X|e] | e] * K0(e) (F-F0)(e) ]
##     = E[ 0 * K0 (F-F0) ] = 0,
##
## so the Owen-EL test based on tilde g is asymptotically chi^2_{p-1}
## without (A7) appearing as an extra hypothesis.

source("sim/bcm_sim.R", chdir = FALSE)

## ---------------------------------------------------------
## Nadaraya-Watson kernel regression  hat m(e) = E[X | e]
## ---------------------------------------------------------
nw_kernel <- function(X, e, h = NULL) {
  n <- length(e)
  if (is.null(h)) {
    h <- 1.06 * sd(e) * n^(-1/5)
  }
  ## Fully vectorised: build n x n weight matrix once
  ##   W[i,j] = K((e_i - e_j)/h)
  D <- outer(e, e, "-") / h
  W <- exp(-0.5 * D * D)                # N(0,1) kernel up to constant
  Wn <- W / pmax(rowSums(W), 1e-12)     # normalise rows
  Wn %*% X                              # n x p
}

## ---------------------------------------------------------
## Centered BJ-residual estimating function and test statistic
## ---------------------------------------------------------
bcm_lr_centered <- function(X, delta, beta0, h = NULL) {
  e   <- -as.numeric(X %*% beta0)
  bj  <- bj_residuals(X, delta, beta0)
  ord <- bj$ord
  hat_eps <- bj$hat_eps                  # sorted by e
  Xs      <- bj$Xs                       # sorted rows of X
  es      <- sort(e)
  mhat    <- nw_kernel(Xs, es, h = h)    # n x p
  g  <- (Xs - mhat) * hat_eps            # element-wise: each row times scalar
  ## scale identification: drop the first coordinate
  if (ncol(g) >= 2L) g <- g[, -1L, drop = FALSE]
  res <- tryCatch(
    emplik::el.test(g, mu = rep(0, ncol(g))),
    error = function(e) list(`-2LLR` = NA_real_)
  )
  list(stat = res[["-2LLR"]], df = ncol(g))
}

## ---------------------------------------------------------
## Calibration: chi^2_{p-1} under H0  (compare to uncentered)
## ---------------------------------------------------------
run_cell_pair <- function(n, beta_true, dgp, reps, seed) {
  set.seed(seed)
  stat_un <- stat_cn <- numeric(reps)
  gbar_un <- gbar_cn <- numeric(reps)
  for (r in seq_len(reps)) {
    d  <- dgp(n, beta_true)
    o1 <- bcm_lr_owen(d$X, d$delta, beta_true)
    o2 <- bcm_lr_centered(d$X, d$delta, beta_true)
    stat_un[r] <- o1$stat
    stat_cn[r] <- o2$stat
    ## also record the (mean of g) for each
    bj <- bj_residuals(d$X, d$delta, beta_true)
    g_un <- bj$g[, -1L, drop = FALSE]
    gbar_un[r] <- mean(g_un)
    e <- sort(-as.numeric(d$X %*% beta_true))
    mhat <- nw_kernel(bj$Xs, e)
    g_cn <- ((bj$Xs - mhat) * bj$hat_eps)[, -1L, drop = FALSE]
    gbar_cn[r] <- mean(g_cn)
  }
  list(stat_un = stat_un, stat_cn = stat_cn,
       gbar_un = gbar_un, gbar_cn = gbar_cn,
       df = length(beta_true) - 1L)
}

summarise_pair <- function(cell) {
  df <- cell$df
  fmt <- function(s) {
    s <- s[!is.na(s) & is.finite(s) & s >= 0]
    list(
      mean = mean(s),
      t1   = mean(s > qchisq(0.95, df)),
      ks_p = suppressWarnings(ks.test(s, "pchisq", df = df))$p.value
    )
  }
  un <- fmt(cell$stat_un)
  cn <- fmt(cell$stat_cn)
  list(
    df = df,
    un = un, cn = cn,
    sqrtn_bias_un = sqrt(length(cell$stat_un)) * mean(cell$gbar_un, na.rm = TRUE),
    sqrtn_bias_cn = sqrt(length(cell$stat_cn)) * mean(cell$gbar_cn, na.rm = TRUE)
  )
}

main_centered <- function(reps = 1000L,
                          out_path = "sim/results/centered_results.rds") {
  cat(sprintf("\nBCM centered vs uncentered Owen-EL: reps=%d per cell\n", reps))
  cat(strrep("=", 70), "\n", sep = "")

  beta_p2 <- c(1.0,  0.5)
  beta_p3 <- c(1.0,  0.5, -0.3)
  cfgs <- list(
    list(name = "logistic_p2_n500",  n = 500L,  beta = beta_p2, dgp = gen_logistic),
    list(name = "logistic_p2_n1000", n = 1000L, beta = beta_p2, dgp = gen_logistic),
    list(name = "logistic_p3_n500",  n = 500L,  beta = beta_p3, dgp = gen_logistic)
  )
  cells <- list()
  for (cfg in cfgs) {
    cat(sprintf("  cell %-22s ... ", cfg$name)); flush.console()
    t0 <- proc.time()
    cell <- run_cell_pair(cfg$n, cfg$beta, cfg$dgp, reps,
                          seed = 9092L + 11L * cfg$n + nchar(cfg$name))
    el <- (proc.time() - t0)[3L]
    s <- summarise_pair(cell)
    cat(sprintf("%6.1fs\n", el))
    cat(sprintf("    uncentered: mean=%.3f  T-I=%.3f  KSp=%.3f  sqrt(n) bias=%.3f\n",
                s$un$mean, s$un$t1, s$un$ks_p, s$sqrtn_bias_un))
    cat(sprintf("    centered:   mean=%.3f  T-I=%.3f  KSp=%.3f  sqrt(n) bias=%.3f\n",
                s$cn$mean, s$cn$t1, s$cn$ks_p, s$sqrtn_bias_cn))
    cells[[cfg$name]] <- list(meta = cfg, raw = cell, summary = s)
  }
  saveRDS(cells, file = out_path)
  invisible(cells)
}

if (sys.nframe() == 0L) {
  main_centered()
}
