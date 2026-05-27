## sim/bcm_rank.R
##
## Rank-objective empirical-likelihood test for multinomial choice
## (Route I of the multinomial extension in Section 7 of the chapter).
##
## Random utility model: U_{ij} = beta_j^T x_i + eps_{ij}, j = 1, ..., J-1;
## U_{i0} = 0; y_i = argmax_j U_{ij}.  Test H0: beta = beta_0.
##
## Construction (Horowitz-smoothed multinomial max-score):
##   tilde mu_i(beta; h) = prod_{l != y_i} Phi( (beta_{y_i} - beta_l)^T x_i / h )
## with beta_0 = 0.  This is the smooth indicator that subject i's
## observed choice has the highest predicted utility.  Its gradient
##   g_i(beta) = grad_beta tilde mu_i(beta; h)
## has zero population mean at beta = beta_0 by the first-order
## condition of the smoothed max-score objective.
##
## Crucially this avoids the within-pair selection bias of
## Section 7.1 because every subject contributes through every
## comparison y_i vs l != y_i; no observation is restricted.
##
## Owen-EL on bar g_n with one coordinate dropped for scale
## identification gives a chi^2_{p(J-1)-1}-calibrated test at
## bandwidth h_n asymp n^{-1/5}.

source("sim/bcm_multinomial.R", chdir = FALSE)

EPS_KER <- 1e-12

## ---------------------------------------------------------
## 1. Smoothed max-score score, vectorised
##    Returns n x (p * (J-1)) matrix g.
## ---------------------------------------------------------
smooth_maxscore_score <- function(X, y, betas, h) {
  ## Returns g_i = nabla_beta tilde_mu_i(beta; h)
  ##              = tilde_mu_i * nabla log tilde_mu_i
  ## with tilde_mu_i(beta; h) = prod_{l != y_i} Phi((beta_{y_i} - beta_l)^T x_i / h),
  ## treating beta_0 = 0 as the base.  The full gradient (not the log)
  ## is the one whose population mean vanishes at beta_0 by the
  ## first-order condition of the smoothed max-score objective.
  n  <- nrow(X); p <- nrow(betas); Jm <- ncol(betas); J <- Jm + 1L
  lp <- X %*% betas                       # n x (J-1)
  lp_full <- cbind(0, lp)                 # n x J

  ## first compute tilde_mu_i for each i
  mu  <- rep(1.0, n)
  for (l_val in 0L:(J - 1L)) {
    for (y_val in 0L:(J - 1L)) {
      if (l_val == y_val) next
      rows_y <- which(y == y_val)
      if (length(rows_y) == 0L) next
      diff_lp <- (lp_full[rows_y, y_val + 1L] - lp_full[rows_y, l_val + 1L]) / h
      Phi_v <- pmax(pnorm(diff_lp), EPS_KER)
      mu[rows_y] <- mu[rows_y] * Phi_v
    }
  }

  ## now build the log-gradient first (per-(y_val, l_val) sums)
  log_grad <- matrix(0, n, p * Jm)
  for (y_val in 0L:(J - 1L)) {
    rows_y <- which(y == y_val)
    if (length(rows_y) == 0L) next
    Xy <- X[rows_y, , drop = FALSE]
    for (l_val in 0L:(J - 1L)) {
      if (l_val == y_val) next
      diff_lp <- (lp_full[rows_y, y_val + 1L] - lp_full[rows_y, l_val + 1L]) / h
      phi_v <- dnorm(diff_lp)
      Phi_v <- pmax(pnorm(diff_lp), EPS_KER)
      ratio <- phi_v / Phi_v / h
      term  <- ratio * Xy
      if (y_val >= 1L) {
        cols <- ((y_val - 1L) * p + 1L):(y_val * p)
        log_grad[rows_y, cols] <- log_grad[rows_y, cols] + term
      }
      if (l_val >= 1L) {
        cols <- ((l_val - 1L) * p + 1L):(l_val * p)
        log_grad[rows_y, cols] <- log_grad[rows_y, cols] - term
      }
    }
  }
  ## true gradient: tilde_mu_i * log_grad_i
  mu * log_grad
}

## ---------------------------------------------------------
## 2. Rank-EL statistic for H0: beta = beta_0
##    Drop one scale-identification coordinate (the first overall).
## ---------------------------------------------------------
rank_el_lr <- function(X, y, betas, h = NULL) {
  ## Default bandwidth: under-smoothed relative to Horowitz's
  ## bias-variance-optimal h ~ n^{-1/5}.  Under-smoothing
  ## (h ~ n^{-1/3}) reduces the asymptotic bias of the smoothed
  ## score from O(h^2 sqrt(n)) = O(n^{1/10}) (Horowitz rate, which
  ## would grow with n) down to O(h^2 sqrt(n)) = O(n^{-1/6})
  ## (vanishing).  The variance inflates by O(h^{-1}) = O(n^{1/3})
  ## per kernel evaluation but cancels in the Owen-EL self-normalisation.
  if (is.null(h)) h <- 1.5 * sd(X) * nrow(X)^(-1/3)
  G <- smooth_maxscore_score(X, y, betas, h)
  G <- G[, -1L, drop = FALSE]              # one overall scale identifier
  res <- tryCatch(
    emplik::el.test(G, mu = rep(0, ncol(G))),
    error = function(e) list(`-2LLR` = NA_real_)
  )
  list(stat = res[["-2LLR"]], df = ncol(G), h = h)
}

## ---------------------------------------------------------
## 3. Calibration cell
## ---------------------------------------------------------
run_rank_cell <- function(n, betas, dgp_args, reps, seed, h = NULL) {
  set.seed(seed)
  stats   <- numeric(reps)
  times   <- numeric(reps)
  for (r in seq_len(reps)) {
    d  <- do.call(gen_multinomial, c(list(n = n, betas = betas), dgp_args))
    t0 <- proc.time()
    out <- rank_el_lr(d$X, d$y, betas, h = h)
    times[r] <- (proc.time() - t0)[3L]
    stats[r] <- out$stat
  }
  list(stats = stats, times = times,
       df = nrow(betas) * ncol(betas) - 1L)
}

summarise_rank <- function(cell) {
  s    <- cell$stats
  df   <- cell$df
  good <- !is.na(s) & is.finite(s) & s >= 0
  s    <- s[good]
  list(
    df          = df,
    n_ok        = sum(good),
    type1_10    = mean(s > qchisq(0.90, df)),
    type1_05    = mean(s > qchisq(0.95, df)),
    type1_01    = mean(s > qchisq(0.99, df)),
    mean_stat   = mean(s),
    median_stat = median(s),
    ks_p        = suppressWarnings(ks.test(s, "pchisq", df = df))$p.value,
    median_time = median(cell$times[good])
  )
}

## ---------------------------------------------------------
## 4. Main driver
## ---------------------------------------------------------
main_rank <- function(reps = 1000L,
                      out_path = "sim/results/rank_results.rds") {
  cat(sprintf("\nRank-EL multinomial: reps=%d per cell\n", reps))
  cat(strrep("=", 70), "\n", sep = "")

  betas_p2_J3 <- cbind(c(1.0, 0.5), c(-0.5, 1.0))     # df = 2*2 - 1 = 3
  betas_p2_J4 <- cbind(c(1.0, 0.5), c(-0.5, 1.0), c(0.7, -0.3))  # df = 5
  betas_p3_J3 <- cbind(c(1.0, 0.5, -0.3), c(-0.5, 1.0, 0.2))     # df = 5

  cfgs <- list(
    list(name = "gumbel_p2_J3_n500",  n = 500L,  betas = betas_p2_J3,
         dgp_args = list(eps_dist = "gumbel")),
    list(name = "gumbel_p2_J3_n1000", n = 1000L, betas = betas_p2_J3,
         dgp_args = list(eps_dist = "gumbel")),
    list(name = "gumbel_p2_J3_n2000", n = 2000L, betas = betas_p2_J3,
         dgp_args = list(eps_dist = "gumbel")),
    list(name = "normal_p2_J3_n1000", n = 1000L, betas = betas_p2_J3,
         dgp_args = list(eps_dist = "normal")),
    list(name = "gumbel_p2_J4_n1000", n = 1000L, betas = betas_p2_J4,
         dgp_args = list(eps_dist = "gumbel")),
    list(name = "gumbel_p3_J3_n1000", n = 1000L, betas = betas_p3_J3,
         dgp_args = list(eps_dist = "gumbel"))
  )

  cells <- list()
  for (cfg in cfgs) {
    cat(sprintf("  cell %-22s ... ", cfg$name)); flush.console()
    t0 <- proc.time()
    cell <- run_rank_cell(cfg$n, cfg$betas, cfg$dgp_args, reps,
                          seed = 42L + 17L * cfg$n + nchar(cfg$name))
    el <- (proc.time() - t0)[3L]
    s  <- summarise_rank(cell)
    cat(sprintf("%5.1fs | df=%d T-I(.05)=%.3f mean=%.2f KSp=%.3f\n",
                el, s$df, s$type1_05, s$mean_stat, s$ks_p))
    cells[[cfg$name]] <- list(meta = cfg, raw = cell, summary = s)
  }
  saveRDS(cells, file = out_path)
  cat("\nSaved to ", out_path, "\n", sep = "")
  invisible(cells)
}

if (sys.nframe() == 0L) {
  main_rank()
}
