UNITEST_kmc <- function(){
  x <- c( 1, 1.5, 2, 3, 4.2, 5.0, 6.1, 5.3, 4.5, 0.9, 2.1, 4.3) # positive time
  d <- c( 1,   1, 0, 1, 0, 1, 1, 1, 1, 0, 0,   1)               # status censored/uncensored
  #### compute e-value and its adjustment ####
  g=list( f=function(x) { x-3.7} )
  result = kmc.solve( x,d,g)
  
  return(sprintf('%0.5f', result["loglik.null"]));
}

UNITEST_kmcbj <- function(){
  library(survival)
  stanford5 <- stanford2[!is.na(stanford2$t5), ]
  
  y=log10(stanford5$time)
  d <- stanford5$status
  oy = order(y,-d)
  d=d[oy]
  y=y[oy]
  x=cbind(1,stanford5$age)[oy,]
  beta0 = c(3.2, -0.015)
  result = kmc.bjtest(y, d, x=x, beta = beta0, 
                      init.st="naive")[["-2LLR"]]
  return(sprintf('%0.5f', result));
}

UNITEST_kmcbcm <- function() {
  set.seed(42)
  n     <- 200; beta0 <- c(1.0, 0.5)
  X     <- matrix(runif(n * 2, -1, 1), n, 2)
  eps   <- rlogis(n)
  delta <- as.integer(as.numeric(X %*% beta0) + eps <= 0)
  ## Test at the true beta: -2LLR should be moderate (one MC realisation)
  o_correct <- kmc.bcm.test(X, delta, beta = beta0, centered = FALSE)
  ## Test at a deliberately wrong beta: -2LLR should be larger
  o_wrong   <- kmc.bcm.test(X, delta, beta = c(1, 2), centered = FALSE)
  list(
    correct = sprintf("%.5f", o_correct[["-2LLR"]]),
    wrong   = sprintf("%.5f", o_wrong[["-2LLR"]]),
    correct_df = o_correct$df,
    correct_p  = o_correct$pvalue,
    correct_conv = o_correct$convergence
  )
}

UNITEST_kmcbcm_centered <- function() {
  set.seed(42)
  n     <- 200; beta0 <- c(1.0, 0.5)
  X     <- matrix(runif(n * 2, -1, 1), n, 2)
  eps   <- rlogis(n)
  delta <- as.integer(as.numeric(X %*% beta0) + eps <= 0)
  o <- kmc.bcm.test(X, delta, beta = beta0, centered = TRUE)
  list(stat = o[["-2LLR"]], df = o$df, h = o$h, conv = o$convergence)
}

test_that("kmc works", {
  expect_equal(UNITEST_kmc(), "-17.51983")

  expect_equal(UNITEST_kmcbj(), "0.20148")
})

test_that("kmc.bcm.test (raw) returns sensible values", {
  res <- UNITEST_kmcbcm()
  ## regression-test reference values from v0.4-4 reference build
  expect_equal(res$correct, "0.28456")
  expect_equal(res$wrong,   "1.85650")
  expect_equal(res$correct_df, 1L)
  expect_true(res$correct_p > 0 && res$correct_p < 1)
  expect_equal(res$correct_conv, 1L)
})

test_that("kmc.bcm.test (centered) runs and returns finite numbers", {
  res <- UNITEST_kmcbcm_centered()
  expect_true(is.finite(res$stat) && res$stat >= 0)
  expect_equal(res$df, 1L)
  expect_true(res$h > 0)
  expect_equal(res$conv, 1L)
})

test_that("kmc.bcm.test argument validation works", {
  set.seed(1)
  X <- matrix(runif(40, -1, 1), 20, 2)
  d <- sample(0:1, 20, replace = TRUE)
  expect_error(kmc.bcm.test(X, d, beta = c(1, 0.5, 0.3)),
               "length\\(beta\\)")
  expect_error(kmc.bcm.test(X, d[-1L], beta = c(1, 0.5)),
               "length\\(delta\\)")
  expect_error(kmc.bcm.test(X, d, beta = c(1, 0.5), drop = 3L),
               "drop")
})

