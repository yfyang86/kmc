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

test_that("kmc works", {
  expect_equal(UNITEST_kmc(), "-17.51983")
  
  expect_equal(UNITEST_kmcbj(), "0.20148")
})

