context("WHInter")

test_that("WHInter accepts sparse matrices as entry", {
  n <- 100
  p <- 150
  
  set.seed(98765)
  Xmat <- Matrix::rsparsematrix(n, p, 0.1, rand.x=function(n) rep(1.0, n)) # only 10% of matrix elements are non-zero
  Yvec <- rnorm(n) # equally split between positives and negatives
  result <- WHInter(Xmat, Yvec)
  
  expect_equal(names(result), c("bias", "beta", "lambda", "dim", "nobs", "offset"))
  
})
