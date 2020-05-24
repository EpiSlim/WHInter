context("WHInter")

test_that("WHInter accepts sparse matrices as entry", {
  require(Matrix)
  
  n <- 100
  p <- 150
  
  set.seed(98765)
  Xmat <- rsparsematrix(n, p, 0.1, rand.x=function(n) rep(1.0, n)) # only 10% of matrix elements are non-zero
  Yvec <- rnorm(n) # equally split between positives and negatives
  
  expect_success(result <- WHInter(Xmat, Yvec))
})
