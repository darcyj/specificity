library(specificity)
library(testthat)

testvar <- 0:9

a <- pairwise_product(testvar)
b <- as.vector(as.dist(outer(X=testvar, Y=testvar, FUN=function(x,y){x*y})))

test_that("pairwise_product matches outer results", { expect_true(all(a==b))})




