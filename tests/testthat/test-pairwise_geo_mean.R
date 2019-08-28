library(specificity)
library(testthat)

testvar <- 0:9

a <- pairwise_geo_mean(testvar)
b <- as.vector(as.dist(outer(X=testvar, Y=testvar, FUN=function(x,y){sqrt(x*y)})))

test_that("pairwise_geo_mean matches outer results", { expect_true(all(a==b))})
