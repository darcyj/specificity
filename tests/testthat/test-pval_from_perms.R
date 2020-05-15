library(specificity)
library(testthat)

set.seed(12345)
perm <- rnorm(1000, mean=50, sd=20)
emp <- 30

p_calc <- (sum(perm < emp) + 1) / length(perm)
p_fun <- pval_from_perms(emp, perm, 1, method="raw")


test_that("pval matches calculated value (left tail, raw)", { expect_identical(p_calc, p_fun) })

emp <- 80
p_calc <- (1+ sum(perm > emp)) / length(perm)
p_fun <- pval_from_perms(emp, perm, 3, method="raw")
test_that("pval matches calculated value (right tail, raw)", { expect_identical(p_calc, p_fun) })

p_fun <- pval_from_perms(emp, perm, 2, method="raw")
test_that("pval matches calculated value (2 tail, raw)", { expect_identical(p_calc, p_fun) })



# testing "gamma_fit" method

set.seed(12345)
perms <- replicate(1000, rgamma(500, shape=10, rate=0.005), simplify=FALSE)
emp <- 1600
parametric_P <- pgamma(q=emp, shape=10, rate=0.005)
calculated_Ps <- simplify2array(mclapply(X=perms,
	FUN=function(x){pval_from_perms(emp, x, 1, method="gamma_fit")},
	mc.cores=10
))
p_fun <- round(mean(calculated_Ps),2)
p_calc <- round(parametric_P, 2)
test_that("pval matches calculated value (1 tail, gamma_fit)", { expect_identical(p_calc, p_fun) })





