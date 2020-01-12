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



# testing "MASS_fit" method

set.seed(12345)
perms <- replicate(1000, rnorm(500, mean=0, sd=1))
emp <- -2.2
parametric_P <- dnorm(x=emp, mean=0, sd=1)
calculated_Ps <- apply(X=perms, MAR=2, FUN=function(x){pval_from_perms(emp, x, 1, method="MASS_fit")})
p_fun <- round(mean(calculated_Ps),2)
p_calc <- round(parametric_P, 2)
test_that("pval matches calculated value (1 tail, MASS_fit)", { expect_identical(p_calc, p_fun) })

