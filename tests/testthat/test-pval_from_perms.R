library(specificity)
library(testthat)

set.seed(12345)
perm <- rnorm(1000, mean=50, sd=20)
emp <- 30

p_calc <- sum(perm < emp) / length(perm)
p_fun <- pval_from_perms(emp, perm, 1)


test_that("pval matches calculated value (left tail)", { expect_identical(p_calc, p_fun) })

emp <- 80
p_calc <- sum(perm > emp) / length(perm)
p_fun <- pval_from_perms(emp, perm, 3)
test_that("pval matches calculated value (right tail)", { expect_identical(p_calc, p_fun) })

p_fun <- pval_from_perms(emp, perm, 2)
test_that("pval matches calculated value (2 tail)", { expect_identical(p_calc, p_fun) })
