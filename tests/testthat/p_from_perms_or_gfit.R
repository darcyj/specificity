library(specificity)
library(testthat)
library(parallel)

set.seed(12345)
perm <- rnorm(1000, mean=50, sd=20)
emp <- 30
ncores <- 2 # for CRAN check

p_calc <- (sum(perm < emp) + 1) / length(perm)
p_fun <- specificity:::p_from_perms_or_gfit(emp, perm, tails=1)
test_that("pval matches calculated value (left tail, perms)", { expect_identical(p_calc, p_fun) })

emp <- 80
p_calc <-  1 - ((sum(perm < emp) + 1) / length(perm))
p_fun <- specificity:::p_from_perms_or_gfit(emp, perm, tails=3)
test_that("pval matches calculated value (right tail, perms)", { expect_identical(p_calc, p_fun) })

p_fun <- specificity:::p_from_perms_or_gfit(emp, perm, tails=2)
test_that("pval matches calculated value (2 tail, raw)", { expect_identical(p_calc, p_fun) })



# testing "gamma_fit" method
set.seed(12345)
perms <- replicate(1000, rgamma(500, shape=10, rate=0.005), simplify=FALSE)
emp <- 1600
parametric_P <- pgamma(q=emp, shape=10, rate=0.005)
calculated_Ps <- simplify2array(mclapply(X=perms,
	FUN=function(x){
		gf <- specificity:::fit_gamma_fwd_rev(x)
		specificity:::p_from_perms_or_gfit(emp, gfit=gf, fallback=1)
	},
	mc.cores=ncores
))
p_fun <- round(mean(calculated_Ps),2)
p_calc <- round(parametric_P, 2)
test_that("pval matches calculated value (1 tail, gamma_fit)", { expect_identical(p_calc, p_fun) })


# testing "gamma_fit" method again, including reversal
set.seed(12345)
perms <- replicate(1000, rgamma(300, shape=runif(1, 2, 14), rate=runif(1, 0.001, 1)), simplify=FALSE)
perms <- lapply(X=perms, FUN=function(x){if(sample(c(T,F), 1)){  min(x) + max(x) - x  }else{x}})
emps <- lapply(X=perms, FUN=function(x){sample(x, size=1)})

perm_ps <- simplify2array(mcmapply(
	FUN=specificity:::p_from_perms_or_gfit,
	emp=emps,
	perm=perms,
	mc.cores=ncores
))

gfs <- mclapply(X=perms, FUN=specificity:::fit_gamma_fwd_rev, mc.cores=10)
gamma_ps <- simplify2array(mcmapply(
	FUN=specificity:::p_from_perms_or_gfit,
	emp=emps,
	gfit=gfs,
	fallback=perm_ps,
	mc.cores=ncores
))

test_that("perm pvals correlate to gamma pvals > 0.95 for p > 0.001", { expect_gt(cor(perm_ps, gamma_ps), 0.95) })
