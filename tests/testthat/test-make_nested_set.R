library(specificity)
library(testthat)

ncores <- ifelse(.Platform$OS.type == "unix", yes=2, no=1) # for CRAN check

# a is a very simple tree
a <- ape::read.tree(text="(((A,B),C),D);")
ans <- make_nested_set(a, n_cores=ncores)

test_that("nested set right order", {
	expect_identical(ans[,1], a$edge[,2])
})

test_that("nested set has congruity for rooted tree", {
	expect_true(all(ans[,4] == 1))
})

test_that("nested set fails for unrooted tree", {
	au <- ape::unroot(a)
	expect_error(make_nested_set(au, n_cores=ncores))
})
