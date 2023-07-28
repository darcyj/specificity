library(specificity)
library(testthat)


test_that("wpd's pd metric matches simple branch sum", {
  # use plant genera tree from specificity::endophyte
  set.seed(12345)
  a <- endophyte$supertree
  s <- sample(a$tip.label, 20)
  wpd_result <- wpd(s,a, metric="PD")
  mat_result <- sum(ape::keep.tip(a,s)$edge.length)
  expect_true(round(wpd_result, 3) == round(mat_result, 3))
})

test_that("wpd Hp matches manual calculation",{
  a <- ape::read.tree(text="((a:2,b:2):1,(c:1,d:1):2);")
  manual_result <- round(2.641777, 3)
  wpd_result <- round(wpd(s=c("a","a","b","c"), a), 3)
  expect_true(manual_result == wpd_result)
})

test_that("wpd Hp matches manual calculation, using w for weights",{
  a <- ape::read.tree(text="((a:2,b:2):1,(c:1,d:1):2);")
  manual_result <- round(2.641777, 3)
  wpd_result <- round(wpd(s=c("a", "b", "c"), a, w=c(2,1,1)), 3)
  expect_true(manual_result == wpd_result)
})


test_that("wpd_table matches manual calculation",{
  a <- ape::read.tree(text="((a:2,b:2):1,(c:1,d:1):2);")
  manual_result <- round(2.641777, 3)
  m <- data.frame(a=c(2,2), b=c(1,1), c=c(1,1))
  wpd_table_result <- round( wpd_table(m, a, ncores=1), 3 )
  expect_true(manual_result == wpd_table_result[1])
})
