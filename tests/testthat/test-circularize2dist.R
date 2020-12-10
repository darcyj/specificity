library(specificity)
library(testthat)



months <- c(1, 4, 11)
cdm_months <- as.matrix(circularize2dist(months, 12))
colnames(cdm_months) <- rownames(cdm_months) <- months

test_that("diag of cdist is 0 when converted to dm", {
	expect_true(all(diag(cdm_months)==0))
})

test_that("circular dist from jan to nov is 2 months", {
	expect_true( cdm_months["11", "1"] == 2 )
})

test_that("circular dist from jan to apr is 3 months", {
	expect_true( cdm_months["4", "1"] == 3 )
})

test_that("circular dist from nov to apr is 5 months", {
	expect_true( cdm_months["4", "11"] == 5 )
})

# check errors

test_that("non-numeric x causes error", {
	expect_error( circularize2dist(letters[1:5], 5) )
})
test_that("non-numeric maxx causes error", {
	expect_error( circularize2dist(1:10, "5") )
})
test_that("non-scalar maxx causes error", {
	expect_error( circularize2dist(1:10, 11:15) )
})
test_that("non-vector x causes error", {
	expect_error( circularize2dist(matrix(1:9, ncol=3), 100) )
})
test_that("maxx < x causes error", {
	expect_error( circularize2dist(1:20, 5) )
})