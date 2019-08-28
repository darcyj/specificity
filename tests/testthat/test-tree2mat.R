library(specificity)
library(testthat)


test_that("tree2mat matches ape::cophenetic.phylo", {
	# use plant genera tree from specificity::endophyte
	set.seed(12345)
	a <- endophyte$supertree
	atips <- sample(a$tip.label, 20)
	cph_mat <- ape::cophenetic.phylo(ape::keep.tip(a, atips))
	# re-order cph_mat to match atips
	cph_mat <- cph_mat[order(match(rownames(cph_mat),atips)), order(match(colnames(cph_mat),atips))]
	cph_dis <- as.dist(cph_mat)
	t2m_dis <- tree2mat(a, atips)
	# round to 4 decimal places gets around C vs R precision balogna
	expect_true(all(round(cph_dis, 4) == round(t2m_dis, 4)))
})

test_that("tree2mat error if tip not in tree", {
	expect_error(tree2mat(endophyte$supertree, c("bob", "Cyanea", "Euphorbia")))
})

test_that("tree2mat error if tip in tree multiple times", {
	a <- endophyte$supertree
	a$tip.label[1:4] <- "bad_label"
	expect_error(tree2mat(a, c("Lycium", "Cyanea", "Euphorbia", "bad_label")))
})
