library(specificity)
library(ape)
library(testthat)

# make ontology with 6 tips and 2 interior nodes
onto <- data.frame(
	l1 = c(9, 9, 9, 9, 9, 9),
	l2 = c(7, 7, 7, 8, 8, 8),
	l3 = c(1, 2, 3, 4, 5, 6)
)

tree <- ape::read.tree(text=onto2nwk(onto))

test_that("test onto makes tree with 6 tips", {
	expect_true(length(tree$tip.label) == 6)
})

test_that("test onto makes rooted tree", {
	expect_true(ape::is.rooted(tree))
})

# make ontology with 6 tips and 2 interior nodes,
	# but some tips are polyphyletic
onto_poly <- data.frame(
	l1 = c(9, 9, 9, 9, 9, 9),
	l2 = c(7, 7, 7, 8, 8, 8),
	l3 = c(1, 2, 3, 4, 5, 1)
)
tree_poly <- ape::read.tree(text=onto2nwk(onto_poly))


test_that("test onto can produce polyphyletic tree", {
	expect_true(sum(tree_poly$tip.label==1)==2)
})


onto_unrooted <- data.frame(
	l2 = c(7, 7, 7, 8, 8, 8),
	l3 = c(1, 2, 3, 4, 5, 6)
)
tree_unrooted_but_not_really <- ape::read.tree(text=onto2nwk(onto_unrooted))
test_that("unrooted ontos are rooted", {
	expect_true(tree_unrooted_but_not_really$Nnode == 3)
})