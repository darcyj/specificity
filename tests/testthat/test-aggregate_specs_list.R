library(specificity)
library(testthat)

# make data
data(endophyte)
otutable <- prop_abund(endophyte$otutable)
otutable <- occ_threshold(otutable, 30)

# define n cores parameter
ncores <- 2 # for CRAN

# make data - use index_rough to make it faster 
specs_list <- list()
specs_list$elevation <- phy_or_env_spec(otutable, endophyte$metadata$Elevation, 
	n_sim=100, n_cores=ncores, denom_type="sim_center")
specs_list$rainfall <- phy_or_env_spec(otutable, endophyte$metadata$Rainfall,
	n_sim=100, n_cores=ncores, denom_type="sim_center")




# easy test - right number of rows as output?
	test_that("aggregate_specs_list returns right number of rows with byFeature=FALSE", {
		ag <- aggregate_specs_list(sl=specs_list, byFeature=FALSE,
			fd=endophyte$taxonomy, fd_id=1)
		expect_true(sum(sapply(specs_list, nrow)) == nrow(ag))
	})
	test_that("aggregate_specs_list returns right number of rows with byFeature=TRUE", {
		ag <- aggregate_specs_list(sl=specs_list, byFeature=TRUE,
			fd=endophyte$taxonomy, fd_id=1)
		expect_true(nrow(specs_list$rainfall) == nrow(ag))
	})


# same as above, but with incomplete objects in spec
	test_that("aggregate_specs_list byFeature=TRUE returns right number of rows when sl is incomplete", {
		sl2 <- specs_list
		sl2$elevation <- sl2$elevation[20:nrow(sl2$elevation),]
		ag <- aggregate_specs_list(sl=sl2, byFeature=TRUE,
			fd=endophyte$taxonomy, fd_id=1)
		expect_true(nrow(sl2$rainfall) == nrow(ag))
	})
	test_that("aggregate_specs_list byFeature=FALSE returns right number of rows when sl is incomplete", {
		sl2 <- specs_list
		sl2$elevation <- sl2$elevation[20:nrow(sl2$elevation),]
		ag <- aggregate_specs_list(sl=sl2, byFeature=FALSE,
			fd=endophyte$taxonomy, fd_id=1)
		expect_true(sum(sapply(sl2, nrow)) == nrow(ag))
	})

# now with incomplete featuredata
	fd2 <- endophyte$taxonomy
	fd2 <- fd2[10:nrow(fd2), ]
	test_that("aggregate_specs_list byFeature=FALSE returns right number of rows when fd is incomplete", {
		ag <- aggregate_specs_list(sl=specs_list, byFeature=FALSE,
			fd=fd2, fd_id=1)
		expect_true(sum(sapply(specs_list, nrow)) == nrow(ag))
	})
	test_that("aggregate_specs_list byFeature=TRUE returns right number of rows when fd is incomplete", {
		ag <- aggregate_specs_list(sl=specs_list, byFeature=TRUE,
			fd=fd2, fd_id=1)
		expect_true(nrow(specs_list$rainfall) == nrow(ag))
	})

# make sure it still works when featuredata is re-ordered
	fd1 <- endophyte$taxonomy
	fd3 <- fd1[sample(1:nrow(fd1)), ]
	sl1 <- specs_list
	test_that("aggregate_specs_list byFeature=TRUE is not sensitive to fd ordering", {
		ag_reg <- aggregate_specs_list(sl=sl1, byFeature=TRUE, fd=fd3, fd_id=1)
		ag_perm <- aggregate_specs_list(sl=sl1, byFeature=TRUE, fd=fd1, fd_id=1)

		expect_true(all(ag_reg$tax == ag_perm$tax))
	})
	test_that("aggregate_specs_list byFeature=FALSE is not sensitive to fd ordering", {
		ag_reg <- aggregate_specs_list(sl=sl1, byFeature=FALSE, fd=fd3, fd_id=1)
		ag_perm <- aggregate_specs_list(sl=sl1, byFeature=FALSE, fd=fd1, fd_id=1)

		expect_true(all(ag_reg$tax == ag_perm$tax))
	})
