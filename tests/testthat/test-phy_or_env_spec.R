library(specificity)
library(testthat)

data(endophyte)
attach(endophyte)

# make some OTUs to test.
# flat = no pattern
z_flat <- rep(1, length(metadata$Elevation))
# lo = VERY specific
z_elev_lo <- rep(0, length(metadata$Elevation))
z_elev_lo[metadata$Elevation > 1000 & metadata$Elevation < 1200] <- 200
# med = KINDA specific
z_elev_med <- rep(0, length(metadata$Elevation))
z_elev_med[metadata$Elevation > 850 & metadata$Elevation < 1350] <- 200
# hi = LESS specific
z_elev_hi <- rep(0, length(metadata$Elevation))
z_elev_hi[metadata$Elevation > 700 & metadata$Elevation < 1500] <- 200
# lo_all = lo but present at all sites. same for med_all and hi_all
z_elev_lo_all <- z_elev_lo + 2
z_elev_med_all <- z_elev_med + 2
z_elev_hi_all <- z_elev_hi + 2

# par(mfrow=c(3,2))
# plot(z_elev_lo      ~ metadata$Elevation, ylim=c(0, 13))
# plot(z_elev_lo_all  ~ metadata$Elevation, ylim=c(0, 13))
# plot(z_elev_med     ~ metadata$Elevation, ylim=c(0, 13))
# plot(z_elev_med_all ~ metadata$Elevation, ylim=c(0, 13))
# plot(z_elev_hi      ~ metadata$Elevation, ylim=c(0, 13))
# plot(z_elev_hi_all  ~ metadata$Elevation, ylim=c(0, 13))
#     

# calculate specificities of fake data above
spec_results <- phy_or_env_spec(
	abunds_mat = data.frame(z_flat, z_elev_lo, z_elev_med, z_elev_hi,
		z_elev_lo_all, z_elev_med_all, z_elev_hi_all),
	env = metadata$Elevation,
	n_sim = 100,
	n_cores = 10,
	verbose=FALSE,
	p_method="raw" , diagnostic = T
)

test_that("flat otu Specificity is zero", {
	expect_identical(0, round(spec_results$Spec[1], 4))
})

test_that("flat otu pval is 1", {
	expect_identical(1, spec_results$Pval[1])
})
	# sim <- spec_results[1, 7:ncol(spec_results)] 
	# emp <- spec_results[1,4]


test_that("Spec follows lo < med < hi sims", { expect_true(all(c(
	spec_results$Spec[2] < spec_results$Spec[3],
	spec_results$Spec[3] < spec_results$Spec[4]
)))})

test_that("Spec follows lo < med < hi sims with 1 added", {expect_true(all(c(
	spec_results$Spec[5] < spec_results$Spec[6],
	spec_results$Spec[6] < spec_results$Spec[7]
)))})

test_that("Spec sensitive to occupancy as expected",{expect_true(all(c(
	spec_results$Spec[2] < spec_results$Spec[5],
	spec_results$Spec[3] < spec_results$Spec[6],
	spec_results$Spec[4] < spec_results$Spec[7]
)))})

test_that("non-square env matrix gives error", {expect_error(
	crap <- phy_or_env_spec(
		abunds_mat = data.frame(z_flat, z_elev_lo, z_elev_med, z_elev_hi,
			z_elev_lo_all, z_elev_med_all, z_elev_hi_all),
		env = cbind(metadata$Elevation, metadata$Elevation), # nx2 matrix, not nxn.
		n_sim = 1000,
		n_cores = 10
	)
)})

test_that("wrong length env vector gives error", {expect_error(
	crap <- phy_or_env_spec(
		abunds_mat = data.frame(z_flat, z_elev_lo, z_elev_med, z_elev_hi,
			z_elev_lo_all, z_elev_med_all, z_elev_hi_all),
		env = c(metadata$Elevation, 0),# length n+1 not n
		n_sim = 1000,
		n_cores = 10
	)
)})

test_that("wrong length env dist gives error", {expect_error(
	crap <- phy_or_env_spec(
		abunds_mat = data.frame(z_flat, z_elev_lo, z_elev_med, z_elev_hi,
			z_elev_lo_all, z_elev_med_all, z_elev_hi_all),
		env = dist(c(metadata$Elevation, 0)),# length n+1 not n
		n_sim = 1000,
		n_cores = 10
	)
)})

test_that("non-numeric env matrix gives error", {expect_error(
	crap <- phy_or_env_spec(
		abunds_mat = data.frame(z_flat, z_elev_lo, z_elev_med, z_elev_hi,
			z_elev_lo_all, z_elev_med_all, z_elev_hi_all),
		env = matrix("A", ncol=length(z_flat), nrow=length(z_flat)),
		n_sim = 1000,
		n_cores = 10
	)
)})

test_that("non-numeric env vector gives error", {expect_error(
		crap <- phy_or_env_spec(
			abunds_mat = data.frame(z_flat, z_elev_lo, z_elev_med, z_elev_hi,
				z_elev_lo_all, z_elev_med_all, z_elev_hi_all),
			env = rep("A", length(z_flat)),
			n_sim = 1000,
			n_cores = 10
		)
)})

test_that("hosts but no phylo gives error", {expect_error(
	crap <- phy_or_env_spec(
		abunds_mat = data.frame(z_flat, z_elev_lo, z_elev_med, z_elev_hi,
			z_elev_lo_all, z_elev_med_all, z_elev_hi_all),
		hosts=metadata$PlantGenus,
		hosts_phylo = NULL,
		n_sim = 1000,
		n_cores = 10
	)
)})

test_that("phylo but no hosts gives error", {expect_error(
	crap <- phy_or_env_spec(
		abunds_mat = data.frame(z_flat, z_elev_lo, z_elev_med, z_elev_hi,
			z_elev_lo_all, z_elev_med_all, z_elev_hi_all),
		hosts=NULL,
		hosts_phylo = supertree,
		n_sim = 1000,
		n_cores = 10
	)
)})

test_that("hosts+env+phylo gives error", {expect_error(
	crap <- phy_or_env_spec(
		abunds_mat = data.frame(z_flat, z_elev_lo, z_elev_med, z_elev_hi,
			z_elev_lo_all, z_elev_med_all, z_elev_hi_all),
		hosts=metadata$PlantGenus,
		hosts_phylo = supertree,
		env=metadata$Elevation,
		n_sim = 1000,
		n_cores = 10
	)
)})

test_that("hosts+env gives error", {expect_error(
	crap <- phy_or_env_spec(
		abunds_mat = data.frame(z_flat, z_elev_lo, z_elev_med, z_elev_hi,
			z_elev_lo_all, z_elev_med_all, z_elev_hi_all),
		hosts=metadata$PlantGenus,
		env=metadata$Elevation,
		n_sim = 1000,
		n_cores = 10
	)
)})

test_that("env+phylo gives error", {expect_error(
	crap <- phy_or_env_spec(
		abunds_mat = data.frame(z_flat, z_elev_lo, z_elev_med, z_elev_hi,
			z_elev_lo_all, z_elev_med_all, z_elev_hi_all),
		hosts_phylo = supertree,
		env=metadata$Elevation,
		n_sim = 1000,
		n_cores = 10
	)
)})

test_that("incomplete phylo gives error", {expect_error(
	crap <- phy_or_env_spec(
		abunds_mat = data.frame(z_flat, z_elev_lo, z_elev_med, z_elev_hi,
			z_elev_lo_all, z_elev_med_all, z_elev_hi_all),
		hosts=metadata$PlantGenus,
		hosts_phylo = ape::drop.tip(supertree, metadata$PlantGenus[1]),
		n_sim = 1000,
		n_cores = 10
	)
)})

test_that("hosts wrong length gives error", {expect_error(
	crap <- phy_or_env_spec(
		abunds_mat = data.frame(z_flat, z_elev_lo, z_elev_med, z_elev_hi,
			z_elev_lo_all, z_elev_med_all, z_elev_hi_all),
		hosts=c(metadata$PlantGenus, metadata$PlantGenus[1]),
		hosts_phylo = supertree,
		n_sim = 1000,
		n_cores = 10
	)
)})

badhosts <- metadata$PlantGenus
badphylo <- supertree
badphylo$tip.label[badphylo$tip.label == badhosts[1]] <- paste(";", badhosts[1], sep="")
badhosts[badhosts == badhosts[1]] <- paste(";", badhosts[1], sep="")
test_that("semicolons in hosts gives error", {expect_error(
	crap <- phy_or_env_spec(
		abunds_mat = data.frame(z_flat, z_elev_lo, z_elev_med, z_elev_hi,
			z_elev_lo_all, z_elev_med_all, z_elev_hi_all),
		hosts=badhosts,
		hosts_phylo = badphylo,
		n_sim = 1000,
		n_cores = 10
	)
)})

