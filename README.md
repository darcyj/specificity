# specificity
R package for analyzing specificity in ecological data

by John L. Darcy
Updated 29 Jan 2020

## What is "specificity"?
The word "specificity" is often used in the context of host-parasite interactions. For example, a parasite may be "specific" to a narrow phylogenetic group of hosts. This is phylogenetic specificity. However, in this R package I define specificity more broadly: a species has specificity to a variable when that species occupies a limited range of that variable. A species may have specificity to elevation if it can only be found within a narrow elevational band. With this software, you can compute specificity for many species at once, to different types of variables including 1-dimensional variables (e.g. elevation), 2-dimensional variables (e.g. geographic distance), and phylogenies. This R package can be used to analyze patterns of ecological specificity, particularly in microbiome data. Please check out the tutorial vignette in the `vignette` folder of this repository, which will walk you through the use of this software. 

## Installation
To install this package, make sure that R has package `devtools` installed. You can install it by running `install.packages("devtools")`. Then, load devtools using `library("devtools")` and then install specificity using `install_github("darcyj/specificity")`. 

## Requirements
specificity needs a few R packages to run. They should be installed automatically when you run the above. They are: `ape`, `geiger`, `fitdistrplus`, and `fields`. You'll also need R, obviously. All tests are run on R 3.6.2, but the code will probably run fine on earlier versions. For parallel processing you'll need to be using Linux (or maybe OS X?), since R's parallelism doesn't work well on Windows. If you're using Windows, just set `n_cores` to 1 for any given function. 

## Documentation
Documentation for specific functions can be obtained within R using `?fun` where fun is the function you're interested in. There is also a full PDF documentation file in this repository named `specificity.pdf`. 

## Contact
Send me an email if you'd like.