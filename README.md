# specificity
R package for analyzing specificity in ecological data

by John L. Darcy
Updated 15 MAY 2020

## What is "specificity"?
The word "specificity" is often used in the context of host-parasite interactions. For example, a parasite may be "specific" to a narrow phylogenetic group of hosts. This is phylogenetic specificity. However, in this R package I define specificity more broadly: a species has specificity to a variable when that species occupies a limited range of that variable. A species may have specificity to elevation if it can only be found within a narrow elevational band. With this software, you can compute specificity for many species at once, to different types of variables including 1-dimensional variables (e.g. elevation), 2-dimensional variables (e.g. geographic distance), and phylogenies. This R package can be used to analyze patterns of ecological specificity, particularly in microbiome data. Please check out the tutorial vignette in the `vignette` folder of this repository, which will walk you through the use of this software. 

## Installation
To install this package, make sure that R has package `remotes` installed. It's part of `devtools` so if you have that, you should be good. You can install it by running `install.packages("remotes")`. Then, load `remotes` using `library("remotes")` and then install `specificity` using `install_github("darcyj/specificity")`. 

## Requirements
`specificity` needs a few R packages to run. They should be installed automatically when you run the above. `Rcpp` will be installed in order to compile the C++ code in `specificity`. For users on Apple OS X, you may be asked by Rstudio to install "command line developer tools". If you're on OS X and not using Rstudio, you'll need to run `xcode-select --install` from terminal beforehand. If installing from Rstudio and asked whether to compile from source, you'll get best results if you say "yes". Packages `ape` and `fields` will also be installed automatically. For parallel processing you'll need to be using Linux or OS X, since R's parallelism doesn't work well on Windows. If you're using Windows, just set `n_cores` to 1 for any given paralellized function. OR if you're on Windows 10, running R through Windows Subsystem for Linux may give better multithreading results - please let me know if you try this.

## Tutorial Vignette
A full tutorial for `specificity` can be found in the `vignette` folder of this repository. To find it, just go to [github.com/darcyj/specificity](https://github.com/darcyj/specificity) and click on the folder called "vignette". In there you'll find both a .pdf (for everybody) and .rmd (for rstudio users) that will walk you through the package.

## Documentation
Documentation for specific functions can be obtained within R using `?fun` where fun is the function you're interested in. There is also a full .pdf documentation file in this repository named `specificity.pdf`, which is a standard glossary of all functions in the package.

## Contact
Send me an email if you'd like.