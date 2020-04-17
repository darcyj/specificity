# RUN FROM THIS FOLDER
packwd <- getwd()


library(devtools)
library(roxygen2)
library(Rcpp)
library(testthat)

usethis::use_package("ape")
usethis::use_package("parallel")
usethis::use_package("Rcpp")
usethis::use_package("fields")
usethis::use_package("fitdistrplus")


# this looks for Rcpp::export in cpp files
Rcpp::compileAttributes("./")

document()

# ADD C++ namespace stuff to NAMESPACE file, because document() above OVERWRITES it:
write("", file="NAMESPACE",append=TRUE)
write("# Rcpp namespace stuff, added by update_r_package.r", file="NAMESPACE", append=TRUE)
write("useDynLib(specificity, .registration=TRUE)", file="NAMESPACE", append=TRUE)
write("exportPattern(\"^[[:alpha:]]+\")", file="NAMESPACE", append=TRUE)
write("importFrom(Rcpp, evalCpp)", file="NAMESPACE", append=TRUE)

# install
setwd("../")
install("specificity")

# test
setwd("specificity/tests/")
### doing tests can take a minute! ###
capture.output(test_check("specificity"), file="../test_results.txt", split=TRUE)

# make pdf manual
# texinfo required: sudo apt-get install texinfo
setwd(packwd)
pack <- "specificity"
path <- find.package(pack)
system("rm specificity.pdf")
# commmand below requires: texlive-fonts-extra, texinfo (install both with apt)
system(paste(shQuote(file.path(R.home("bin"), "R")), "CMD", "Rd2pdf", shQuote(path)))


