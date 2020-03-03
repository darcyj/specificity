# RUN FROM THIS FOLDER
packwd <- getwd()


library(devtools)
library(roxygen2)
library(Rcpp)
library(testthat)

usethis::use_package("ape")
usethis::use_package("geiger")
usethis::use_package("parallel")
usethis::use_package("Rcpp")
usethis::use_package("fields")

# this looks for Rcpp::export in cpp files
Rcpp::compileAttributes("./")

document()

# ADD C++ namespace stuff to NAMESPACE file, because document() above OVERWRITES it:
write("", file="NAMESPACE",append=TRUE)
write("# Rcpp namespace stuff, added by update_r_package.r", file="NAMESPACE", append=TRUE)
write("useDynLib(specificity, .registration=TRUE)", file="NAMESPACE", append=TRUE)
write("exportPattern(\"^[[:alpha:]]+\")", file="NAMESPACE", append=TRUE)
write("importFrom(Rcpp, evalCpp)", file="NAMESPACE", append=TRUE)


setwd("../")
install("specificity")
setwd("specificity/tests/")
### doing tests can take a minute! ###
capture.output(test_check("specificity"), file="../test_results.txt", split=TRUE)

# make pdf manual
# texinfo required: sudo apt-get install texinfo
setwd(packwd)
pack <- "specificity"
path <- find.package(pack)
system("rm specificity.pdf")
system(paste(shQuote(file.path(R.home("bin"), "R")), "CMD", "Rd2pdf", shQuote(path)))


