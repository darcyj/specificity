#include <Rcpp.h>
using namespace Rcpp;
//
// THE FOLLOWING DOCUMENTATION IS AVAILABLE IN ../man/pairwise_geo_mean.Rd
// pairwise_geo_mean
//
// Calculates pairwise geometric means from unique 2-element
// combinations of vector x. Written in C++ because R slow. 
// The output vector is the same length and same order
// as a lower triangle of matrix with rows and columns x.
//
// @author John L. Darcy
// 
// @param x numeric vector.
//
// @return vector of pairwise geometric means, of length (l*l-l)/2,
//   where l=length(x).

// [[Rcpp::export]]
NumericVector pairwise_geo_mean(const NumericVector x) {
	// get number of pairwise points (npp) and initialize output vector
	int n = x.size();
	int npp = ((n * n) - n) / 2;
	NumericVector output(npp);
	//Counter just counts where we are in output.
	int counter = 0;
	int nm1 = n - 1;
	int jp1 = 0;
	for(int j = 0; j < nm1; j++){ 
		jp1 = j + 1;
		for(int i = jp1; i < n; i++){
			//calculate geometric mean of x[i] and x[j], put it in output[counter]
			output[counter] = sqrt(x[i] * x[j]);
			// increment counter so the next time the loop goes round,
			// it will write to the next spot.
			counter++;
		}
	}
	return output;
}
