// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::depends("RcppParallel")]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::export]]
double weighted_mean(const arma::vec x, const arma::vec w) {
	// x : vector of distances
	// w : vector of weights for x, same length as x
	// use armadillo (BLAS) element-wise multiplication (% operator)
	// use armadillo sum function as well
	return( sum(x % w) / sum(w) );
}

// [[Rcpp::export]]
arma::vec pairwise_geo_mean(const arma::vec x) {
	// get number of pairwise points (npp) and initialize output vector
	int n = x.size();
	int npp = ((n * n) - n) / 2;
	arma::vec output(npp);
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

// [[Rcpp::export]]
double dist_spec(const arma::vec w, const arma::vec x_pair) {
	// this function computes specificity of abundances w to distances x_pair.
	// x_pair MUST BE of length ((n^2)/2)-n, where n=w.size().
	// THIS IS NOT CHECKED for performance reasons.
	// the ORDER of x_pair and w is important as well. Since x_pair
	// is basically a lower triangle, the rows of that lower triangle
	// must correspond to w!!!! No way to check this, just don't be dumb.
	// The inputs are column vectors, but via RcppArmadillo magic, they
	// can 1D numeric arrays. Tested: numeric vector, dist.
	// first, convert w to pairwise weights:
	arma::vec w_pair = pairwise_geo_mean(w);
	// now calculate weighted mean x_pair,w_pair
	return(weighted_mean(x_pair, w_pair));
}

struct Speccer : public Worker {
	// source data
	const RMatrix<double> tmp_mat_w;
	const RVector<double> tmp_vec_x;
	// destination data
	RVector<double> tmp_vec_out;
	// dimensions are passed to the worker so they only get calculated once
	std::size_t nrow_w, ncol_w, size_x;



	// initialize stuff, not sure how this works or what it really does.
	Speccer(const NumericMatrix& mat_w, const NumericVector& vec_x, NumericVector& vec_out,
			std::size_t nrow_w, std::size_t ncol_w, std::size_t size_x )
		: tmp_mat_w(mat_w), tmp_vec_x(vec_x), tmp_vec_out(vec_out),
			nrow_w(nrow_w), ncol_w(ncol_w), size_x(size_x) {}

	// convert RVector/RMatrix into arma type for Rcpp function
	// and the follwing arma data will be shared in parallel computing
	arma::mat RMat2ArmaMat(){
		RMatrix<double> cmat_w = tmp_mat_w;
		arma::mat outmat(cmat_w.begin(), nrow_w, ncol_w, false);
		return outmat;
	}
	arma::vec RVec2ArmaVec(){
		RVector<double> cvec_x = tmp_vec_x;
		arma::vec outvec(cvec_x.begin(), size_x, false);
		return outvec;
	}


	// function call to do the math
	void operator()(std::size_t begin, std::size_t end) {
		//convert (pointers only, not copying anything since copy_aux_mem=false)
		arma::mat MAT_W = RMat2ArmaMat();
		arma::vec VEC_X = RVec2ArmaVec();
		for(std::size_t j = begin; j < end; j++){
			tmp_vec_out[j] = dist_spec(MAT_W.col(j), VEC_X); //temp guts
		}


	}
};

// [[Rcpp::export]]
NumericVector parallelSpecificity(NumericMatrix mat_w, NumericVector vec_x) {
	// calculate dimensions
	std::size_t nrow_w = mat_w.nrow();
	std::size_t ncol_w = mat_w.ncol();
	std::size_t size_x = vec_x.size();


	// allocate the output matrix
	NumericVector vec_out(mat_w.ncol());
	// create the worker (lower case)
	Speccer speccer(mat_w, vec_x, vec_out, nrow_w, ncol_w, size_x);
	// call parallelFor to do the work
	parallelFor(0, mat_w.ncol(), speccer);
	// return the output matrix
	return vec_out;
}
