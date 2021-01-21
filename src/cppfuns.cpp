#include <Rcpp.h>
#include <random>
using namespace Rcpp;

//
// pairwise_product
//
// Calculates products means from unique 2-element
// combinations of vector x. Written in C++ because R slow. 
// The output vector is the same length and same order
// as a lower triangle of matrix with rows and columns x.
// 
// [[Rcpp::export]]
NumericVector pairwise_product(const NumericVector& x) {
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
			output[counter] = x[i] * x[j];
			// increment counter so the next time the loop goes round,
			// it will write to the next spot.
			counter++;
		}
	}
	return output;
}

// shuffle function to replace std::random_shuffle because R CMD check doesn't
// like rand. lines where std::random_shuffle is used are commented out for easy
// reversion. This function copied *almost* directly from stackoverflow, thx
// justin: https://stackoverflow.com/questions/50243461/stl-random-shuffle-generates-highly-correlated-sequences
// modified by me to have default behavior of only setting seed if seed > 0, 
// otherwise just nondeterministic random
void shuffle_ints_so(IntegerVector &x, int seed=0)
{
	thread_local std::mt19937 engine;

	if(seed > 0){
		engine.seed(seed);
	}

	std::shuffle(x.begin(), x.end(), engine);
}



// rao1sp
//
// Calculates empirical rao for one species. 
//
// [[Rcpp::export]]
float rao1sp(const NumericVector& p, const NumericVector& D, bool perm=false, int seed=0){
	// integer for how big p is (or rather, what the final index of p is)
	int npm1 = p.size() - 1;
	// pidx is a vector of indices for p, so that p doesn't need to
	// be manipulated (copied) in this function. integers are much 
	// more light weight than a vector of floats
	IntegerVector pidx = seq(0, npm1);
	// make rao output
	float rao = 0;
	// if permutation requested, permute pidx
	if(perm){
		// std::random_shuffle(pidx.begin(), pidx.end(); // old way, R CMD check didn't like it
		shuffle_ints_so(pidx, seed);
	}
	// loop through p and D and create rao output
	// i: row of D (as if it were a matrix)
	// j: col of D (as if it were a matrix)
	// k: index of vectorized D
	int k = 0;
	for(int j = 0; j < npm1; j++){
		for(int i = j+1; i <= npm1; i++){
			rao += p[pidx[i]] * p[pidx[j]] * D[k];
			k++;
		}
	}
	return(rao);
}


// rao1sp
//
// Calculates many sim (i.e. null) rao for one species. 
//
// [[Rcpp::export]]
NumericVector raoperms(const NumericVector& p, const NumericVector& D, const int n_sim=1000, int seed=12345){
	NumericVector output(n_sim);
	for(int i = 0; i<n_sim; i++){
		// nobody likes srand :( alternative solution maybe from stackoverflow:
		// https://stackoverflow.com/questions/50243461/stl-random-shuffle-generates-highly-correlated-sequences
		// see comment by Justin and response by Baum mit Augen
		// srand(seed);
		// output[i] = rao1sp(p, D, true);
		output[i] = rao1sp(p, D, true, seed);
		// increment seed to change it up for next time :)
		// only do this if seed > 0. 0 means nondeterministic random.
		if(seed > 0){
			seed++;
		}
	}
	return(output);
}


// spec_core is an internal function used by rao_genetic_max(). It used to be
// used by phy_or_env_spec() too, but now that's done by the much-improved
// rao1sp() and raoperms(). 
// p is a mxn matrix where rows are samples and cols are species (or perms)
// D is a vector of distances
// before use, make SURE that D is length (m(m-1))/2
NumericVector spec_core(const NumericMatrix& p, const NumericVector& D) {
	const int ncol = p.ncol();
	NumericVector output(ncol);

	// for each column of p:
	for(int col = 0; col < ncol; col++){
		NumericVector rowi = p(_,col);
		NumericVector P2 = pairwise_product(rowi);
		// calculate specificity using Rao
		// output[col] = sum(P2 * D) / sum(P2);
		output[col] = sum(P2 * D);
	}
	return output;
}


// rao_sort_max 
//
// approximates the maximum possible rao value given a 
// species weights vector p and a distance vector D. 
//
// [[Rcpp::export]]
long double rao_sort_max(const NumericVector& p, const NumericVector& D) {
	long double out;
	NumericVector P2 = pairwise_product(p);
	std::sort(P2.begin(), P2.end());
	NumericVector D2 = clone(D);
	std::sort(D2.begin(), D2.end());
	out = sum(P2 * D2);
	return out;
}


// function that just swaps two items at random within a vector
inline NumericVector swapcol (const NumericVector x){
	NumericVector xout = clone(x);
	IntegerVector xpos = seq(0, xout.size() - 1);
	shuffle_ints_so(xpos); // seed argument ommitted, default is just to not set a seed
	//std::random_shuffle(xpos.begin(), xpos.end() ); // old way, R CMD check didn't like it
	std::swap(xout[xpos[0]], xout[xpos[1]]);
	return xout;
}

// function to find which top n
inline IntegerVector which_top_n (const NumericVector x, const int n){
	int xs = x.size();
	LogicalVector used(xs);
	IntegerVector out(n);
	long double m = min(x);
	// for each of the n positions we want...
	for(int i=0; i<n; i++){
		long double best = 0 + m;
		// for each value of x, check if its the highest.
		for(int j=0; j<xs; j++){
			if(x[j] > best && used[j] == 0){
				out[i] = j;
				best = x[j];
			}
		}
		// now that we have the best position, cross it off in used
		used[out[i]] = 1;
	}
	return out;

}

// function that just repeats a sequence of ints from:to until it's length len
inline IntegerVector seqrep (const int from, const int to, const int len){
	IntegerVector out(len);
	out[0] = from;
	for(int i=1; i<len; i++){
		if(out[i-1] >= to){
			out[i] = from;
		}else{
			out[i] = out[i-1] + 1;
		}
	}
	return out;
}

// function for testing equality of doubles
inline bool double_equals (const long double x1, const long double x2, const long double eps=0.001){
	return std::abs(x1 - x2) < eps;
}

// genetic optimization algorithm for finding order of p that maximizes Rao
// [[Rcpp::export]]
List rao_genetic_max (const NumericVector p, const NumericVector D, 
	const int term_cycles=10, const int maxiters = 400, 
	const int popsize=300, const int keep=5, const long double prc = 0.001){
	// check to make sure p and D are correct lengths
	int ps = p.size();
	int p_expctd = (ps * (ps - 1)) / 2;
	if(D.size() != p_expctd){
		throw std::invalid_argument("p and D incompatible lengths.");
	}
	// propagate initial matrix; each col is a p vector (a population of p vecs)
	NumericMatrix pop = NumericMatrix(ps, popsize);
	for(int j=0; j<popsize; j++){
		pop(_,j) = swapcol(p); 
	}

	// set up counters and outputs for generations
	long double last_best_rao = 0.0;
	NumericVector iter_best_raos(maxiters, NA_REAL);
	int n_gens_no_improvement = 0;
	int iter = 0;
	bool stopper = 0;

	// some vectors that can be declared here and constantly re-written since 
	// they never change size:
	IntegerVector best_n_pos(keep);                   // stores which cols of pop are best
	NumericVector pop_raos(popsize);                  // rao vals for each col of pop
	NumericMatrix keepmat = NumericMatrix(ps, keep);  // stores best keep cols of pop while pop is overwritten

	// loop through iters
	while(stopper == 0){
		// calculate rao vals for each 
		NumericVector pop_raos = spec_core(pop, D);

		// find best n=keep 
		best_n_pos = which_top_n(pop_raos, keep);

		// record best rao for posterity
		iter_best_raos[iter] = max(pop_raos);
		
		// copy keeper columns into a new object so they don't get overwritten below
		int j = 0;
		for(int k = 0; k < keep; k++){       // k is the column of keepmat
			j = best_n_pos[k];               // j is the column of pop
			keepmat(_,k) = pop(_,j);      // put kth best column of pop into keepmat as column k
		}

		// start filling up pop with new columns
		// do keepers first
		for(int j=0; j<keep; j++){
			pop(_,j) = keepmat(_,j);
		}
		// now do swap cols
		IntegerVector col2do = seqrep(0, keep - 1, popsize);
		NumericVector tempcol(ps);
		for(int j=keep; j<popsize; j++){
			pop(_,j) = swapcol(keepmat(_,col2do[j]));
		}

		// check stop conditions and modify stopper
		if(double_equals(iter_best_raos[iter], last_best_rao, prc)){
			n_gens_no_improvement ++;
		}

		if(n_gens_no_improvement >= term_cycles){
			stopper = 1;
		}

		if(iter >= maxiters){
			stopper = 1;
		}

		// set last_best_rao
		last_best_rao = iter_best_raos[iter];

		// increment iter - MUST BE LAST IN WHILE LOOP!!!
		iter ++;

	}

	// make output list object
	IntegerVector iterations = seq(0, maxiters - 1);
	NumericVector best_p = pop(_,0);
	List ret;
	ret["best_rao"] = last_best_rao;
	ret["iter_raos"] = iter_best_raos;
	ret["iterations"] = iterations;
	ret["best_p"] = best_p;

	return ret;

}
