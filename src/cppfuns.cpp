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
void shuffle_ints_so(IntegerVector &x, int seed=0){
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


// raoperms
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


// makes two different random integers in range [0, max)
// i.e. max-1 is the largest number it can get
// think of this as sampling two positions from an array of length max
// max cannot be less than 2, obviously.
// R package checks don't like rand(), so using R::runif() instead
inline IntegerVector twoDifferentRandomInts(int max){
	if(max < 2){
		stop("max must be >= 2!");
	}
	IntegerVector out = IntegerVector(2);
	out[0] = floor(R::runif(0,1) * max);
	out[1] = floor(R::runif(0,1) * max);
	// commented out because R checks don't like rand()
	// out[0] = rand() % max;
	// out[1] = rand() % max;
	while(out[1] == out[0]){
		// out[1] = rand() % max;
		out[1] = floor(R::runif(0,1) * max);
	}
	return(out);
}


// function that just swaps two items at random within a vector
// there used to be another version but it wasn't used so bye-bye
// the reason this doesn't edit in place is because in the genetic algorithm,
// unused columns in `pop` get overwritten. 
inline NumericVector swapcol2 (const NumericVector& x){
	NumericVector xout = clone(x);
	IntegerVector xpos = twoDifferentRandomInts(xout.size());
	std::swap(xout[xpos[0]], xout[xpos[1]]);
	return xout;
}

// like above, but may do more than one swap
// number of swaps is drawn at random from an vector of ints ns
// the reason this doesn't edit in place is because in the genetic algorithm,
// unused columns in `pop` get overwritten. 
inline NumericVector swapcoln (const NumericVector& x, const IntegerVector ns){
	NumericVector xout = clone(x);
	int nspos = floor(R::runif(0,1) * ns.size());
	for(int i=0; i<ns[nspos]; i++){
		IntegerVector xpos = twoDifferentRandomInts(x.size());
		std::swap(xout[xpos[0]], xout[xpos[1]]);
	}
	return(xout);
}

// function that randomly permutes a vector
inline NumericVector permcol (const NumericVector& x){
	int xsize = x.size();
	NumericVector xout(xsize);
	IntegerVector xpos = seq(0, xsize - 1);
	shuffle_ints_so(xpos); // seed argument ommitted, default is just to not set a seed
	for(int i=0; i<xsize; i++){
		xout[i] = x[xpos[i]];
	}
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
List rao_genetic_max (const NumericVector& p, const NumericVector& D, 
	const IntegerVector swap_freq,
	int term_cycles=10, int maxiters = 400, 
	int popsize_perm=150, int popsize_swap=150, 
	int keep=5, double prc = 0.001){
	// check to make sure p and D are correct lengths
	int ps = p.size();
	int p_expctd = (ps * (ps - 1)) / 2;
	if(D.size() != p_expctd){
		throw std::invalid_argument("p and D incompatible lengths.");
	}
	// calculate total popsize
	int popsize = popsize_perm + popsize_swap + 1; // +1 is for keeping original p
	// check that all integer args are reasonable
	if( (term_cycles < 1) | (maxiters < 1) | (popsize < 1) | (keep < 1) ){
		throw std::invalid_argument("One or more integer arguments is too low.");
	}
	// check that keep < popsize
	if(keep >= popsize){
		throw std::invalid_argument("keep cannot be >= popsize.");
	}

	// propagate initial matrix; each col is a p vector (a population of p vecs)
	NumericMatrix pop = NumericMatrix(ps, popsize);
	// do perm cols first
	for(int j=0; j<popsize_perm; j++){
		pop(_,j) = permcol(p); 
	}
	// now do swap cols
	for(int j=popsize_perm; j<(popsize_perm+popsize_swap); j++){
		pop(_,j) = swapcoln(p, swap_freq); 
	}
	// final column is initialized with p
	// this is done so that if p is already the best state, it is kept
	// and we don't have to hope we find it later
	pop(_,(popsize_perm+popsize_swap)) = p;

	// set up counters and outputs for generations
	double last_best_rao = 0.0;
	NumericVector iter_best_raos(maxiters, NA_REAL);
	int n_gens_no_improvement = 1;
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
		NumericVector pop_raos(popsize);
		for(int j=0; j<popsize; j++){
			pop_raos[j] = rao1sp(pop(_,j), D);
		}

		// find best n=keep 
		best_n_pos = which_top_n(pop_raos, keep);

		// record best rao for posterity
		iter_best_raos[iter] = max(pop_raos);
		
		// copy keeper columns into a new object so they don't get overwritten below
		int j = 0;
		for(int k = 0; k < keep; k++){       // k is the column of keepmat
			j = best_n_pos[k];               // j is the column of pop
			keepmat(_,k) = pop(_,j);  // put kth best column of pop into keepmat as column k
		}

		// start filling up pop with new columns
		// do keepers first
		for(int j=0; j<keep; j++){
			pop(_,j) = keepmat(_,j);
		}
		// now do swap cols
		IntegerVector col2do = seqrep(0, keep - 1, popsize); // vector of positions in keep range
		for(int j=keep; j<popsize; j++){
			pop(_,j) = swapcoln(keepmat(_,col2do[j]), swap_freq);
		}

		// check stop conditions and modify stopper
		if(double_equals(iter_best_raos[iter], last_best_rao, prc)){
			n_gens_no_improvement ++;
		}else{
			n_gens_no_improvement = 1;
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

	// warning if stopped because of maxiters
	if(iter >= maxiters){
		warning("Reached maxiters");
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

