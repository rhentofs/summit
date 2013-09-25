#include "np_math.h"

/*
 * Round v to nearest int
 */
int np_round(double v)
{
	return((int)(v+0.5));
}

/*
 * Sample len ints from {from,from+by,...,to-by,to}
 * store in int vector vec (of size len)
 */
void np_sample(int seed, int from, int to, int by, int len, int* vec)
{
	int steps = (to - from) / by;

	gsl_rng* r;
	
	r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, seed);
    
	// Sample the data
	vector<int> tmp(len, 0);
	for (int i = 0; i < len; i++)
		tmp[i] = (int)gsl_rng_uniform_int(r, steps) * by + from;
	gsl_rng_free(r);
	
	// Sort data in increasing order
	sort(tmp.begin(), tmp.end());
	for (int i = 0; i < len; i++)
		vec[i] = tmp[i];
}
		
