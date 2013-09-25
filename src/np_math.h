#ifndef NP_MATH_H
#define NP_MATH_H

#include <gsl/gsl_rng.h>
#include <math.h>
#include <vector>
#include <algorithm>

using namespace std;

/*
 * Round v to nearest int
 */
int np_round(double v);

/*
 * Sample len ints from {from,from+by,...,to-by,to}
 * store in int vector vec (of size len)
 */
void np_sample(int seed, int from, int to, int by, int len, int* vec);

#endif //NP_MATH_H
