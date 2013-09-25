#ifndef NP_UTILS_H
#define NP_UTILS_H

const int CHUNK_MULTI = 64;

#include <iostream>
#include <fstream>
#include <math.h>

#include "np_math.h"
#include "np_predictor.h"
#include <utils.h>
#include <bitset>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

extern int vlevel;
#define vcerr(verboselevel) if (verboselevel > vlevel) ; else std::cerr

/*
 * Top level function for doing MCMC estimation of parameters
 * infF:        file with forward read counts
 * infR:        file with forward read counts
 * outf:        file where results are stored
 * minD:        min allowed distance between forward and reverse reads to define a nucleosome center
 * maxD:        max allowed distance between forward and reverse reads to define a nucleosome center
 * sampleSize:  number of observations from the window sum data to base inference on
 * burnins:     number of burnin iterations
 * iterations:  number of iterations to base estimation on
 * verbose:     level of information printed to STDERR, 0: none, 3: alot
 */

int estimateParameters(ifstream* infF, ifstream* infR, ofstream* outf,
					   int minD, int maxD, int lambdaD,
					   double* pF, double* pR, double* lambdaF, double* lambdaR,
					   int sampleSize, int burnins, int iterations,
					   int seed, int minPos, int maxPos,
					   int minR, int minF,
					   int maxF, int maxR,
		       int truncValF, int truncValR);

int checkFileLengths(ifstream* infF, ifstream* infR,
					 int* minPos, int* maxPos,
					 int* minR, int* minF,
					 int* maxF, int* maxR);

int calculateTruncVal(ifstream* infF, ifstream* infR, ofstream* outf,
					  int* truncValF, int* truncValR, double truncLimit,
					  int minPos, int maxPos,
					  int minR, int minF,
					  int maxF, int maxR);

int truncData(storageType* data, int size, int truncVal);

int estimateDistance(ifstream* infF, ifstream* infR, ofstream* outf,
					 int* minD, int* maxD, int* lambdaD,
					 int minLag, int maxLag,
					 int minPos, int maxPos,
					 int minR, int minF,
					 int maxF, int maxR,
					 int truncValF, int truncValR,double* distDens,
					 double dens);

int windowSums(ifstream* infF, ifstream* infR, int* centers,
			   int sampleSize, int maxPos, int minPos,
			   int minF, int minR, int minD, int maxD,
			   int* maxValF, int* maxValR,
			   int thetaF1, int thetaF2, int thetaR1, int thetaR2,
			   int* F, int* R,
			   double* meanF, double* meanR,
			   double* varF, double* varR,
			   double* covFR,
			   int truncValF, int truncValR);
			   
int gibbs(int* F, int* R, bool* ZF, bool* ZR,
		  double* pF, double* pR,
		  double* lambdaF, double* lambdaR,
		  double* deltaF, double* deltaR,
		  double alphaF, double alphaR, double betaF, double betaR,
		  int maxValF, int maxValR,
		  int sampleSize, int iterations, int burnins, int seed,
		  ofstream* outf);

int readParameters(ifstream* ipff, double* pF, double* pR, double* lambdaF, double* lambdaR,int* lambdaD,int* minD,int* maxD,
				   int* estLambdaD,int* estMinD,int* estMaxD,vector<string>* corrv);

int writeParameter(ofstream* opff,string,string);
		   


int writeParameters(ofstream* opff, double* pF, double* pR, double* lambdaF, double* lambdaR, int lambdaD,int maxD,int minD,
		    int estLambdaD,int estMaxD, int extMinD,double* corr);


int centerPositionOdds(ifstream* infF, ifstream* infR,
		       ofstream* outf,
		       double* pF, double* pR, double* lambdaF, double* lambdaR,
		       int minD, int maxD,
		       int* maxValF, int* maxValR,
		       int startCrd, int endCrd,
		       int minPos, int maxPos,
		       int minR, int minF,
		       int maxF, int maxR,
		       int truncValF, int truncValR);
  

int predictBoundaries(ifstream* infF, ifstream* infR,
					  ifstream* infO, ofstream* outf,
					  double* pF, double* pR, double* lambdaF, double* lambdaR,
					  int lambdaD, int minD, int maxD,
					  int minLag, int maxLag,
					  int maxValF, int maxValR,
					  int startCrd,int endCrd,
					  int minPos, int maxPos,
					  int minR, int minF,
					  int maxF, int maxR,
					  int truncValF, int truncValR,
					  double* distDens);
					   

#endif //NP_UTILS_H
