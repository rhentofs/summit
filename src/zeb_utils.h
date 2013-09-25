#ifndef ZEB_UTILS_H
#define ZEB_UTILS_H

const int CHUNK_MULTI = 64;

#include <iostream>
#include <fstream>
#include <math.h>
#include <dirent.h>

#include "np_math.h"
#include <utils.h>
#include <bitset>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

extern int vlevel;
#define vcerr(verboselevel) if (verboselevel > vlevel) ; else std::cerr

int countFiles(string infPath);
int openFiles(string infPath, ifstream* iff, ifstream* ifr, string* fileNames);
int openFilesNames(string infPath, ifstream* iff, ifstream* ifr, string* fileNames, int cnt);

int estimateDistance(ifstream* infF, ifstream* infR, ofstream* outf,
					 int* minD, int* maxD, int* lambdaD,
					 int minLag, int maxLag,
					 int minPos, int maxPos,
					 int minR, int minF,
					 int maxF, int maxR,
					 int truncValF, int truncValR, double* distDens,
					 double dens);

int estimateParameters(ifstream* infF, ifstream* infR, ofstream* outf,
					   int minD, int maxD, int lambdaD,
					   double* pF, double* pR, double* lambdaF, double* lambdaR,
					   int sampleSize, int burnins, int iterations,
					   int seed, int* minPos, int* maxPos,
					   int* minR, int* minF,
					   int* maxF, int* maxR,
					   int truncValF, int truncValR, int fCnt);

int checkFileLengths(ifstream* infF, ifstream* infR,
					 int* minPos, int* maxPos,
					 int* minR, int* minF,
					 int* maxF, int* maxR,
					 int fCnt);
int checkFileLengths(ifstream* infF1, ifstream* infF2,
					 ifstream* infR1, ifstream* infR2,
					 int* minPos, int* maxPos,
					 int* minR1, int* minR2,
					 int* minF1, int* minF2,
					 int* maxF1, int* maxF2,
					 int* maxR1, int* maxR2,
					 int fCnt);

int calculateTruncVal(ifstream* infF, ifstream* infR, ofstream* outf,
					  int* truncValF, int* truncValR, double truncLimit,
					  int* minPos, int* maxPos,
					  int* minR, int* minF,
					  int* maxF, int* maxR,
					  int fCnt);

int truncData(storageType* data, int size, int truncVal);

int windowSums(ifstream* infF, ifstream* infR, int* centers,
			   int* sampleSizes, int* maxPos, int* minPos,
			   int* minF, int* minR, int minD, int maxD,
			   int* maxValF, int* maxValR,
			   int thetaF1, int thetaF2, int thetaR1, int thetaR2,
			   int* F, int* R,
			   double* meanF, double* meanR,
			   double* varF, double* varR,
			   double* covFR,
			   int truncValF, int truncValR, int fCnt);
			   
int gibbs(int* F, int* R, bool* ZF, bool* ZR,
		  double* pF, double* pR,
		  double* lambdaF, double* lambdaR,
		  double* deltaF, double* deltaR,
		  double alphaF, double alphaR, double betaF, double betaR,
		  int maxValF, int maxValR,
		  int sampleSize, int iterations, int burnins, int seed,
		  ofstream* outf);

int readParameters(ifstream* ipff, double* pF, double* pR, double* lambdaF, double* lambdaR,int* lambdaD,
				   int* truncValF, int* truncValR,
				   int* minD,int* maxD,
				   int* estLambdaD,int* estMinD,int* estMaxD,vector<string>* corrv);

int writeParameter(ofstream* opff,string,string);
		   
int writeParameters(ofstream* opff, double* pF, double* pR, double* lambdaF, double* lambdaR, int lambdaD,
					int truncValF, int truncValR,
					int maxD,int minD,
					int estLambdaD,int estMaxD, int extMinD,double* corr);

int changeOdds(ifstream* infF1, ifstream* infF2, ifstream* infR1, ifstream* infR2,
			   ofstream* outo, ofstream* outos, ofstream* outoe, ofstream* outd, ofstream* outs,
			   bool bodds,
			   double* pF1, double* pF2, double* pR1, double* pR2,
			   double* lambdaF1, double* lambdaF2, double* lambdaR1, double* lambdaR2,
			   int minD, int maxD,
			   int* minPos, int* maxPos,
			   int* minR1, int* minR2, int* minF1, int* minF2,
			   int* maxF1, int* maxF2, int* maxR1, int* maxR2,
			   int truncValF1, int truncValF2, int truncValR1, int truncValR2,
			   int fCnt, string* fileNames);

int centerPositionOdds(ifstream* infF, ifstream* infR,
					   ofstream* outo,
					   ofstream* outos,
					   ofstream* outoe,
					   ofstream* outf,
					   ofstream* outs,
					   bool bodds,
					   double* pF, double* pR, double* lambdaF, double* lambdaR,
					   int minD, int maxD,
					   int* maxValF, int* maxValR,
					   int* minPos, int* maxPos,
					   int* minR, int* minF,
					   int* maxF, int* maxR,
					   int truncValF, int truncValR,
					   int fCnt, string* fileNames);

int centerPositionOddsNeg(ifstream* infF, ifstream* infR,
						  ofstream* outo,
						  ofstream* outos,
						  ofstream* outoe,
						  ofstream* outf,
						  ofstream* outs,
						  bool bodds,
						  double* pF, double* pR, double* lambdaF, double* lambdaR,
						  int minD, int maxD,
						  int* maxValF, int* maxValR,
						  int* minPos, int* maxPos,
						  int* minR, int* minF,
						  int* maxF, int* maxR,
						  int truncValF, int truncValR,
						  int fCnt, string* fileNames);

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
