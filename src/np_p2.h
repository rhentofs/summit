#ifndef NP_PREDICTOR_H
#define NP_PREDICTOR_H

#include <map>
#include <set>
#include <fstream>
#include <gsl/gsl_math.h>
#include <iostream>
#include "info.h"
#include <gsl/gsl_randist.h>

using namespace std;

class NPPredictor
{

public:
	NPPredictor(ofstream *outf, int ws, int minD, int maxD, int minLag, int maxLag,
				ifstream *infF, ifstream *infR, int minPos, int maxPos, int intWs, double* distDens,
				double* pF, double* pR, double* lambdaF, double* lambdaR);
	~NPPredictor();

	int addPrediction(double predStart, double predEnd,
					  double postIO, int pos, int verbose); //
	
private:
	void flushHigh(int verbose = 0); //
	
	void flush(int verbose = 0); //
	
	void forcedFlush(int predStart,int verbose = 0); //
	
	void writePrediction(double start, double end, const char *type, double *score, double startSd, double endSd, int support); //
	
	pair<double,Info*> createMapPair(double key, double first, double second, double third, double fourth=0.0); //

	void predict(set<double> *predictions, bool split, int verbose); //

	double distance(Info *first, Info *second); //

	bool member(Info *info, double minStart, double maxStart, double minEnd, double maxEnd); //

	void populateDistMap(pair<double,Info*> p); //

	void printDistMap(); //

	void removeEntry(double id, int verbose=0); //
	
	ofstream *outf;
	int ws;
	int minD;
	int maxD;
	int minLag;
	int maxLag;
	ifstream *infF, *infR;
	int minPos;
	int maxPos;
	int intWs;
	int scnt;
	int ecnt;	
	int id;
	double* distDens;
	double* pF;
	double* pR;
	double* lambdaF;
	double* lambdaR;
	
	struct comp {
		bool operator() (const double& lhs, const double& rhs) const
			{return lhs<rhs;}
	};
	
	// end -> <start,score,id>
	multimap<double,Info*,comp> ends;
	
	struct distComp {
		bool operator() (const double& lhs, const double& rhs) const
			{return lhs<rhs;}
	};

	struct scoreComp {
		bool operator() (const double& lhs, const double& rhs) const
			{return lhs>rhs;}
	};

	//id -> <Info*, <dist -> id>>
	map<double, pair<Info*, multimap<double, double, distComp>*>, distComp> distMap;

	//score -> id
	multimap<double, double, scoreComp> scoreMap;
};

#endif
