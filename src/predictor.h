#ifndef PREDICTOR_H
#define PREDICTOR_H

#include <map>
#include <set>
#include <fstream>
#include <gsl/gsl_math.h>
#include <iostream>
#include "info.h"

using namespace std;

class Predictor
{

public:
	Predictor(ofstream *outf, int lambdaD, int ws, int minD, int maxD);
	virtual ~Predictor();

	virtual int addPrediction(double predStart, double predEnd, double predStartSd,
							  double predEndSd, double postIO, int pos, int verbose);
	
protected:
	ofstream *outf;
	int lambdaD;
	int ws;
	int minD;
	int maxD;

	struct comp {
		bool operator() (const double& lhs, const double& rhs) const
			{return lhs<rhs;}
	};
	
	// start -> <end,score,startsd,endsd,id>
	multimap<double,Info*,comp> starts;
	// end -> <start,score,startsd,endsd,id>
	multimap<double,Info*,comp> ends;

	virtual void flush(int verbose = 0) = 0;
	virtual void forcedFlush(int predStart,int verbose = 0) = 0;

	
	virtual void writePrediction(double start, double end, const char *type, double score, double startSd, double endSd, int support);
	
	pair<double,Info*> createMapPair(double key, double first, double second, double third, double fourth, double fifth, double sixth=0.0);

	virtual void predict(set<double> *predictions, bool split, int verbose);
	
private:
	int scnt;
	int ecnt;
};

#endif
