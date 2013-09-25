#ifndef CLUSTERPREDICTOR_H
#define CLUSTERPREDICTOR_H

#include <iostream>
#include <map>
#include <fstream>
#include <gsl/gsl_math.h>
#include "predictor.h"

using namespace std;

class ClusterPredictor : public Predictor
{
public:
	ClusterPredictor(ofstream *outf, int lambdaD, int ws, int minD, int maxD);
	~ClusterPredictor();

	int addPrediction(double predStart, double predEnd, double predStartSd,
					  double predEndSd, double postIO, int pos, int verbose);
	
protected:

	virtual void predict(set<double> *predictions, bool split, int verbose);

	virtual void flush(int verbose=0);
	
	virtual void forcedFlush(int predStart,int verbose=0);
	
	virtual double distance(Info *first, Info *second) = 0;

	virtual bool member(Info *info, double minStart, double maxStart, double minEnd, double maxEnd) = 0;
	
	void removeEntry(double id, int verbose=0);

	void printDistMap();
	
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
	
private:
	int id;

	void populateDistMap(pair<double,Info*> p);
};
	
#endif
