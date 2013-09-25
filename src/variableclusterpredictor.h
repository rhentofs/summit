#ifndef VARIABLECLUSTERPREDICTOR_H
#define VARIABLECLUSTERPREDICTOR_H

#include "clusterpredictor.h"
#include <fstream>
#include <gsl/gsl_randist.h>

using namespace std;

class VariableClusterPredictor : public ClusterPredictor
{
public:
	VariableClusterPredictor(ofstream *outf, int lambdaD, int ws, int minD, int maxD,
							 ifstream *infF, ifstream *infR, int minPos, int maxPos, int intWs);
	~VariableClusterPredictor();
private:
	ifstream *infF, *infR;
	int minPos;
	int maxPos;
	int intWs;
	
	double distance(Info *first, Info *second);
	bool member(Info *info, double minStart, double maxStart, double minEnd, double maxEnd);
	void predict(set<double> *predictions, bool split, int verbose);
};
	
#endif
