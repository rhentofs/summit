#include "variableclusterpredictor.h"

VariableClusterPredictor::VariableClusterPredictor(ofstream *outf, int lambdaD, int ws, int minD, int maxD,
												   ifstream *infF, ifstream *infR, int minPos, int maxPos, int intWs) :
	ClusterPredictor(outf, lambdaD, ws, minD, maxD), minPos(minPos), maxPos(maxPos), intWs(intWs)
{
	this->infF = infF;
	this->infR = infR;
}

VariableClusterPredictor::~VariableClusterPredictor()
{
  int debug = 0;
  // Check if all regions have been flushed
  
  if(debug) cerr << "In destructor" << endl;
  if(debug) cerr << "distMap: " << distMap.size() << endl;
  if(debug) cerr << "scoreMap: " << scoreMap.size() << endl;
  if(debug) cerr << "ends: " << ends.size() << endl;
  
  if (!scoreMap.empty())
    this->flush();
  
  if(debug) cerr << "In destructor" << endl;
  if(debug) cerr << "distMap: " << distMap.size() << endl;
  if(debug) cerr << "scoreMap: " << scoreMap.size() << endl;
  if(debug) cerr << "ends: " << ends.size() << endl;
}

double VariableClusterPredictor::distance(Info *first, Info *second)
{
	// Distance is 0 - the overlap between the regions giving values in interval [-1.0,0.0]
	double olen = (min(first->get(1),second->get(1)) - max(first->get(0), second->get(0)) + 1);
	double len = (max(first->get(1),second->get(1)) - min(first->get(0), second->get(0)) + 1);
	if (len <= 0.0 || olen <= 0.0)
		return(0.0);
	return(0.0 - (olen / len));
}

bool VariableClusterPredictor::member(Info *info, double minStart,
									  double maxStart, double minEnd, double maxEnd)
{
	// Member if prediction overlap with cluster region [maxStart,minEnd]
	//double olen = (min(info->get(1),minEnd) - max(info->get(0), maxStart) + 1);
	//return(olen > 0.0);

	// Member if old cluster center is contained within new cluster region [maxStart+ws-intWs,minEnd-ws+intWs]
	double center = 0.5 * (maxStart + minEnd);
	return((center < min(minEnd,info->get(1))-ws+intWs) && (center > max(maxStart,info->get(0))+ws-intWs));
}

void VariableClusterPredictor::predict(set<double> *predictions, bool split, int verbose)
{
	map<double, pair<Info*, multimap<double, double, distComp>*>, distComp>::iterator distIt;
	set<double>::iterator predIt;
	int nPreds;
	double minStart, maxStart, minEnd, maxEnd;

	nPreds = (int)predictions->size();
	predIt = predictions->begin();
	
	distIt = distMap.find((*predIt));
	minStart = maxStart = (*distIt).second.first->get(0);
	minEnd = maxEnd = (*distIt).second.first->get(1);

	removeEntry((*predIt),verbose);
	predictions->erase(predIt);

	predIt++;
	
	// Go through predictions
	if (verbose)
		cerr << "Merging predictions" << endl;

	while (predIt != predictions->end())
	{
	  distIt = distMap.find((*predIt));

	  minStart = min(minStart, (*distIt).second.first->get(0));
	  maxStart = max(maxStart, (*distIt).second.first->get(0));
	  minEnd = min(minEnd, (*distIt).second.first->get(1));
	  maxEnd = max(maxEnd, (*distIt).second.first->get(1));
	  
	  removeEntry((*predIt),verbose);
	  predictions->erase(predIt);

	  predIt++;
	}
	
	int size = (int)maxStart-(int)minStart+ws-intWs;

	//cerr << "BOUNDARIES: " << minStart << " " << maxStart << " " << minEnd << " " << maxEnd << " " << size << " " << nPreds << endl;
	//cerr << "BOUNDARY start: " << (int)minStart << ", " << (int)minStart+size << endl;
	//cerr << "BOUNDARY end: " << (int)minEnd - ws + intWs << ", " <<  (int)minEnd - ws + intWs + size << endl;
	//cerr << "minPos: " << minPos << endl;

	int h = 5;
	int lowLag = h, highLag = h;
	lowLag = (minStart - minPos - lowLag < 0 ? minStart - minPos : lowLag);
	highLag = (maxPos - maxEnd - highLag < 0 ? maxPos - maxEnd : highLag);

	size = size + lowLag + highLag;
	
	unsigned short *dataF = new unsigned short[size];
	unsigned short *dataR = new unsigned short[size];
	
	infF->seekg(sizeof(unsigned short) * ((int)minStart - minPos - lowLag + 1) + sizeof(int), ios::beg);
	infF->read((char*)dataF, sizeof(unsigned short) * size);
	
	infR->seekg(sizeof(unsigned short) * ((int)minEnd - ws + intWs - minPos - lowLag + 1) + sizeof(int), ios::beg);
	infR->read((char*)dataR, sizeof(unsigned short) * size);

	// Parzen's window kernel density estimation using a gaussian kernel and bandwith h (5)
	double* Fdens = new double[size];
	double* Rdens = new double[size];
	double Fsum;
	double Rsum;
	int Fcnt = 0;
	int Rcnt = 0;

	double sumF = 0;
	double sumR = 0;
	int cntF = 0;
	int cntR = 0;
	
	for (int i = 0; i < size; i++)
	{
		Fsum = 0.0;
		Rsum = 0.0;
		Fcnt = 0.0;
		Rcnt = 0.0;

		if ((i >= lowLag) && (i < size - highLag))
		{
			sumF += (double)dataF[i] * ((double)(i-lowLag)/minStart + 1.0);
			sumR += (double)dataR[i] * ((double)(i-lowLag)/(double)((int)minEnd-ws+intWs) + 1.0);
			cntF += (int)dataF[i];
			cntR += (int)dataR[i];
		}
		
		for (int j = max(0,i-h); j < min(i+h,size); j++)
		{
			if ((int)dataF[j] > 0)
			{
				Fsum += (double)dataF[j] * gsl_ran_gaussian_pdf((double)(i-j), (double)h);
				Fcnt++;
			}

			if ((int)dataR[j] > 0)
			{
				Rsum += (double)dataR[j] * gsl_ran_gaussian_pdf((double)(i-j), (double)h);
				Rcnt++;
			}
		}
		Fdens[i] = (Fcnt > 0 ? Fsum / (double)(Fcnt * h) : 0.0);
		Rdens[i] = (Rcnt > 0 ? Rsum / (double)(Rcnt * h) : 0.0);
	}

	double prob = 1.0;
	
	double aveF = minStart * sumF / (double)cntF;
	double aveR = (double)((int)minEnd-ws+intWs) * sumR / (double)cntR;

	double sumSdS, sumSdE;
	sumSdS = 0.0;
	sumSdE = 0.0;
	
	for (int i = lowLag; i < size - highLag; i++)
	{
		sumSdS += (double)dataF[i] * pow((double)(i+(int)minStart) - aveF, 2);
		sumSdE += (double)dataR[i] * pow((double)(i+(int)minEnd-ws+intWs) - aveR, 2);
	}

	//cerr << aveF << " " << aveR << " " << sumSdS << " " << sumSdE << " " << cntF << " " << cntR << endl;
	
	double maxHigh = GSL_NAN;
	double maxLow = GSL_NAN;
	double best = GSL_NEGINF;
	double val;
	int nobests = 0;
	int dist;
	// Find the most probable (first occurrence) low and high pair
	bool noF = true;
	bool noR = true;
	for (int i = lowLag; i < size-highLag; i++)
	{
		if (Fdens[i] > 0.0)
		{
			noF = false;
			for (int j = lowLag; j < size-highLag; j++)
			{
				if (Rdens[j] > 0.0)
				{
					noR = false;
					dist = (j-lowLag+(int)minEnd-ws+intWs)-(i-lowLag+(int)minStart);
					
					//cerr << (i-lowLag+(int)minStart) << ", " << (j-lowLag+(int)minEnd-ws+intWs);
					
					if (dist < minD)
					{
						//cerr << " lower" << endl;
						continue;
					}
					else if (dist > maxD)
					{
						//cerr << " higher" << endl;
						break;
					}
								
					// Keep track of 'best' pair for maxHigh and maxLow
					val = Fdens[i] * Rdens[j] * gsl_ran_poisson_pdf(dist, lambdaD);

					//cerr << ", " << Fdens[i] << ", " << Rdens[j] << ", " << dataF[i] << ", " << dataR[j] << ", "<< val;
					
					if (val > 0.0)
					{
						// Best candidate so far?
						if (val > best)
						{
							best = val;
							maxLow = (double)(i-lowLag+(int)minStart);
							maxHigh = (double)(j-lowLag+(int)minEnd-ws+intWs);
							nobests = 1;

							//cerr << " best";
						}
						// Keep track of the number of equal best candidates
						else if (val == best)
						{
							nobests++;
							//cerr << " best";
						}
						//cerr << endl;
					}
				}
			}
		}
	}

	//if (nobests > 1)
	//	cerr << nobests << " equal best candidates, using the first one only" << endl;
	
	double predS = maxLow;
	double predE = maxHigh;
	
	delete[] dataF;
	delete[] dataR;
	delete[] Fdens;
	delete[] Rdens;

	// Write prediction
	if (verbose)
		cerr << "Writing prediction" << endl;

	//if (!gsl_isnan(predS) && !gsl_isnan(predE))
	writePrediction(predS, predE, (split ? "D": "C"),
					prob,
					sqrt(sumSdS / (double)cntF),
					sqrt(sumSdE / (double)cntR),
					nPreds);
}
