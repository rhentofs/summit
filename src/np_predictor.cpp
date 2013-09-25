#include "np_predictor.h"
#include <vector>

NPPredictor::NPPredictor(ofstream *outf, int ws,
						 int minD, int maxD, int minLag, int maxLag,
						 ifstream *infF, ifstream *infR,
						 int minPos, int maxPos, int intWs, double* distDens,
						 double* pF, double* pR, double* lambdaF, double* lambdaR) :
	ws(ws), minD(minD), maxD(maxD), minLag(minLag),
	maxLag(maxLag), minPos(minPos), maxPos(maxPos), intWs(intWs)
{
	this->infF = infF;
	this->infR = infR;
	this->outf = outf;
	this->distDens = distDens;
	this->id = 0;
	this->pF = pF;
	this->lambdaF = lambdaF;
	this->pR = pR;
	this->lambdaR = lambdaR;	
}

NPPredictor::~NPPredictor()
{
	int debug = 0;
  
	// Check if all regions have been flushed
	if (!scoreMap.empty())
		this->flush();

	// This is not needed but what the heck...
	if (!ends.empty())
	{
		multimap<double,Info*,comp>::iterator tmpIt;
	  
		for (tmpIt = ends.begin(); tmpIt != ends.end(); tmpIt++)
			delete (*tmpIt).second;

		ends.clear();
	}
}

pair<double,Info*> NPPredictor::createMapPair(double key, double first, double second, double third, double fourth)
{
	pair<double,Info*> p;

	p.first = key;
	p.second = new Info(first, second, third, fourth);
	
	return(p);
}

void NPPredictor::populateDistMap(pair<double,Info*> p)
{
	// Do we need to have a pointer to the pair?
	pair<double,pair<Info*,multimap<double, double, distComp>* > > entry;
	
	entry.first = p.first;
	entry.second.first = p.second;
	entry.second.second = new multimap<double, double, distComp>;
	entry.second.second->clear();

	// Non-empty distmap?
	if (distMap.empty())
	{
		distMap.insert(entry);
		
		//printDistMap();
	}
	else
	{
		// Populate new entry with distances to the existing entries from the new entry
		// and add distances to the new entry from the existing entries
		map<double, pair<Info*, multimap<double, double, distComp>*>, distComp>::iterator distIt;
		double dist;
		for (distIt = distMap.begin(); distIt != distMap.end(); distIt++)
		{
			dist = distance(entry.second.first, (*distIt).second.first);
			entry.second.second->insert(pair<double,double>(dist,(*distIt).first));
			(*distIt).second.second->insert(pair<double,double>(dist,p.first));
		}
		distMap.insert(entry);
		//printDistMap();
	}
}

int NPPredictor::addPrediction(double predStart, double predEnd,
							   double postIO, int pos, int verbose)
{

	Info *info = new Info(predStart,predEnd,postIO,double(id));
	// Only add if proper prediction (non-nan)  
	if (!gsl_isnan(predStart) && !gsl_isnan(predEnd) && !gsl_isnan(postIO))
    { 
		//cerr<<"ClusterPredictor::aP:\t"<<postIO<<"\t"<<pos<<endl;
		if (verbose)
		{
			cerr << "Adding prediction" << endl;
			cerr << pos << " " <<predStart << " " << predEnd << endl;
			cerr << "Before: distMap.size() == " << distMap.size() << " ends.size() == " << ends.size() << " scoreMap.size() == " << scoreMap.size() << endl;
		}
      
		// Empty multimap: just add prediction
		if (scoreMap.empty())
		{
			scoreMap.insert(pair<double,double>(postIO,(double)id));
			ends.insert(createMapPair(predEnd, predStart, postIO, (double)id));
			populateDistMap(createMapPair((double)id, predStart, predEnd, postIO, (double)id));
			//printDistMap();
			id++;
		}
		// Non-empty multimap
		else
		{
			map<double, pair<Info*, multimap<double, double, distComp>*>, distComp>::iterator entry;
			
			double sid = (*scoreMap.begin()).second;
			double score = (*scoreMap.begin()).first;

			entry = distMap.find(sid);
			double start = (*entry).second.first->get(0);
			double end = (*entry).second.first->get(1);

			// New prediction sufficiently far away from the one with highest score?
			bool safe = (score > postIO) && !member(info,start,start,end,end);
			while (safe)
			{
				if (verbose)
					cerr << "pos - end >= min(ws,maxD): " << pos << " - " << end << ">= " << min(ws,maxD) << endl;
				this->flushHigh(verbose);
	      
				safe = false;
				if (!scoreMap.empty())
				{
					sid = (*scoreMap.begin()).second;
					score = (*scoreMap.begin()).first;
					entry = distMap.find(sid);
					start = (*entry).second.first->get(0);
					end = (*entry).second.first->get(1);
					safe = (score > postIO) && !member(info,start,start,end,end);
				}
			}
	  
			// 			  multimap<double,Info*,comp>::reverse_iterator end = ends.rbegin();
	  
			// 			// New prediction sufficiently far away from the last old one? 
			// 			if (pos - (*end).first >= min(ws,maxD))
			// 			{
			// 				if (verbose)
			// 					cerr << "pos - (*end).first >= min(ws,maxD): " << pos << " - " << (*end).first << ">= " << min(ws,maxD) << endl;
			// 				this->flush(verbose);
			// 			}
			// 			else
			// 			{
			// 				int max_dist = 2*maxD; // completly as hoc!
			// 				multimap<double,Info*,comp>::iterator start = ends.begin();	  
			// 				// New prediction sufficiently far away from the many of the old ones? Shaky...
			// 				if (pos - (*start).first >= max_dist) // distance 'pos' to 'end'-coordinate in first prediction.    
			// 				{
			// 					if (verbose)
			// 					{
			// 						cerr << "pos - (*start).first >= max_dist: " << pos << " - " << (int)(*start).first;
			// 						cerr << "( = "<< pos - (int)(*start).first<< ") >= " << max_dist << "\t"<< ends.size()<<endl;
			// 					}
			// 					this->forcedFlush(predStart,verbose);
			// 				}
			// 			}
	  
			scoreMap.insert(pair<double,double>(postIO,(double)id));
			ends.insert(createMapPair(predEnd, predStart, postIO, (double)id));
			populateDistMap(createMapPair((double)id, predStart, predEnd, postIO, (double)id));
			id++;
		}
      
		if (verbose)
			cerr << "After: distMap.size() == " << distMap.size() << " ends.size() == " << ends.size() << " scoreMap.size() == " << scoreMap.size() << endl;
    }
  
	delete info;
  
	return(0);
}

void NPPredictor::flushHigh(int verbose)
{
	double id;
	Info *info;
	multimap<double, double, distComp>::iterator distIt;
	map<double, pair<Info*, multimap<double, double, distComp>*>, distComp>::iterator entry;
	multimap<double, double, distComp> *dist;
	set<double> preds;
	
	double minStart, maxStart, minEnd, maxEnd;
	bool split;

	//printDistMap();
	
	// Flush highest scoring prediction
	if (verbose)
		cerr << "Flushing predictions" << endl;

	id = (*scoreMap.begin()).second;
	entry = distMap.find(id);
	if (verbose)
		cerr << id << (entry != distMap.end() ? "found" : "not found") << endl;
	minStart = maxStart = (*entry).second.first->get(0);
	minEnd = maxEnd = (*entry).second.first->get(1);
	
	preds.insert(id);
	split = false;
	dist = (*entry).second.second;
	if (verbose)
		cerr << "Looping" << endl;
	for (distIt = dist->begin(); distIt != dist->end(); distIt++)
	{
		// Valid member to current cluster?
		info = (*distMap.find((*distIt).second)).second.first;
		if (member(info, minStart, maxStart, minEnd, maxEnd))
		{
			// Update min and max coordinates
			minStart = min(minStart, info->get(0));
			maxStart = max(maxStart, info->get(0));
			minEnd = min(minEnd, info->get(1));
			maxEnd = max(maxEnd, info->get(1));
			
			preds.insert(info->get(3));
		}
		else
		{
			split = (info->get(0) < maxEnd) && (info->get(1) > minStart);
			break;
		}
	}
	if (verbose)
		cerr << "Making prediction from predictions (" << preds.size() << ", " << scoreMap.size() << ")" << endl;
	
	// Make prediction from predictions
	unsigned int pSize = preds.size();
	unsigned int sSize = scoreMap.size();
	unsigned int dSize = distMap.size();

	predict(&preds, split, verbose);

	if (sSize - pSize != scoreMap.size())
	{
	    cerr << "scoreMap: " << sSize -pSize << " " << scoreMap.size() << endl;
	    cin.get();
	}
	if (dSize -pSize != distMap.size())
	{
	    cerr << "distMap: " << dSize -pSize << " " << distMap.size() << endl;
	    cin.get();
	}
	preds.clear();
}

double NPPredictor::distance(Info *first, Info *second)
{
	// Distance is 0 - the overlap between the regions giving values in interval [-1.0,0.0]
	double olen = (min(first->get(1),second->get(1)) - max(first->get(0), second->get(0)) + 1);
	double len = (max(first->get(1),second->get(1)) - min(first->get(0), second->get(0)) + 1);
	if (len <= 0.0 || olen <= 0.0)
		return(0.0);
	return(0.0 - (olen / len));
}

bool NPPredictor::member(Info *info, double minStart,
						 double maxStart, double minEnd, double maxEnd)
{
	//Member if center distance between the cluster and the new prediction is less than maxD
	double center = 0.5 * (maxStart + minEnd);
	double infoCenter = 0.5 * (info->get(0) + info->get(1));

	return(abs((int)center - (int)infoCenter) <= maxD);
}

void NPPredictor::predict(set<double> *predictions, bool split, int verbose)
{
	map<double, pair<Info*, multimap<double, double, distComp>*>, distComp>::iterator distIt;
	set<double>::iterator predIt;
	int nPreds;
	double minStart, maxStart, minEnd, maxEnd;
	double score[3];
	vector<int> startD;
	vector<int> endD;
	vector<double> scoreD;

	nPreds = (int)predictions->size();
	predIt = predictions->begin();

	distIt = distMap.find((*predIt));
	minStart = maxStart = (*distIt).second.first->get(0);
	minEnd = maxEnd = (*distIt).second.first->get(1);
	// DEBUG: keep track of all starts & ends in this prediction
	startD.push_back(minStart);
	endD.push_back(minEnd);

	score[0] = (*distIt).second.first->get(2);
	score[1] = (*distIt).second.first->get(2);
	score[2] = (*distIt).second.first->get(2);
	
	scoreD.push_back(score[0]);

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

		// DEBUG: keep track of all starts & ends in this prediction
		startD.push_back((*distIt).second.first->get(0));
		endD.push_back((*distIt).second.first->get(1));
		scoreD.push_back((*distIt).second.first->get(2));
	  
		score[0] = min(score[0], (*distIt).second.first->get(2));
		score[1] += (*distIt).second.first->get(2);
		score[2] = max(score[2], (*distIt).second.first->get(2));
	  
		removeEntry((*predIt),verbose);
		predictions->erase(predIt);

		predIt++;
	}
	score[1] /= (double)nPreds;
	
	int sizeF = (int)maxStart-(int)minStart+1+ws-intWs;
	int sizeR = (int)maxEnd-(int)minEnd+1+ws-intWs;

	int h = 5;
	int lowLagF = h, lowLagR = h, highLagF = h, highLagR = h;
	lowLagF = ((int)minStart - minPos - lowLagF < 0 ? (int)minStart - minPos : lowLagF);
	highLagF = (maxPos - (int)maxStart - highLagF < 0 ? maxPos - (int)maxStart : highLagF);
	lowLagR = ((int)minEnd - minPos - lowLagR < 0 ? (int)minEnd - minPos : lowLagR);
	highLagR = (maxPos - (int)maxEnd - highLagR < 0 ? maxPos - (int)maxEnd : highLagR);

	sizeF = sizeF + lowLagR + highLagF;
	sizeR = sizeR + lowLagR + highLagR;
	
	unsigned short *dataF = new unsigned short[sizeF];
	unsigned short *dataR = new unsigned short[sizeR];
	
	infF->seekg(sizeof(unsigned short) * ((int)minStart - minPos - lowLagF + 1) + sizeof(int), ios::beg);
	infF->read((char*)dataF, sizeof(unsigned short) * sizeF);
	
	infR->seekg(sizeof(unsigned short) * ((int)minEnd - ws + intWs - minPos - lowLagR + 1) + sizeof(int), ios::beg);
	infR->read((char*)dataR, sizeof(unsigned short) * sizeR);

	// Parzen's window kernel density estimation using a gaussian kernel and bandwith h (5)
	double* Fdens = new double[sizeF];
	double* Rdens = new double[sizeR];
	double Fsum;
	double Rsum;
	int Fcnt = 0;
	int Rcnt = 0;

	double sumF = 0;
	double sumR = 0;
	int cntF = 0;
	int cntR = 0;
	
	for (int i = 0; i < sizeF; i++)
	{
		Fsum = 0.0;
		Fcnt = 0.0;

		if ((i >= lowLagF) && (i < sizeF - highLagF))
		{
			sumF += (double)dataF[i] * ((double)(i-lowLagF)/minStart + 1.0);
			cntF += (int)dataF[i];
		}
		
		for (int j = max(0,i-h); j < min(i+h,sizeF); j++)
		{
			if ((int)dataF[j] > 0)
			{
				Fsum += (double)dataF[j] * gsl_ran_gaussian_pdf((double)(i-j), (double)h);
				Fcnt++;
			}
		}
		Fdens[i] = (Fcnt > 0 ? Fsum / (double)(Fcnt * h) : 0.0);
	}

	for (int i = 0; i < sizeR; i++)
	{
		Rsum = 0.0;
		Rcnt = 0.0;

		if ((i >= lowLagR) && (i < sizeR - highLagR))
		{
			sumR += (double)dataR[i] * ((double)(i-lowLagR)/(double)((int)minEnd-ws+intWs) + 1.0);
			cntR += (int)dataR[i];
		}
		
		for (int j = max(0,i-h); j < min(i+h,sizeR); j++)
		{
			if ((int)dataR[j] > 0)
			{
				Rsum += (double)dataR[j] * gsl_ran_gaussian_pdf((double)(i-j), (double)h);
				Rcnt++;
			}
		}
		Rdens[i] = (Rcnt > 0 ? Rsum / (double)(Rcnt * h) : 0.0);
	}
	
	double aveF = minStart * sumF / (double)cntF;
	double aveR = (double)((int)minEnd-ws+intWs) * sumR / (double)cntR;

	double sumSdS, sumSdE;
	sumSdS = 0.0;
	sumSdE = 0.0;
	
	for (int i = lowLagF; i < sizeF - highLagF; i++)
		sumSdS += (double)dataF[i] * pow((double)(i+(int)minStart) - aveF, 2);

	for (int i = lowLagR; i < sizeR - highLagR; i++)
		sumSdE += (double)dataR[i] * pow((double)(i+(int)minEnd-ws+intWs) - aveR, 2);

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

	bool debug = false;
		
	for (int i = lowLagF; i < sizeF-highLagF; i++)
	{
		if (Fdens[i] > 0.0)
		{
			noF = false;
			for (int j = lowLagR; j < sizeR-highLagR; j++)
			{
				if (Rdens[j] > 0.0)
				{
					noR = false;
					dist = (j-lowLagR+(int)minEnd-ws+intWs)-(i-lowLagF+(int)minStart);

					if (debug)
						cerr << i << ", " << j << ", " << (i-lowLagF+(int)minStart) << ", " << (j-lowLagR+(int)minEnd-ws+intWs);
					
					if (dist < minD)
					{
						if (debug)
							cerr << " lower" << endl;
						continue;
					}
					else if (dist > maxD)
					{
						if (debug)
							cerr << " higher" << endl;
						break;
					}
								
					// Keep track of 'best' pair for maxHigh and maxLow
					//Old:
					//val = Fdens[i] * Rdens[j] * gsl_ran_poisson_pdf(dist, lambdaD);
					val = Fdens[i] * Rdens[j] * distDens[dist-minLag];

					if (debug)
						cerr << ", " << Fdens[i] << ", " << Rdens[j] << ", " << dataF[i] << ", " << dataR[j] << ", " << val;
					
					if (val > 0.0)
					{
						// Best candidate so far?
						if (debug)
							cerr << ", val: " << val;
						if (val > best)
						{
							best = val;
							maxLow = (double)(i-lowLagF+(int)minStart);
							maxHigh = (double)(j-lowLagR+(int)minEnd-ws+intWs);
							nobests = 1;

							if (debug)
								cerr << " best";
						}
						// Keep track of the number of equal best candidates
						else if (val == best)
						{
							nobests++;
							if (debug)
								cerr << " best";
						}
					}
					if (debug)
						cerr << endl;
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

	int tmpC = (int)((predS+predE)*0.5 + 0.5);
	if(tmpC == 814006)
	{
	    cout<<endl<<"Starts,ends & score: "<<endl;
	    for (unsigned int i = 0;i< startD.size();i++)
			cout<<startD[i]<<"\t"<<endD[i]<<"\t"<<scoreD[i]<<endl;
	    cout<<endl;
	}
	writePrediction(predS, predE, (split ? "D": "C"),
					score,
					sqrt(sumSdS / (double)cntF),
					sqrt(sumSdE / (double)cntR),
					nPreds);
}

// void NPPredictor::predict(set<double> *predictions, bool split, int verbose)
// {
// 	map<double, pair<Info*, multimap<double, double, distComp>*>, distComp>::iterator distIt;
// 	set<double>::iterator predIt;
// 	int nPreds;
// 	double minStart, maxStart, minEnd, maxEnd;

// 	nPreds = (int)predictions->size();
// 	predIt = predictions->begin();
	
// 	distIt = distMap.find((*predIt));
// 	minStart = maxStart = (*distIt).second.first->get(0);
// 	minEnd = maxEnd = (*distIt).second.first->get(1);

// 	removeEntry((*predIt),verbose);
// 	predictions->erase(predIt);

// 	predIt++;
	
// 	// Go through predictions
// 	if (verbose)
// 		cerr << "Merging predictions" << endl;

// 	while (predIt != predictions->end())
// 	{
// 	  distIt = distMap.find((*predIt));

// 	  minStart = min(minStart, (*distIt).second.first->get(0));
// 	  maxStart = max(maxStart, (*distIt).second.first->get(0));
// 	  minEnd = min(minEnd, (*distIt).second.first->get(1));
// 	  maxEnd = max(maxEnd, (*distIt).second.first->get(1));
	  
// 	  removeEntry((*predIt),verbose);
// 	  predictions->erase(predIt);

// 	  predIt++;
// 	}
	
// 	int size = (int)maxStart-(int)minStart+ws-intWs;

// 	//cerr << "BOUNDARIES: " << minStart << " " << maxStart << " " << minEnd << " " << maxEnd << " " << size << " " << nPreds << endl;
// 	//cerr << "BOUNDARY start: " << (int)minStart << ", " << (int)minStart+size << endl;
// 	//cerr << "BOUNDARY end: " << (int)minEnd - ws + intWs << ", " <<  (int)minEnd - ws + intWs + size << endl;
// 	//cerr << "minPos: " << minPos << endl;

// 	int h = 5;
// 	int lowLag = h, highLag = h;
// 	lowLag = (minStart - minPos - lowLag < 0 ? minStart - minPos : lowLag);
// 	highLag = (maxPos - maxEnd - highLag < 0 ? maxPos - maxEnd : highLag);

// 	size = size + lowLag + highLag;
	
// 	unsigned short *dataF = new unsigned short[size];
// 	unsigned short *dataR = new unsigned short[size];
	
// 	infF->seekg(sizeof(unsigned short) * ((int)minStart - minPos - lowLag + 1) + sizeof(int), ios::beg);
// 	infF->read((char*)dataF, sizeof(unsigned short) * size);
	
// 	infR->seekg(sizeof(unsigned short) * ((int)minEnd - ws + intWs - minPos - lowLag + 1) + sizeof(int), ios::beg);
// 	infR->read((char*)dataR, sizeof(unsigned short) * size);

// 	// Parzen's window kernel density estimation using a gaussian kernel and bandwith h (5)
// 	double* Fdens = new double[size];
// 	double* Rdens = new double[size];
// 	double Fsum;
// 	double Rsum;
// 	int Fcnt = 0;
// 	int Rcnt = 0;

// 	double sumF = 0;
// 	double sumR = 0;
// 	int cntF = 0;
// 	int cntR = 0;
	
// 	for (int i = 0; i < size; i++)
// 	{
// 		Fsum = 0.0;
// 		Rsum = 0.0;
// 		Fcnt = 0.0;
// 		Rcnt = 0.0;

// 		if ((i >= lowLag) && (i < size - highLag))
// 		{
// 			sumF += (double)dataF[i] * ((double)(i-lowLag)/minStart + 1.0);
// 			sumR += (double)dataR[i] * ((double)(i-lowLag)/(double)((int)minEnd-ws+intWs) + 1.0);
// 			cntF += (int)dataF[i];
// 			cntR += (int)dataR[i];
// 		}
		
// 		for (int j = max(0,i-h); j < min(i+h,size); j++)
// 		{
// 			if ((int)dataF[j] > 0)
// 			{
// 				Fsum += (double)dataF[j] * gsl_ran_gaussian_pdf((double)(i-j), (double)h);
// 				Fcnt++;
// 			}

// 			if ((int)dataR[j] > 0)
// 			{
// 				Rsum += (double)dataR[j] * gsl_ran_gaussian_pdf((double)(i-j), (double)h);
// 				Rcnt++;
// 			}
// 		}
// 		Fdens[i] = (Fcnt > 0 ? Fsum / (double)(Fcnt * h) : 0.0);
// 		Rdens[i] = (Rcnt > 0 ? Rsum / (double)(Rcnt * h) : 0.0);
// 	}
	
// 	double aveF = minStart * sumF / (double)cntF;
// 	double aveR = (double)((int)minEnd-ws+intWs) * sumR / (double)cntR;

// 	double sumSdS, sumSdE;
// 	sumSdS = 0.0;
// 	sumSdE = 0.0;
	
// 	for (int i = lowLag; i < size - highLag; i++)
// 	{
// 		sumSdS += (double)dataF[i] * pow((double)(i+(int)minStart) - aveF, 2);
// 		sumSdE += (double)dataR[i] * pow((double)(i+(int)minEnd-ws+intWs) - aveR, 2);
// 	}

// 	//cerr << aveF << " " << aveR << " " << sumSdS << " " << sumSdE << " " << cntF << " " << cntR << endl;
	
// 	double maxHigh = GSL_NAN;
// 	double maxLow = GSL_NAN;
// 	double best = GSL_NEGINF;
// 	double val;
// 	int nobests = 0;
// 	int dist;
// 	// Find the most probable (first occurrence) low and high pair
// 	bool noF = true;
// 	bool noR = true;
// 	for (int i = lowLag; i < size-highLag; i++)
// 	{
// 		if (Fdens[i] > 0.0)
// 		{
// 			noF = false;
// 			for (int j = lowLag; j < size-highLag; j++)
// 			{
// 				if (Rdens[j] > 0.0)
// 				{
// 					noR = false;
// 					dist = (j-lowLag+(int)minEnd-ws+intWs)-(i-lowLag+(int)minStart);
					
// 					//cerr << (i-lowLag+(int)minStart) << ", " << (j-lowLag+(int)minEnd-ws+intWs);
					
// 					if (dist < minD)
// 					{
// 						//cerr << " lower" << endl;
// 						continue;
// 					}
// 					else if (dist > maxD)
// 					{
// 						//cerr << " higher" << endl;
// 						break;
// 					}
								
// 					// Keep track of 'best' pair for maxHigh and maxLow
// 					//Old:
// 					//val = Fdens[i] * Rdens[j] * gsl_ran_poisson_pdf(dist, lambdaD);
// 					val = Fdens[i] * Rdens[j] * distDens[dist-minLag];

// 					//cerr << ", " << Fdens[i] << ", " << Rdens[j] << ", " << dataF[i] << ", " << dataR[j] << ", "<< val;
					
// 					if (val > 0.0)
// 					{
// 						// Best candidate so far?
// 						if (val > best)
// 						{
// 							best = val;
// 							maxLow = (double)(i-lowLag+(int)minStart);
// 							maxHigh = (double)(j-lowLag+(int)minEnd-ws+intWs);
// 							nobests = 1;

// 							//cerr << " best";
// 						}
// 						// Keep track of the number of equal best candidates
// 						else if (val == best)
// 						{
// 							nobests++;
// 							//cerr << " best";
// 						}
// 						//cerr << endl;
// 					}
// 				}
// 			}
// 		}
// 	}

// 	//if (nobests > 1)
// 	//	cerr << nobests << " equal best candidates, using the first one only" << endl;

// 	double predS = maxLow;
// 	double predE = maxHigh;
	
// 	// TODO:
// 	sumF = 0.0;
// 	sumR = 0.0;
// 	double oddsF = 0.5;
// 	double oddsR = 0.5;
// 	int mid = (int)((predS+predE)*0.5 + 0.5);
		
// 	//double oddsF = log(max(gsl_ran_poisson_pdf(i, lambdaF[0]) * pF[0], DBL_MIN)) -
// 	//	log(max(gsl_ran_poisson_pdf(i, lambdaF[1]) * pF[1], DBL_MIN));
// 	//double oddsR = log(max(gsl_ran_poisson_pdf(i, lambdaR[0]) * pR[0], DBL_MIN)) -
// 	//	log(max(gsl_ran_poisson_pdf(i, lambdaR[1]) * pR[1], DBL_MIN));
		
// 	delete[] dataF;
// 	delete[] dataR;
// 	delete[] Fdens;
// 	delete[] Rdens;

// 	// Write prediction
// 	if (verbose)
// 		cerr << "Writing prediction" << endl;

// 	//if (!gsl_isnan(predS) && !gsl_isnan(predE))
// 	writePrediction(predS, predE, (split ? "D": "C"),
// 					oddsF+oddsR,
// 					sqrt(sumSdS / (double)cntF),
// 					sqrt(sumSdE / (double)cntR),
// 					nPreds);
// }

void NPPredictor::printDistMap()
{
	map<double, pair<Info*, multimap<double, double, distComp>*>, distComp>::iterator distIt;
	multimap<double, double, distComp>::iterator it;
	for (distIt = distMap.begin(); distIt != distMap.end(); distIt++)
	{
		(*distIt).second.first->print();
		for (it = (*distIt).second.second->begin(); it !=(*distIt).second.second->end(); it++)
		{
			cerr << (*it).first << "\t" << (*it).second << endl;
		}
	}
}

void NPPredictor::removeEntry(double id, int verbose)
{
	map<double, pair<Info*, multimap<double, double, distComp>*>, distComp>::iterator distIt;
	map<double, pair<Info*, multimap<double, double, distComp>*>, distComp>::iterator eraseIt;
	pair<multimap<double, double, scoreComp>::iterator, multimap<double, double, scoreComp>::iterator> range;
	pair<multimap<double,Info*,comp>::iterator, multimap<double,Info*,comp>::iterator> rangeEnds;
	multimap<double, double, distComp>::iterator it;
	multimap<double, double, scoreComp>::iterator scoreIt;
	multimap<double,Info*,comp>::iterator endsIt;
	double score;
	double end;

	// Remove entries from distMap:
	if (verbose)
		cerr << " Removing entries from distMap (" << distMap.size() << ")" << endl;
	for (distIt = distMap.begin(); distIt != distMap.end(); distIt++)
	{
		// Entry matches id, remove the whole entry
		if ((*distIt).first == id)
		{
			if (verbose)
				cerr << "match" << endl;
			// Clear the info and distances to the other entries within this entry
			score = (*distIt).second.first->get(2);
			end = (*distIt).second.first->get(1);
			//if (verbose)
			//	cerr << " Deleting info";
			//(*distIt).second.first->print();
			delete (*distIt).second.first;
			(*distIt).second.second->clear();
			if (verbose)
				cerr << " Deleting distmultimap";
			delete (*distIt).second.second;

			eraseIt = distIt;
		}
		// Entry doesn't match id, remove the distance to the specific entry from the entry
		else
		{
			//if (verbose)
			//	cerr << "no match for id " << (*distIt).first << endl;
			for (it = (*distIt).second.second->begin(); it != (*distIt).second.second->end(); it++)
			{
				if ((*it).second == id)
				{
					(*distIt).second.second->erase(it);
					break;
				}
			}
		}
	}

	if (verbose)
		cerr << "Deleting entry " << (*eraseIt).first << " (" << id << ") from distMap" << endl;
	distMap.erase(eraseIt);
	
	// Remove entry from scoreMap:
	if (verbose)
		cerr << " Removing entries from scoreMap";
	range = scoreMap.equal_range(score);
	for (scoreIt = range.first; scoreIt != range.second; scoreIt++)
	{
		if ((*scoreIt).second == id)
		{
			scoreMap.erase(scoreIt);
			break;
		}
	}

	// Remove entry from ends
	if (verbose)
		cerr << " Removing entries from ends";
	rangeEnds = ends.equal_range(end);
	for (endsIt = rangeEnds.first; endsIt != rangeEnds.second; endsIt++)
	{
		if ((*endsIt).second->get(2) == id)
		{
			delete (*endsIt).second;
			ends.erase(endsIt);
			break;
		}
	}
	
}

void NPPredictor::flush(int verbose)
{
	double id;
	Info *info;
	multimap<double, double, distComp>::iterator distIt;
	map<double, pair<Info*, multimap<double, double, distComp>*>, distComp>::iterator entry;
	multimap<double, double, distComp> *dist;
	set<double> preds;
	
	double minStart, maxStart, minEnd, maxEnd;
	bool split;

	//printDistMap();
	
	// Flush all predictions
	if (verbose)
		cerr << "Flushing predictions" << endl;
	while (!scoreMap.empty())
	{
		id = (*scoreMap.begin()).second;
		entry = distMap.find(id);
		if (verbose)
			cerr << id << (entry != distMap.end() ? "found" : "not found") << endl;
		minStart = maxStart = (*entry).second.first->get(0);
		minEnd = maxEnd = (*entry).second.first->get(1);

		preds.insert(id);
		split = false;
		dist = (*entry).second.second;
		if (verbose)
			cerr << "Looping" << endl;
		for (distIt = dist->begin(); distIt != dist->end(); distIt++)
		{
			// Valid member to current cluster?
			info = (*distMap.find((*distIt).second)).second.first;
			if (member(info, minStart, maxStart, minEnd, maxEnd))
			{
				// Update min and max coordinates
				minStart = min(minStart, info->get(0));
				maxStart = max(maxStart, info->get(0));
				minEnd = min(minEnd, info->get(1));
				maxEnd = max(maxEnd, info->get(1));

				preds.insert(info->get(3));
			}
			else
			{
				split = (info->get(0) < maxEnd) && (info->get(1) > minStart);
				break;
			}
		}
		if (verbose)
			cerr << "Making prediction from predictions (" << preds.size() << ", " << scoreMap.size() << ")" << endl;
		
		// Make prediction from predictions
		predict(&preds, split, verbose);
		
		preds.clear();
	}
}

void NPPredictor::forcedFlush(int predStart,int verbose)
{
	if(verbose) cerr<<"Forced Flush! \t";
	double id;
	Info *info;
	multimap<double, double, distComp>::iterator distIt;
	map<double, pair<Info*, multimap<double, double, distComp>*>, distComp>::iterator entry;
	multimap<double, double, distComp> *dist;
	set<double> preds;
  
	double minStart, maxStart, minEnd, maxEnd;
	bool split;
	bool exit = false;

	//printDistMap();
	if(verbose) cerr<<"Before: "<<scoreMap.size()<<"\t";
	// Flush the first set of predictions (until first divided-point)
	if (verbose)
		cerr << "Flushing predictions" << endl;
	while (!scoreMap.empty() && !exit)
    {
		id = (*scoreMap.begin()).second;
		entry = distMap.find(id);
		if (verbose)
			cerr << id << (entry != distMap.end() ? "found" : "not found") << endl;
		minStart = maxStart = (*entry).second.first->get(0);
		minEnd = maxEnd = (*entry).second.first->get(1);
      
		preds.insert(id);
		split = false;
		dist = (*entry).second.second;
		if (verbose)
			cerr << "Looping" << endl;
		for (distIt = dist->begin(); distIt != dist->end(); distIt++)
		{
			// Valid member to current cluster?
			info = (*distMap.find((*distIt).second)).second.first;
			if (member(info, minStart, maxStart, minEnd, maxEnd))
			{
				// Update min and max coordinates
				minStart = min(minStart, info->get(0));
				maxStart = max(maxStart, info->get(0));
				minEnd = min(minEnd, info->get(1));
				maxEnd = max(maxEnd, info->get(1));
	      
				preds.insert(info->get(3));
			}
			else
			{
				split = (info->get(0) < maxEnd) && (info->get(1) > minStart);
				exit = true;
				break;
			}
		}
		if (verbose)
			cerr << "Making prediction from predictions (" << preds.size() << ", " << scoreMap.size() << ")" << endl;
      
		// Make prediction from predictions
		predict(&preds, split, verbose);
      
		preds.clear();
    }
    if(verbose)cerr<<"After: "<<scoreMap.size()<<endl;
}

void NPPredictor::writePrediction(double start, double end, const char *type, double* score, double startSd, double endSd, int support)
{
	int debug = 0;
	// Center coordinate
	(*outf) << (int)((start+end)*0.5 + 0.5) << "\t";
	if(debug) cerr << (int)((start+end)*0.5 + 0.5) << "\t";
  
	// Prediction type
	(*outf) << type << "\t";
	if(debug) cerr << type << "\t";
  
	// Start coordiante
	(*outf) << (int)(start + 0.5) << "\t";
	if(debug) cerr << (int)(start + 0.5) << "\t";
  
	// End coordiante
	(*outf) << (int)(end + 0.5) << "\t";
	if(debug) cerr << (int)(end + 0.5) << "\t";
  
	// Score
	(*outf) << score[0] << "\t";
	if(debug) cerr << score[0] << "\t";
	(*outf) << score[1] << "\t";
	if(debug) cerr << score[1] << "\t";
	(*outf) << score[2] << "\t";
	if(debug) cerr << score[2] << "\t";
  
	// Start standard deviation
	(*outf) << startSd << "\t";
	if(debug) cerr << startSd << "\t";
  
	// End standard deviation
	(*outf) << endSd << "\t";
	if(debug) cerr << endSd << "\t";
  
	// Support
	(*outf) << support << endl;
	if(debug) cerr << support << endl;
}
