#include "predictor.h"

Predictor::Predictor(ofstream *outf, int lambdaD, int ws, int minD, int maxD) :
	lambdaD(lambdaD), ws(ws), minD(minD), maxD(maxD)
{
	this->outf = outf;
	this->scnt = this->ecnt = 0;
}

Predictor::~Predictor()
{
	// Check if all regions have been flushed
	if (!starts.empty() || !ends.empty())
	{
	  multimap<double,Info*,comp>::iterator tmpIt;
	  
	  for (tmpIt = starts.begin(); tmpIt != starts.end(); tmpIt++)
	    delete (*tmpIt).second;
	  
	  for (tmpIt = ends.begin(); tmpIt != ends.end(); tmpIt++)
	    delete (*tmpIt).second;
	  
	  starts.clear();
	  ends.clear();
	}
}

int Predictor::addPrediction(double predStart, double predEnd, double predStartSd,
							 double predEndSd, double postIO, int pos, int verbose)
{
  // Only add if proper prediction (non-nan)
  if (!gsl_isnan(predStart) && !gsl_isnan(predEnd) && !gsl_isnan(postIO))
    {
      // Empty multimap: just add prediction
      if (starts.empty())
	{
	  starts.insert(createMapPair(predStart, predEnd, postIO, predStartSd, predEndSd, (double)scnt));
	  ends.insert(createMapPair(predEnd, predStart, postIO, predStartSd, predEndSd, (double)ecnt));
	  scnt++;
	  ecnt++;
	}
      // Non-empty multimap
      else
	{
	  multimap<double,Info*,comp>::reverse_iterator end = ends.rbegin();
	  
	  // New prediction sufficiently far away from the old ones?
	  if (pos - end->first >= min(ws,maxD))
	    {
	      this->flush(verbose);
	    }
	  starts.insert(createMapPair(predStart, predEnd, postIO, predStartSd, predEndSd, (double)scnt));
	  ends.insert(createMapPair(predEnd, predStart, postIO, predStartSd, predEndSd, (double)ecnt));
	  scnt++;
	  ecnt++;
	}
    }
  
  return(0);
}

pair<double,Info*> Predictor::createMapPair(double key, double first, double second, double third, double fourth, double fifth, double sixth)
{
	pair<double,Info*> p;

	p.first = key;
	p.second = new Info(first, second, third, fourth, fifth, sixth);
	
	return(p);
}

void Predictor::writePrediction(double start, double end, const char *type, double score, double startSd, double endSd, int support)
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
  (*outf) << score << "\t";
  if(debug) cerr << score << "\t";
  
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

void Predictor::predict(set<double> *predictions, bool split, int verbose)
{
	multimap<double,Info*,comp>::iterator startIt, endIt, tmpIt;
	set<double> predsTmpStart, predsTmpEnd;
	multiset<double> predStarts;
	multiset<double> predEnds;
	multiset<double>::iterator it;
	int nPreds;

	double sumS, sumE, sumPredS, sumPredE, avgS;
	double sumSdS, sumSdE;
	double predS, predE;

	double maxPostIO = log(DBL_MAX)-log(DBL_MIN);

	nPreds = (int)predictions->size();
	sumS = 0.0;
	sumE = 0.0;
	sumPredS = 0.0;
	sumPredE = 0.0;
	avgS = 0.0;
	
	predsTmpStart = (*predictions);
	predsTmpEnd = (*predictions);
	
	// Go through starts
	for (startIt = starts.begin(); (!predsTmpStart.empty() && startIt != starts.end());)
	{
		// Prediction in set?
		if (predictions->count((*startIt).second->get(4)) > 0)
		{
			sumPredS += ((*startIt).second->get(1)/maxPostIO)*(*startIt).first;
			sumS += (*startIt).second->get(1)/maxPostIO;
			
			sumPredE += ((*startIt).second->get(1)/maxPostIO)*(*startIt).second->get(0);
			sumE += (*startIt).second->get(1)/maxPostIO;
			
			avgS += (*startIt).second->get(1)/maxPostIO;
			
			predictions->erase((*startIt).second->get(4));
			
			predStarts.insert((*startIt).first);
			predEnds.insert((*startIt).second->get(0));
			
			tmpIt = startIt;
			startIt++;
			predsTmpStart.erase((*tmpIt).second->get(4));
			
			delete (*tmpIt).second;
			starts.erase(tmpIt);
		}
		else
		{
			if (predsTmpStart.count((*startIt).second->get(4)) > 0)
			{
				tmpIt = startIt;
				startIt++;
				predsTmpStart.erase((*tmpIt).second->get(4));
				
				delete (*tmpIt).second;
				starts.erase(tmpIt);
			}
			else
				startIt++;
		}
	}
	
	// Go through ends
	for (endIt = ends.begin(); (!predsTmpEnd.empty() && endIt != ends.end());)
	{
		// Prediction in set?
		if (predictions->count((*endIt).second->get(4)) > 0)
		{
			sumPredE += ((*endIt).second->get(1)/maxPostIO)*(*endIt).first;
			sumE += (*endIt).second->get(1)/maxPostIO;
				
			sumPredS += ((*endIt).second->get(1)/maxPostIO)*(*endIt).second->get(0);
			sumS += (*endIt).second->get(1)/maxPostIO;
				
			avgS += (*endIt).second->get(1)/maxPostIO;
				
			predictions->erase((*endIt).second->get(4));
				
			predEnds.insert((*endIt).first);
			predStarts.insert((*endIt).second->get(0));

			tmpIt = endIt;
			endIt++;
			predsTmpEnd.erase((*tmpIt).second->get(4));

			delete (*tmpIt).second;
			ends.erase(tmpIt);
		}
		else
		{
			if (predsTmpEnd.count((*endIt).second->get(4)) > 0)
			{
				tmpIt = endIt;
				endIt++;
				predsTmpEnd.erase((*tmpIt).second->get(4));

				delete (*tmpIt).second;
				ends.erase(tmpIt);
			}
			else
				endIt++;
		}
	}
				
	predS = sumPredS / sumS;
	predE = sumPredE / sumE;
		
	sumSdS = 0.0;
	sumSdE = 0.0;

	if (verbose)
		cerr << "Calculating standard deviation" << endl;
		
	// Estimate the standard deviation on each side
	for (it = predStarts.begin(); it != predStarts.end(); it++)
		sumSdS += pow((*it) - predS, 2);
	for (it = predEnds.begin(); it != predEnds.end(); it++)
		sumSdE += pow((*it) - predE, 2);
		
	// Write prediction
	if (verbose)
		cerr << "Writing prediction" << endl;
	writePrediction(predS, predE, (split ? "D" : "C"),
					maxPostIO*avgS/(double)nPreds,
					sqrt(sumSdS / (double)predStarts.size()),
					sqrt(sumSdE / (double)predEnds.size()),
					nPreds);
}
