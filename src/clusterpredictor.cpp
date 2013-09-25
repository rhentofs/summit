#include "clusterpredictor.h"

ClusterPredictor::ClusterPredictor(ofstream *outf, int lambdaD, int ws, int minD, int maxD) :
	Predictor(outf, lambdaD, ws, minD, maxD)
{
	this->id = 0;
}

ClusterPredictor::~ClusterPredictor()
{
}

int ClusterPredictor::addPrediction(double predStart, double predEnd, double predStartSd,
				    double predEndSd, double postIO, int pos, int verbose)
{

  // Only add if proper prediction (non-nan)
  
  if (!gsl_isnan(predStart) && !gsl_isnan(predEnd) && !gsl_isnan(postIO))
    { 
      //cerr<<"ClusterPredictor::aP:\t"<<postIO<<"\t"<<pos<<endl;
      if (verbose)
	{
	  cerr << "Adding prediction" << endl;
	  cerr << pos << " " <<predStart << " " << predEnd << " " << predStartSd << " " <<predEndSd << endl;
	  cerr << "Before: distMap.size() == " << distMap.size() << " ends.size() == " << ends.size() << " scoreMap.size() == " << scoreMap.size() << endl;
	}
      
      // Empty multimap: just add prediction
      if (scoreMap.empty())
	{
	  scoreMap.insert(pair<double,double>(postIO,(double)id));
	  ends.insert(createMapPair(predEnd, predStart, postIO, predStartSd, predEndSd, (double)id));
	  populateDistMap(createMapPair((double)id, predStart, predEnd, postIO, predStartSd, predEndSd, (double)id));
	  //printDistMap();
	  id++;
	}
      // Non-empty multimap
      else
	{
	  multimap<double,Info*,comp>::reverse_iterator end = ends.rbegin();
	  
	  
	  
	  // New prediction sufficiently far away from the last old one? 
	  if (pos - (*end).first >= min(ws,maxD))
	    {
	      if (verbose)
		cerr << "pos - (*end).first >= min(ws,maxD): " << pos << " - " << (*end).first << ">= " << min(ws,maxD) << endl;
	      this->flush(verbose);
	    }
	  else
	    {
	      int max_dist = 2*maxD; // completly as hoc!
	      multimap<double,Info*,comp>::iterator start = ends.begin();	  
	      // New prediction sufficiently far away from the many of the old ones? Shaky...
	      if (pos - (*start).first >= max_dist) // distance 'pos' to 'end'-coordinate in first prediction.    
		{
		  if (verbose)
		    {
		      cerr << "pos - (*start).first >= max_dist: " << pos << " - " << (int)(*start).first;
		      cerr << "( = "<< pos - (int)(*start).first<< ") >= " << max_dist << "\t"<< ends.size()<<endl;
		    }
		  this->forcedFlush(predStart,verbose);
		}
	    }


	  scoreMap.insert(pair<double,double>(postIO,(double)id));
	  ends.insert(createMapPair(predEnd, predStart, postIO, predStartSd, predEndSd, (double)id));
	  populateDistMap(createMapPair((double)id, predStart, predEnd, postIO, predStartSd, predEndSd, (double)id));
	  id++;
	}
      
      if (verbose)
      	cerr << "After: distMap.size() == " << distMap.size() << " ends.size() == " << ends.size() << " scoreMap.size() == " << scoreMap.size() << endl;
    }
  return(0);
}

void ClusterPredictor::populateDistMap(pair<double,Info*> p)
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

void ClusterPredictor::printDistMap()
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

void ClusterPredictor::removeEntry(double id, int verbose)
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
		if ((*endsIt).second->get(4) == id)
		{
			delete (*endsIt).second;
			ends.erase(endsIt);
			break;
		}
	}
	
}

void ClusterPredictor::flush(int verbose)
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

				preds.insert(info->get(5));
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

void ClusterPredictor::forcedFlush(int predStart,int verbose)
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
	      
	      preds.insert(info->get(5));
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

void ClusterPredictor::predict(set<double> *predictions, bool split, int verbose)
{
	map<double, pair<Info*, multimap<double, double, distComp>*>, distComp>::iterator distIt;
	multiset<double> predStarts;
	multiset<double> predEnds;
	multiset<double>::iterator it;
	set<double>::iterator predIt;
	int nPreds;

	double sumPredS, sumPredE, avgS;
	double sumSdS, sumSdE;
	double predS, predE;

	double maxPostIO = log(DBL_MAX)-log(DBL_MIN);

	nPreds = (int)predictions->size();
	sumPredS = 0.0;
	sumPredE = 0.0;
	avgS = 0.0;
	
	// Go through predictions
	if (verbose)
	  cerr << "Merging predictions" << endl;
	for (predIt = predictions->begin(); predIt != predictions->end(); predIt++)
	{
	  distIt = distMap.find((*predIt));
	  
	  sumPredS += ((*distIt).second.first->get(2)/maxPostIO) * (*distIt).second.first->get(0);
	  
	  sumPredE += ((*distIt).second.first->get(2)/maxPostIO) * (*distIt).second.first->get(1);
	  
	  //cerr<<(int)((*distIt).second.first->get(2) > 0);
	  avgS += (*distIt).second.first->get(2)/maxPostIO;
	  
	  predStarts.insert((*distIt).second.first->get(0));
	  predEnds.insert((*distIt).second.first->get(1));
	  
	  removeEntry((*predIt),verbose);
	  predictions->erase(predIt);
	}	
	//cerr<<endl<<"clsuterPredictor::predict:\t"<<nPreds<<"\t"<<sumPredS<<"\t"<<sumPredE<<"\t"<<avgS<<"\t"<<sumPredS/avgS<<"\t"<<sumPredE/avgS<<endl;			
	predS = sumPredS / avgS;
	predE = sumPredE / avgS;
		
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
