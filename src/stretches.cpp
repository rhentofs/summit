#include "zeb_utils.h"

int identify_stretches(ifstream *ifi, ofstream *ofb, ofstream *ofs, string chrom, int gap, int ms, int odds)
{
	int minPos, maxPos;
	getMinMaxF(ifi, &minPos, &maxPos);

	int read = 0;
	int chunk = CHUNK_MULTI*MAXNR;

	// if necessary, reduce to file size.
	chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);

	storageType *tmp = new storageType[chunk];
	storageType *tmp2 = new storageType[chunk];

	int maxStretch = 1000;
	int stretches[maxStretch];
	for (int i = 0; i < maxStretch; i++)
		stretches[i] = 0;
	
	int start = 0;
	int end = 0;
	int gapCnt = 0;
	bool lastZero = true;
	int gapstart = 0;
	bool gapobserved = false;
	int score = 0;
	int tmpstart = 0;
	int tmpend = 0;
	int len = 0;
	bool split = false;
	
	while ((maxPos - minPos + 1 - read) >= chunk && chunk > 0)
	{
		cerr << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)read / (double)(maxPos - minPos + 1) << "  \t% complete.\r";

		ifi->seekg(sizeof(storageType) * (read) + sizeof(int), ios::beg);
		ifi->read((char*)tmp, sizeof(storageType) * chunk);

		int iter = 0;
		while (iter < chunk)
		{
			split = false;
			if ((int)tmp[iter] >= odds)
			{
				gapCnt = 0;
				// first in strech?
				if(lastZero)
				{
					start = minPos + read + iter;
					end = minPos + read + iter;
					lastZero = false;
					gapobserved = false;
				}
				else
				{
					end = minPos + read + iter;
				}
			}
			else
			{
				if (gapobserved && ((end - start + 1) >= ms))
					split = true;
				if ((gapCnt > gap) || split)
				{
					gapCnt = 0;
					if(!lastZero || split) //strechEnd 
					{
						tmpstart = start;
						if (split)
							tmpend = gapstart - 1;
						else
							tmpend = end;
						
						score = 0;
						// Read odds data
						ifi->seekg(sizeof(storageType) * (tmpstart-minPos) + sizeof(int), ios::beg);
						ifi->read((char*)tmp2, sizeof(storageType) * (tmpend-tmpstart+1));
						for (int j = 0; j < (tmpend-tmpstart+1); j++)
						{
							score = (tmp2[j] > score ? tmp2[j] : score);
						}

						len = tmpend-tmpstart;
						stretches[(len > maxStretch ? maxStretch-1 : len)]++;
						if (split)
							(*ofb) << chrom << "\t" << tmpstart - 1 << "\t" << tmpend << "\tsplit\t";
						else
							if (gapCnt > 0)
								(*ofb) << chrom << "\t" << tmpstart - 1 << "\t" << tmpend << "\tmerge\t";
							else
								(*ofb) << chrom << "\t" << tmpstart - 1 << "\t" << tmpend << "\tstretch\t";
						(*ofb) << score << endl;
					}

					if (split)
					{
						//cerr << "reset\t";
						iter = gapstart - minPos - read;
						if (iter < 0)
						{
							read += iter;
							if (((maxPos - minPos + 1 - read) < chunk))
								chunk = (maxPos - minPos + 1) - read;
							iter = 0;
							ifi->seekg(sizeof(storageType) * (read) + sizeof(int), ios::beg);
							ifi->read((char*)tmp, sizeof(storageType) * chunk);
						}
						start = 0;
						end = 0;
						gapCnt = 0;
					}
					
					lastZero = true;
					gapobserved = false;
				}
				else
				{
					gapCnt++;
					if (gapCnt == 1)
					{
						gapobserved = true;
						gapstart = minPos + read + iter;
					}
				}
			}			
			iter++;
			//cerr << iter << "\t" << read << endl;
		}
		
		read += chunk;
		
		// Check if we can read in a full chunk of data or simply the rest of the file. 
		if (((maxPos - minPos + 1 - read) < chunk))
			chunk = (maxPos - minPos + 1) - read;
	}

	if (!lastZero)
	{
		tmpstart = start;
		tmpend = end;
		score = 0;
		// Read odds data
		ifi->seekg(sizeof(storageType) * (tmpstart-minPos) + sizeof(int), ios::beg);
		ifi->read((char*)tmp2, sizeof(storageType) * (tmpend-tmpstart+1));
		for (int j = 0; j < (tmpend-tmpstart+1); j++)
		{
			score = (tmp2[j] > score ? tmp2[j] : score);
		}
		len = tmpend-tmpstart;
		stretches[(len > maxStretch ? maxStretch-1 : len)]++;
		if (split)
			(*ofb) << chrom << "\t" << tmpstart - 1 << "\t" << tmpend << "\tsplit\t";
		else
			if (gapCnt > 0)
				(*ofb) << chrom << "\t" << tmpstart - 1 << "\t" << tmpend << "\tmerge\t";
			else
				(*ofb) << chrom << "\t" << tmpstart - 1 << "\t" << tmpend << "\tstretch\t";
		(*ofb) << score << endl;
	}
	
	if (ifi->eof())
		ifi->clear();

	for(int i = 0; i < maxStretch; i++)
		(*ofs) << i+1 << "\t" << stretches[i] << endl;
	
	cerr << "\t\t" << setprecision(3) << setw(4) << 100 << "  \t% complete." << endl;
	
	delete[] tmp;
	delete[] tmp2;
	
	return(1);
}

int vlevel = 3;
int main(int argc, char* argv[]) 
{
	string inf   = "";
	string outf  = "";
	string chrom = "";
	string outs  = "";
	int gap      = 0;
	int ms       = 130;
	int odds     = 1;
	
	string errorLine =  "usage " + 
		string(argv[0]) + 
		" [parameters] \n" +
		"\t-i\t<binary infile>\n" +
		"\t-o\t<bed outfile>\n" +
		"\t-s\t<stats outfile>\n" +
		"\t-g\t<max gap allowed, default = 0> \n" +
		"\t-m\t<max strech allowed, default = 130> \n" +
		"\t-t\t<min odds (threshold) value considered, default = 1> \n" +
		"\t-c\t<chromosome name>";

	vlevel = 3;

	bool fail = false;
	string failmessage = "";
	
	for (int i = 1; i < argc; i++)
	{
	    if(strcmp(argv[i],"-i") == 0)
			inf.assign(argv[++i]);
	    else if(strcmp(argv[i],"-o") == 0)
			outf.assign(argv[++i]);
		else if(strcmp(argv[i],"-s") == 0)
			outs.assign(argv[++i]);
		else if(strcmp(argv[i],"-c") == 0)
			chrom.assign(argv[++i]);
		else if(strcmp(argv[i],"-g") == 0)
			gap = atoi(argv[++i]);
		else if(strcmp(argv[i],"-m") == 0)
			ms = atoi(argv[++i]);
		else if(strcmp(argv[i],"-t") == 0)
			odds = atoi(argv[++i]);
		else
		{
			failmessage.assign("Unknown argument: ");
			failmessage.append(argv[i]);
			failmessage.append("\n");
			fail = true;
		}
	}

	if (strcmp(inf.c_str(), "") == 0)
	{
		failmessage.append("-i must be specified.\n");
		fail = true;
	}

	if (strcmp(outf.c_str(), "") == 0)
	{
		failmessage.append("-o must be specified.\n");
		fail = true;
	}

	if (strcmp(outs.c_str(), "") == 0)
	{
		failmessage.append("-s must be specified.\n");
		fail = true;
	}

	if (strcmp(chrom.c_str(), "") == 0)
	{
		failmessage.append("-c must be specified.\n");
		fail = true;
	}

	if (gap < 0)
	{
		failmessage.append("-g must be non-negative.\n");
		fail = true;
	}

	if (ms < 0)
	{
		failmessage.append("-m must be non-negative.\n");
		fail = true;
	}

	if (odds <= 0)
	{
		failmessage.append("-o must be positive.\n");
		fail = true;
	}
	
	if (fail)
	{
		cerr << endl << failmessage.c_str() << endl << errorLine << endl;
		return(-1);
	}

	ifstream ifi;
	ofstream ofb;
	ofstream ofs;

	ifi.open(inf.c_str(),ios::binary);
	ofb.open(outf.c_str(),ios::trunc);
	ofs.open(outs.c_str(),ios::trunc);

	if (ifi.fail())
	{
		failmessage.append("ERROR: Input file \"");
		failmessage.append(inf.c_str());
		failmessage.append("\" could not be opened, aborting.\n");
		fail = true;
	}
	
	if (ofb.fail())
	{
		failmessage.append("ERROR: Output file \"");
		failmessage.append(outf.c_str());
		failmessage.append("\" could not be created, aborting.\n");
		fail = true;
	}

	if (ofs.fail())
	{
		failmessage.append("ERROR: Output file \"");
		failmessage.append(outs.c_str());
		failmessage.append("\" could not be created, aborting.\n");
		fail = true;
	}
	
	if (fail)
	{
		cerr << endl << failmessage.c_str() << endl << errorLine << endl;
		return(-1);
	}

	ofs << "! Command:";
	for (int i = 0; i < argc; i++)
		ofs << " " << argv[i];
	ofs << endl;

	ofb << "track name=\"" << outf << "\" description=\"\" visibility=2" << endl;

	identify_stretches(&ifi,&ofb,&ofs,chrom,gap,ms,odds);

	ifi.close();
	ofb.close();
	ofs.close();
	
	return(1);
}


	
		

	
