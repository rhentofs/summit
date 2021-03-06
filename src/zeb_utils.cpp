#include "zeb_utils.h"

int countFiles(string infPath)
{
	ifstream tmpf;
	DIR *d = opendir(infPath.c_str());
	struct dirent* drnt;
	int infFCnt = 0;
	if (d)
	{
		while ((drnt = readdir(d)))
		{
			string file = drnt->d_name;
			string id = file.substr(0,1);
			
			if (strcmp(id.c_str(),"F") == 0)
			{
				string filename = infPath + file;
				tmpf.open(filename.c_str(),ios::binary);
				if (!tmpf.fail())
				{
					tmpf.close();
					file.replace(0,1,"R");
					filename = infPath + file;
					tmpf.open(filename.c_str(),ios::binary);
					if (!tmpf.fail())
					{
						tmpf.close();
						infFCnt++;
					}
				}
			}
		}
		closedir(d);
	}

	return(infFCnt);
}

int openFiles(string infPath, ifstream* iff, ifstream* ifr, string* fileNames)
{
	int iter = 0;
	ifstream tmpf;
	DIR* d = opendir(infPath.c_str());
	struct dirent* drnt;
	if (d)
	{
		while ((drnt = readdir(d)))
		{
			string file = drnt->d_name;
			string id = file.substr(0,1);
			
			if (strcmp(id.c_str(),"F") == 0)
			{
				string filename = infPath + file;
				tmpf.open(filename.c_str(),ios::binary);
				if (!tmpf.fail())
				{
					tmpf.close();
					file.replace(0,1,"R");
					string filename2 = infPath + file;
					tmpf.open(filename.c_str(),ios::binary);
					if (!tmpf.fail())
					{
						tmpf.close();
						iff[iter].open(filename.c_str(), ios::binary);
						ifr[iter].open(filename2.c_str(), ios::binary);
						file.replace(0,2,"");
						file.replace(file.size()-4,4,"");
						fileNames[iter] = file;
						iter++;
						
					}
				}
			}
		}
		closedir(d);
	}

	return(iter);
}

int openFilesNames(string infPath, ifstream* iff, ifstream* ifr, string* fileNames, int cnt)
{
	for (int i = 0; i < cnt; i++)
	{
		string fileF = infPath + "F_" + fileNames[i] + ".bin";
		string fileR = infPath + "R_" + fileNames[i] + ".bin";

		iff[i].open(fileF.c_str(),ios::binary);
		ifr[i].open(fileR.c_str(),ios::binary);

		if (iff[i].fail() || ifr[i].fail())
			return(-1);
	}

	return(cnt);
}

int estimateDistance(ifstream* infF, ifstream* infR, ofstream* outf,
					 int* minD, int* maxD, int* lambdaD,
					 int minLag, int maxLag,
					 int minPos, int maxPos,
					 int minR, int minF,
					 int maxF, int maxR,
					 int truncValF, int truncValR, double* distDens,
					 double dens)
{
	int chunk = CHUNK_MULTI*MAXNR;
	// if necessary, reduce to file size.
	chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);

	storageType* tmpF = new storageType[chunk];
	storageType* tmpR = new storageType[chunk];
	
	double meanF = 0.0;
	double meanR = 0.0;
	double sumF = 0.0;
	double sumR = 0.0;

	int read = 0;

	while ((maxPos - minPos + 1) - read >= chunk && chunk > 0)
	{
		// Read data from F file
		infF->seekg(sizeof(storageType) * (read + (minPos-minF)) + sizeof(int), ios::beg);
		infF->read((char*)tmpF, sizeof(storageType) * chunk);
		// Read data from R file
		infR->seekg(sizeof(storageType) * (read + (minPos-minR)) + sizeof(int), ios::beg);
		infR->read((char*)tmpR, sizeof(storageType) * chunk);

		// Truncate the data to truncValF and truncValR
		truncData(tmpF, chunk, truncValF);
		truncData(tmpR, chunk, truncValR);
		
		sumF = 0.0;
		sumR = 0.0;
		
		//Loop all positions
		for (int i = 0; i < chunk; i++)
		{
			sumF += (double)tmpF[i];
			sumR += (double)tmpR[i];
		}

		meanF = (meanF * ((double)read / (double)(read + chunk))) + (sumF / (double)(read + chunk));
		meanR = (meanR * ((double)read / (double)(read + chunk))) + (sumR / (double)(read + chunk));
		
		read += chunk;
		
		// Check if we can read in a full chunk of data or simply the rest of the file. 
		if (((maxPos - minPos + 1 - read) < chunk))
			chunk = (maxPos - minPos + 1) - read;
	}

	// Check if we have read past the EOF so the eofbit is set
	if (infF->eof())
		infF->clear();
	if (infR->eof())
		infR->clear();
	
	chunk = CHUNK_MULTI*MAXNR;
	// if necessary, reduce to file size.
	chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);

	read = 0;

	int maxDist = 2*((maxLag-minLag)+np_round((double)minLag/2.0))+1;
	long cnt[maxDist];

	for (int dist = 0; dist < maxDist; dist++)
		cnt[dist] = 0;
	
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";

	bool posPos[maxDist];

	double avg = 0.0;
	long avgcnt = 0;

	long odds = 0;
	long evens = 0;
	
	while (((maxPos - minPos + 1) - read >= chunk) && (chunk - maxLag >= 0))
	{
		// Read data from F file
		infF->seekg(sizeof(storageType) * (read + (minPos-minF)) + sizeof(int), ios::beg);
		infF->read((char*)tmpF, sizeof(storageType) * chunk);
		// Read data from R file
		infR->seekg(sizeof(storageType) * (read + (minPos-minR)) + sizeof(int), ios::beg);
		infR->read((char*)tmpR, sizeof(storageType) * chunk);

		// Truncate the data to truncValF and truncValR
		truncData(tmpF, chunk, truncValF);
		truncData(tmpR, chunk, truncValR);
		
		//Loop all F positions
		bool fabove = false;
		bool rabove = false;
		int lowPos = -1;
		int highPos = -1;
		int lowNeg = -1;
		int tmp = 0;
		
		for (int j = 0; j <= maxDist; j++)
			posPos[j] = false;
		
		for (int i = 0; i <= chunk - maxLag; i++)
		{
			fabove = (int)tmpF[i] > 0;

			if (fabove)
			{
				if (lowPos == -1)
				{
					lowPos = i;
					highPos = i;
					posPos[0] = true;
				}

				if ((lowNeg != -1) && ((lowNeg - i) < minLag))
				{
					for (int j = 0; j < maxDist; j++)
					{
						if (posPos[j])
						{
							tmp = 2*(highPos-(lowPos+j))+(lowNeg-highPos);
							if (((lowNeg-highPos) % 2) == 0)
								evens++;
							else
								odds++;
							avg = ((avg * (double)avgcnt) + ((double)tmp * (double)tmpF[lowPos+j])) /
								(avgcnt + tmpF[lowPos+j]);
							avgcnt += tmpF[lowPos+j];
							cnt[tmp] += tmpF[lowPos+j];
							//cnt[tmp][tmpF[lowPos+j]] += 1;
							posPos[j] = false;
						}
					}
						
					lowNeg = -1;
					lowPos = i;
				}
				highPos = i;
				posPos[i-lowPos] = true;

				if (lowNeg == -1)
				{
					for (int lag = minLag; lag <= maxLag; lag++)
					{
						rabove = (int)tmpR[i+lag] > 0;
						if (rabove)
						{
							lowNeg = i+lag;
							break;
						}
					}
					if (lowNeg == -1)
						lowPos = -1;
				}
			}
		}

		if (lowNeg != -1)
		{
			for (int j = 0; j < maxDist; j++)
			{
				if (posPos[j])
				{
					tmp = 2*(highPos-(lowPos+j))+(lowNeg-highPos);
					if (((lowNeg-highPos) % 2) == 0)
						evens++;
					else
						odds++;

					avg = ((avg * (double)avgcnt) + ((double)tmp * (double)tmpF[lowPos+j])) /
						(avgcnt + tmpF[lowPos+j]);
					avgcnt += tmpF[lowPos+j];
					cnt[tmp] += tmpF[lowPos+j];
					//cnt[tmp][tmpF[lowPos+j]] += 1;
				}
			}
		}
		
		read += chunk-maxLag+1;
		
		// Check if we can read in a full chunk of data or simply the rest of the file. 
		if (((maxPos - minPos + 1 - read) < chunk) && ((maxPos - minPos + 1 - read) >= maxLag))
			chunk = (maxPos - minPos + 1) - read;
		
		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)read / (double)(maxPos - minPos + 1) << "  \t% complete.\r";
	}
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;

	cerr << "Average fragment length: " << avg << endl;
	cerr << "Sense read counts: " << avgcnt << endl;

	cerr << "lowNegs " << evens << "\t" << odds << endl; 

	
	// Check if we have read past the EOF so the eofbit is set
	if (infF->eof())
		infF->clear();
	if (infR->eof())
		infR->clear();
	
	*(outf) << "lag\t" << "cnt\t" << endl;
	
	for (int i = 0; i < maxDist; i++)
	{
		*(outf) << i << "\t" << cnt[i] << endl;
		//*(outf) << i << "\t";
		//for (int j = 0; j < truncValF-1; j++)
		//	*(outf) << cnt[i][j] << "\t";
		//*(outf) << cnt[i][truncValF-1] << endl;
	}

	/*
	  	double sum = 0.0;
	for (int lag = minLag; lag <= maxLag; lag++)
		sum += cnt[lag];
	
	double maxVal = 0.0;
	for (int lag = minLag; lag <= maxLag; lag++)
	{
		distDens[lag-minLag] = cnt[lag] / sum;
		maxVal = max(maxVal, distDens[lag-minLag]);
	}
	
	double steps[1000];
	for (int i = 0; i < 1000; i++)
		steps[i] = maxVal - ((double)(i-1) * maxVal) / 1000.0;
	
	double peak = 0.0;
	int peakId = -1;
	int start = -1;
	int end = -1;

	double tmpPeak = 0.0;
	int tmpPeakId = -1;
	int tmpStart = -1;
	int tmpEnd = -1;

	for (int i = 0; i < 1000; i++)
	{
		peak = tmpPeak = sum = 0.0;
		peakId = tmpPeakId = start = tmpStart = end = tmpEnd = -1;

		for (int lag = minLag; lag <= maxLag; lag++)
		{
			if (distDens[lag-minLag] >= steps[i])
			{
				if (sum == 0.0)
				{
					sum = distDens[lag-minLag];
					tmpPeak = distDens[lag-minLag];
                    tmpPeakId = lag;
                    tmpStart = lag;
                    tmpEnd = lag;
                }
                else
                {
                    sum += distDens[lag-minLag];
                    tmpEnd = lag;
                    if (distDens[lag-minLag] > tmpPeak)
                    {
                        tmpPeak = distDens[lag-minLag];
                        tmpPeakId = lag;
                    }
                }
            }
            else
            {
                if (tmpPeakId != -1 && sum > 0.0)
                {
                    if (sum >= dens)
                    {
                        if (tmpPeak > peak)
                        {
                            peak = tmpPeak;
                            peakId = tmpPeakId;
                            start = tmpStart;
                            end = tmpEnd;
                        }
                    }
                }
                else
                {
                    if (sum >= dens)
                    {
                        peak = tmpPeak;
                        peakId = tmpPeakId;
                        start = tmpStart;
                        end = tmpEnd;
                    }
                }
                sum = 0.0;
            }
        }
        if (peakId != -1)
            break;
    }
	
	*minD = start;
	*maxD = end;
	*lambdaD = peakId;

	vcerr(3) << "\t\testimated minD: " << *minD << endl;
	vcerr(3) << "\t\testimated maxD: " << *maxD << endl;
	vcerr(3) << "\t\testimated lambdaD: " << *lambdaD << endl;
	*/
	
	delete[] tmpF;
	delete[] tmpR;
	
	return(1);
}


int checkFileLengths(ifstream* infF, ifstream* infR,
					 int* minPos, int* maxPos,
					 int* minR, int* minF,
					 int* maxF, int* maxR,
					 int fCnt)
{
	vcerr(2) << "*** Checking file lengths ***" << endl;

	for (int i=0; i<fCnt; i++)
	{
		getMinMaxF(&infF[i], &minF[i], &maxF[i]);
		getMinMaxF(&infR[i], &minR[i], &maxR[i]);
		
		if ((minF[i] != minR[i]) || (maxF[i] != maxR[i])) // Input files have not the same coordinates
		{
			//vcerr(3) << "\tFiles start at different positions." << endl;
			//vcerr(3) << "Largest common region: ";
			
			minPos[i] = (minF[i] > minR[i] ? minF[i] : minR[i]);
			maxPos[i] = (maxF[i] < maxR[i] ? maxF[i] : maxR[i]);
			
			//vcerr(3) << "[ " << minPos[i] << ", " << maxPos[i] << " ]." << endl;
		}
		else // input files have the exact same coordinates
		{
			minPos[i] = minF[i];
			maxPos[i] = maxF[i];
		}
	}
	
	return(1);
}

int checkFileLengths(ifstream* infF1, ifstream* infF2,
					 ifstream* infR1, ifstream* infR2,
					 int* minPos, int* maxPos,
					 int* minR1, int* minR2,
					 int* minF1, int* minF2,
					 int* maxF1, int* maxF2,
					 int* maxR1, int* maxR2,
					 int fCnt)
{
	vcerr(2) << "*** Checking file lengths ***" << endl;

	for (int i=0; i<fCnt; i++)
	{
		getMinMaxF(&infF1[i], &minF1[i], &maxF1[i]);
		getMinMaxF(&infR1[i], &minR1[i], &maxR1[i]);
		getMinMaxF(&infF2[i], &minF2[i], &maxF2[i]);
		getMinMaxF(&infR2[i], &minR2[i], &maxR2[i]);
		
		if ((minF1[i] != minR1[i]) ||
			(maxF1[i] != maxR1[i]) ||
			(minF2[i] != minR2[i]) ||
			(maxF2[i] != maxR2[i]) ||
			(minF1[i] != minF2[i]) ||
			(minR1[i] != minR2[i]) ||
			(maxF1[i] != maxF2[i]) ||
			(maxR1[i] != maxR2[i])) // Input files have not the same coordinates
		{
			//vcerr(3) << "\tFiles start at different positions." << endl;
			//vcerr(3) << "Largest common region: ";

			int minP = minF1[i];
			int maxP = maxF1[i];

			if (minF2[i] > minP)
				minP = minF2[i];
			if (minR1[i] > minP)
				minP = minR1[i];
			if (minR2[i] > minP)
				minP = minR2[i];

			if (maxF2[i] < maxP)
				maxP = maxF2[i];
			if (maxR1[i] < maxP)
				maxP = maxR1[i];
			if (maxR2[i] < maxP)
				maxP = maxR2[i];
			
			minPos[i] = minP;
			maxPos[i] = maxP;
			
			//vcerr(3) << "[ " << minPos[i] << ", " << maxPos[i] << " ]." << endl;
		}
		else // input files have the exact same coordinates
		{
			minPos[i] = minF1[i];
			maxPos[i] = maxF1[i];
		}
	}
	
	return(1);
}


int calculateTruncVal(ifstream* infF, ifstream* infR, ofstream* outf,
					  int* truncValF, int* truncValR, double truncLimit,
					  int* minPos, int* maxPos,
					  int* minR, int* minF,
					  int* maxF, int* maxR,
					  int fCnt)
{
	
	int maxValF = 0;
	int maxValR = 0;
	int read = 0;
	int chunk = CHUNK_MULTI*MAXNR;
	storageType *tmpF, *tmpR;
	
	double totalLength = 0;
	for (int i = 0; i < fCnt; i++)
		totalLength += (maxPos[i] - minPos[i] + 1) / 1000.0;

	double totalRead = 0.0;

	vcerr(2) << "\t* " << "Looking up max values" << endl;
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";
	for (int iter = 0; iter < fCnt; iter++)
	{
		chunk = CHUNK_MULTI*MAXNR;
		// if necessary, reduce to file size.
		chunk = ((maxPos[iter] - minPos[iter] + 1) < chunk ? (maxPos[iter] - minPos[iter] + 1) : chunk);
		
		tmpF = new storageType[chunk];
		tmpR = new storageType[chunk];
		
		read = 0;
		
		// Identify max values in F and R files
		while ((maxPos[iter] - minPos[iter] + 1 - read) >= chunk && chunk > 0)
		{
			// Read data from F file
			infF[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minF[iter])) + sizeof(int), ios::beg);
			infF[iter].read((char*)tmpF, sizeof(storageType) * chunk);
			// Read data from R file
			infR[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minR[iter])) + sizeof(int), ios::beg);
			infR[iter].read((char*)tmpR, sizeof(storageType) * chunk);
			
			//Loop all positions
			for (int i = 0; i < chunk; i++)
			{
				maxValF = ((int)tmpF[i] > maxValF ? (int)tmpF[i] : maxValF);
				maxValR = ((int)tmpR[i] > maxValR ? (int)tmpR[i] : maxValR);
			}
			
			read += chunk;
			
			// Check if we can read in a full chunk of data or simply the rest of the file. 
			if (((maxPos[iter] - minPos[iter] + 1 - read) < chunk))
				chunk = (maxPos[iter] - minPos[iter] + 1) - read;
		}
		
		totalRead += (double)read / 1000.0;
		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * totalRead / totalLength << "  \t% complete.\r";

		// Check if we have read past the EOF so the eofbit is set
		if (infF[iter].eof())
			infF[iter].clear();
		if (infR[iter].eof())
			infR[iter].clear();

		delete[] tmpF;
		delete[] tmpR;
	}
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;
	vcerr(3) << "\t\tmaxValF: " << maxValF << " maxValR: " << maxValR << endl;

	// count vectors
	long countsF[maxValF+1];
	long countsR[maxValR+1];
	long countF = 0;
	long countR = 0;

	for (int i = 0; i < maxValF; i++)
		countsF[i] = 0;
	for (int i = 0; i < maxValR; i++)
		countsR[i] = 0;

	totalRead = 0.0;
	
	vcerr(2) << "\t* " << "Counting value ocurrences" << endl;
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";
	for (int iter = 0; iter < fCnt; iter++)
	{
		chunk = CHUNK_MULTI*MAXNR;
		// if necessary, reduce to file size.
		chunk = ((maxPos[iter] - minPos[iter] + 1) < chunk ? (maxPos[iter] - minPos[iter] + 1) : chunk);
		
		tmpF = new storageType[chunk];
		tmpR = new storageType[chunk];
		
		read = 0;
	
		// Populate count vectors
		while ((maxPos[iter] - minPos[iter] + 1) - read >= chunk && chunk > 0)
		{
			// Read data from F file
			infF[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minF[iter])) + sizeof(int), ios::beg);
			infF[iter].read((char*)tmpF, sizeof(storageType) * chunk);
			// Read data from R file
			infR[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minR[iter])) + sizeof(int), ios::beg);
			infR[iter].read((char*)tmpR, sizeof(storageType) * chunk);
			
			//Loop all positions
			for (int i = 0; i < chunk; i++)
			{
				countsF[(int)tmpF[i]]++;
				countsR[(int)tmpR[i]]++;
				countF += (long)tmpF[i];
				countR += (long)tmpR[i];
			}
			
			read += chunk;
			
			// Check if we can read in a full chunk of data or simply the rest of the file. 
			if (((maxPos[iter] - minPos[iter] + 1 - read) < chunk))
				chunk = (maxPos[iter] - minPos[iter] + 1) - read;
		}

		totalRead += (double)read / 1000.0;
		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * totalRead / totalLength << "  \t% complete.\r";
		
		// Check if we have read past the EOF so the eofbit is set
		if (infF[iter].eof())
			infF[iter].clear();
		if (infR[iter].eof())
			infR[iter].clear();

		delete[] tmpF;
		delete[] tmpR;
	}
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;
	
	if(truncLimit == 1)
	{
	    *truncValF = maxValF;
	    *truncValR = maxValR;
	}
	else
	{
		// Identify truncate values
		vcerr(2) << "\t* " << "Identifying truncate values" << endl;
	  
		double tmp = 0.0;
		*truncValF = 0;
		for (int i = 1; i <= maxValF; i++)
	    {
			tmp = ((tmp * (double)countF)  + ((double)i * (double)countsF[i])) / (double)countF;
			if (tmp > truncLimit)
				break;
			(*truncValF)++;
	    }
		
		*truncValR = 0;
		tmp = 0.0;
		for (int i = 1; i <= maxValR; i++)
	    {
			tmp = ((tmp * (double)countR)  + ((double)i * (double)countsR[i])) / (double)countR;
			if (tmp > truncLimit)
				break;
			(*truncValR)++;
	    }
		vcerr(3) << "\t\t" << "truncValF: " << (*truncValF) << " truncValR: " << (*truncValR) << endl;
	}
	
	// Write counts to outf
	for (int i = 1; i <= max(maxValF, maxValR); i++)
		(*outf) << i << "\t" << (i > maxValF ? 0 : countsF[i]) << "\t" << (i > maxValR ? 0 : countsR[i]) << endl;

	return(1);
}


int truncData(storageType* data, int size, int truncVal)
{
	storageType tv = (storageType)truncVal;
	for (int i = 0; i < size; i++)
		data[i] = (data[i] > tv ? tv : data[i]);
	return(1);
}

int estimateParameters(ifstream* infF, ifstream* infR, ofstream* outf,
					   int minD, int maxD, int lambdaD,
					   double* pF, double* pR, double* lambdaF, double* lambdaR,
					   int sampleSize, int burnins, int iterations,
					   int seed, int* minPos, int* maxPos,
					   int* minR, int* minF,
					   int* maxF, int* maxR,
					   int truncValF, int truncValR, int fCnt)
{
	/*
	  Assign the data range variables
	*/

	// Window size of reads
	int w = ceil((double)(maxD - minD) / 2.0);
	
	// Offsets from center position
	int thetaF1, thetaF2, thetaR1, thetaR2;
	thetaF1 = thetaR2 = ceil((double)(maxD - 1) / 2.0);
	thetaF2 = thetaR1 = thetaF1 - w + 1;
	maxD = thetaF1 + thetaR2 + 1;

	bool wholeSample[fCnt];
		
	// Data range
	int dataRange[fCnt];
	// Maximum number of non-overlapping maxDs
	int maxDs[fCnt];
	int maxDsTotal = 0;
	int sampleSizeTotal = 0;
	int sampleSizes[fCnt];
	bool validSize[fCnt];
	for (int iter = 0; iter < fCnt; iter++)
	{
		dataRange[iter] = maxPos[iter] - minPos[iter] + 1;
		maxDs[iter] = dataRange[iter] / (maxD + 1);
		maxDsTotal += maxDs[iter];
		// Valid sample size?
		wholeSample[iter] = (sampleSize < 1);
		sampleSizes[iter] = sampleSize;
		if (wholeSample[iter])
			sampleSizes[iter] = maxDs[iter];
		validSize[iter] = (sampleSizes[iter] <= maxDs[iter]);
		if (!validSize[iter])
		{
			sampleSizes[iter] = maxDs[iter];
			wholeSample[iter] = true;
		}
		sampleSizeTotal += sampleSizes[iter];
	}
	
	vcerr(3) << "\t* Settings" << endl;
	vcerr(3) << setprecision(3);
	vcerr(3) << "\t\tRead count window length: " << w << endl;
	vcerr(3) << "\t\tMaximum number of non-overlapping regions: " << maxDsTotal << endl;
	for (int iter = 0; iter < fCnt; iter++)
		if (!validSize[iter])
			vcerr(3) << "\t\tSample size too big, reset to whole sample." << endl;
	vcerr(3) << "\t\tNumber of non-overlapping regions in sample: " << sampleSizeTotal << endl;
	//vcerr(3) << " " << (wholeSample ? "(whole sample)" : "(user specified)") << endl;
	vcerr(3) << "\t\tOffset forward reads: [" << -thetaF1 << "," << -thetaF2 << "]" << endl;
	vcerr(3) << "\t\tOffset reverse reads: [" << thetaR1 << "," << thetaR2 << "]" << endl;
	
	/*
	  Sample the data used for estimation of parameters and calculate sample statistics
	*/
	
	vcerr(2) << "\t* Sampling data" << endl;
	
	// Set seed to random int if not specified (-1)
	if (seed == -1)
	{
		srand(time(NULL));
		seed = rand();
	}
	
	// Allocate vector of positions on which to base parameter estimation
	int* centers = new int[sampleSizeTotal];
	int centersIter = 0;
	for (int iter = 0; iter < fCnt; iter++)
	{
		if (!wholeSample[iter])
		{
			vcerr(3) << "\t\tSampling center positions" << endl;
			int* centerstmp = new int[sampleSizes[iter]];
			np_sample(seed, minPos[iter]+thetaF1, maxPos[iter]-thetaR2, maxD+1, sampleSizes[iter], centerstmp);
			for (int i = 0; i < sampleSizes[iter]; i++)
			{
				centers[centersIter] = centerstmp[i];
				centersIter++;
			}
			delete[] centerstmp;
		}
		else
		{
			for (int i=0; i < sampleSizes[iter]; i++)
			{
				centers[centersIter] = i * (maxD+1) + minPos[iter] + thetaF1;
				centersIter++;
			}
		}
	}
	
	vcerr(3) << "\t\tCalculating window sums" << endl;
	int* F = new int[sampleSizeTotal];
	int* R = new int[sampleSizeTotal];
	double meanF, meanR, varF, varR, covFR;
	int maxValF, maxValR;

 	if (windowSums(infF, infR, centers, sampleSizes, maxPos, minPos,
 				   minF, minR, minD, maxD,
 				   &maxValF, &maxValR,
 				   thetaF1, thetaF2, thetaR1, thetaR2,
 				   F, R, &meanF, &meanR, &varF, &varR, &covFR,
 				   truncValF, truncValR, fCnt) < 0)
 	{
 		delete[] centers;
 		delete[] F;
 		delete[] R;
 		cerr << endl << "Error: windowSums failed, aborting." << endl;
 		return(-1);
 	}
		
	vcerr(3) << "\t\tSampling statistics:" << endl;
	vcerr(3) << "\t\tF: mean = " << meanF << " variance = " << varF << endl;
	vcerr(3) << "\t\tR: mean = " << meanR << " variance = " << varR << endl;
	vcerr(3) << "\t\tCovariance = " << covFR << endl;
	vcerr(3) << "\t\tF: max = " << maxValF << " R: max = " << maxValR << endl;
	
	delete[] centers;

 	/*
 	  Estimate parameters through MCMC sampling
	*/

 	vcerr(2) << "\t* Estimating parameters" << endl;
		
 	// Allocations
 	bool* ZF = new bool[sampleSizeTotal];
 	bool* ZR = new bool[sampleSizeTotal];
	
 	// Parameters
 	//double alphaF = pow(meanF,2) / (varF - pow(meanF, 2));
 	//double alphaR = pow(meanR,2) / (varR - pow(meanR, 2));
 	double alphaF = pow(meanF,2) / (varF - meanF);
 	double alphaR = pow(meanR,2) / (varR - meanR);
 	double betaF = alphaF / meanF;
 	double betaR = alphaR / meanR;
 	double deltaF[2] = {1.0,1.0};
 	double deltaR[2] = {1.0,1.0};
	
 	// Non-informative priors
 	deltaF[0] = 1.0;
 	deltaR[0] = 1.0;
 	deltaF[1] = 1.0;
 	deltaR[1] = 1.0;

 	// Informative priors
 	/*
	  deltaF[0] = 1.0;
	  deltaR[0] = 1.0;
	  deltaF[1] = np_round(((double)lambdaD * (2.0 - 0.8)) / (double)w) - 1;
	  deltaR[1] = np_round(((double)lambdaD * (2.0 - 0.8)) / (double)w) - 1;
 	*/
	
 	// Initiate allocations with naive guess based on F and R means
 	for (int i = 0; i < sampleSize; i++)
 	{
 		ZF[i] = (F[i] > meanF * 3.0);
 		ZR[i] = (R[i] > meanR * 3.0);
 	}

 	vcerr(2) << "\t\t Gibbs sampling" << endl;
 	// Gibbs sampling
 	if (gibbs(F, R, ZF, ZR,
 			  pF, pR,
 			  lambdaF, lambdaR,
 			  deltaF, deltaR,
 			  alphaF, alphaR, betaF, betaR,
 			  maxValF, maxValR,
 			  sampleSizeTotal, iterations, burnins, seed,
 			  outf) < 0)
 	{
 		delete[] F;
 		delete[] R;
 		delete[] ZF;
 		delete[] ZR;

 		cerr << endl << "ERROR: gibbs failed, aborting." << endl;
 		return(-1);
 	}

 	delete[] F;
 	delete[] R;
	delete[] ZF;
	delete[] ZR;
	
	return(1);
}

int windowSums(ifstream* infF, ifstream* infR, int* centers,
			   int* sampleSizes, int* maxPos, int* minPos,
			   int* minF, int* minR, int minD, int maxD,
			   int* maxValF, int* maxValR,
			   int thetaF1, int thetaF2, int thetaR1, int thetaR2,
			   int* F, int* R,
			   double* meanF, double* meanR,
			   double* varF, double* varR,
			   double* covFR,
			   int truncValF, int truncValR, int fCnt)
{
	int chunk;
	storageType* tmpF;
	storageType* tmpR;

	*meanF = 0.0;
	*meanR = 0.0;

	*maxValF = 0;
	*maxValR = 0;
	
	int sampleCount;
	int totalSampleCount = 0;
	int read = 0;
	int w = thetaF1 - thetaF2 + 1;

	int totalSampleSize = 0;
	for (int i = 0; i < fCnt; i++)
		totalSampleSize += sampleSizes[i];
	
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";
	for (int iter = 0; iter < fCnt; iter++)
	{
		chunk = CHUNK_MULTI*MAXNR;
		// if necessary, reduce to file size.
		chunk = ((maxPos[iter] - minPos[iter] + 1) < chunk ? (maxPos[iter] - minPos[iter] + 1) : chunk);

		tmpF = new storageType[chunk];
		tmpR = new storageType[chunk];

		read = 0;
		sampleCount = 0;
		
		// Assign Fi and Ri
		while (((maxPos[iter] - minPos[iter] + 1) - read >= chunk) && (chunk - maxD >= 0) && (sampleSizes[iter] > sampleCount))
		{
			if (centers[totalSampleCount] + thetaR2 <= minPos[iter] + read + chunk)
			{				
				// Read data from F file
				infF[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minF[iter])) + sizeof(int), ios::beg);
				infF[iter].read((char*)tmpF, sizeof(storageType) * chunk);
				// Read data from R file
				infR[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minR[iter])) + sizeof(int), ios::beg);
				infR[iter].read((char*)tmpR, sizeof(storageType) * chunk);
				
				// Truncate the data to truncValF and truncValR
				truncData(tmpF, chunk, truncValF);
				truncData(tmpR, chunk, truncValR);
				
				while ((centers[totalSampleCount] < minPos[iter] + read + chunk) && (sampleSizes[iter] > sampleCount) && (centers[totalSampleCount] + thetaR2 <= minPos[iter] + read + chunk))
				{
					F[totalSampleCount] = 0;
					R[totalSampleCount] = 0;
					
					for (int j = 0; j < w; j++)
					{
						F[totalSampleCount] += tmpF[centers[totalSampleCount] - thetaF1 + j - read - minPos[iter]];
						R[totalSampleCount] += tmpR[centers[totalSampleCount] + thetaR1 + j - read - minPos[iter]];
					}
					
					if (F[totalSampleCount] > (*maxValF))
						*maxValF = F[totalSampleCount];
					
					if (R[totalSampleCount] > (*maxValR))
						*maxValR = R[totalSampleCount];
					
					// Update the mean
					*meanF = ((*meanF) * ((double)totalSampleCount / (double)(totalSampleCount+1))) +
						(F[totalSampleCount] / (double)(totalSampleCount+1));
					*meanR = ((*meanR) * ((double)totalSampleCount / (double)(totalSampleCount+1))) +
						(R[totalSampleCount] / (double)(totalSampleCount+1));
					
					sampleCount++;
					totalSampleCount++;
				}
			}
			read += chunk-thetaF1-1;

			// Check if we can read in a full chunk of data or simply the rest of the file. 
			if (((maxPos[iter] - minPos[iter] + 1 - read) < chunk) && ((maxPos[iter] - minPos[iter] + 1 - read) >= maxD))
				chunk = (maxPos[iter] - minPos[iter] + 1) - read;
			
			vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)totalSampleCount / (double)totalSampleSize << "  \t% complete.\r";
		}
		delete[] tmpF;
		delete[] tmpR;
		
		// Check if we have read past the EOF so the eofbit is set
		if (infF[iter].eof())
			infF[iter].clear();
		if (infR[iter].eof())
			infR[iter].clear();
	}
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;
	
	// Calculate variance and covariance
	*varF = 0.0;
	*varR = 0.0;
	*covFR = 0.0;
	double mult;
	for (int i = 0; i < totalSampleSize; i++)
	{
		mult = ((double)i / (double)(i+1));
		*varF = ((*varF) * mult) + (pow((double)F[i] - (*meanF), 2) / (double)(i+1));
		*varR = ((*varR) * mult) + (pow((double)R[i] - (*meanF), 2) / (double)(i+1));
		*covFR = ((*covFR) * mult) + (((double)F[i] - (*meanF)) * ((double)R[i] - (*meanR)) / (double)(i+1));
	}
	
	return(1);
}

int gibbs(int* F, int* R, bool* ZF, bool* ZR,
		  double* pF, double* pR,
		  double* lambdaF, double* lambdaR,
		  double* deltaF, double* deltaR,
		  double alphaF, double alphaR, double betaF, double betaR,
		  int maxValF, int maxValR,
		  int sampleSize, int iterations, int burnins, int seed,
		  ofstream* outf)
{
	gsl_rng* r;
	
	r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, seed);

	long nFS = 0;
	long nFNotS = 0;
	long nRE = 0;
	long nRNotE = 0;

	long FSumS = 0.0;
	long FSumNotS = 0.0;
	long RSumE = 0.0;
	long RSumNotE = 0.0;
	
	for (int i = 0; i < sampleSize; i++)
	{
		if (ZF[i])
		{
			nFS++;
			FSumS += (long)F[i];
		}
		else
		{
			nFNotS++;
			FSumNotS += (long)F[i];
		}

		if (ZR[i])
		{
			nRE++;
			RSumE += (long)R[i];
		}
		else
		{
			nRNotE++;
			RSumNotE += (long)R[i];
		}
	}

	pF[0] = 0.0;
	pF[1] = 0.0;
	pR[0] = 0.0;
	pR[1] = 0.0;
	lambdaF[0] = 0.0;
	lambdaF[1] = 0.0;
	lambdaR[0] = 0.0;
	lambdaR[1] = 0.0;

	double pFWork[2];
	double pRWork[2];
	double lambdaFWork[2];
	double lambdaRWork[2];
	
	double tmp[2];
	double mult;
	bool criteriaMet;

	double posteriorFS[maxValF];
	double posteriorFNotS;
	double posteriorRE[maxValR];
	double posteriorRNotE;
	
	(*outf) << "pS\tpNotS\tpE\tPNotE\tlambdaS\tlambdaNotS\tlambdaE\tlambdaNotE" << endl;

	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";
	// Iterate burnins and iterations
	for (int iter = 0; iter < burnins + iterations; iter++)
	{
	  vcerr(3) << iter << " ";
	    //vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)iter / (double)(burnins + iterations) << "  \t% complete.\r";
		/*
		 *  1. Update of unknowns conditional on allocations
		 */

		//cerr << endl << 1 << endl;
		
		// Update of mixing proportions
		tmp[0] = deltaF[0] + (double)nFS;
		tmp[1] = deltaF[1] + (double)nFNotS;
		gsl_ran_dirichlet(r, 2, tmp, pFWork);

		// Require same mixing proportions for F and R
		pRWork[0] = pFWork[0];
		pRWork[1] = pFWork[1];

		//tmp[0] = deltaR[0] + (double)nRE;
		//tmp[1] = deltaR[1] + (double)nRNotE;
		//gsl_ran_dirichlet(r, 2, tmp, pR);

		// Update of lambdas
		criteriaMet = false;
		// TODO: Make sure that this loop don't get stuck somehow
		while (!criteriaMet)
		{
			lambdaFWork[0] = gsl_ran_gamma(r, alphaF + (double)FSumS, 1.0 / (betaF + (double)nFS));
			lambdaFWork[1] = gsl_ran_gamma(r, alphaF + (double)FSumNotS, 1.0 / (betaF + (double)nFNotS));

			criteriaMet = (lambdaFWork[0] > lambdaFWork[1]);
			if (!criteriaMet)
			  {
			    double tmp = lambdaFWork[0];
			    lambdaFWork[0] = lambdaFWork[1];
			    lambdaFWork[1] = tmp;
			    criteriaMet = true;
			  }
		}

		criteriaMet = false;
		// TODO: Make sure that this loop don't get stuck somehow
		while (!criteriaMet)
		{
			lambdaRWork[0] = gsl_ran_gamma(r, alphaR + (double)RSumE, 1.0 / (betaR + (double)nRE));
			lambdaRWork[1] = gsl_ran_gamma(r, alphaR + (double)RSumNotE, 1.0 / (betaR + (double)nRNotE));

			criteriaMet = (lambdaRWork[0] > lambdaRWork[1]);
			if (!criteriaMet)
			  {
			    double tmp = lambdaRWork[0];
			    lambdaRWork[0] = lambdaRWork[1];
			    lambdaRWork[1] = tmp;
			    criteriaMet = true;
			  }
		}

		/*
		 *  2. Update of allocations conditional on unknowns
		 */

		//cerr << 2 << endl;

		for (int i = 0; i < maxValF; i++)
		{
			posteriorFS[i] = gsl_ran_poisson_pdf(i, lambdaFWork[0]) * pFWork[0];
			posteriorFNotS = gsl_ran_poisson_pdf(i, lambdaFWork[1]) * pFWork[1];

			if ((posteriorFS[i] == 0) && (posteriorFNotS == 0))
			{
				if (i > lambdaFWork[0])
					posteriorFS[i] = 1.0;
			}
			else
				posteriorFS[i] = posteriorFS[i] / (posteriorFS[i] + posteriorFNotS);

			posteriorFS[i] = max(min(posteriorFS[i], 1.0 - DBL_MIN), DBL_MIN);
		}

		for (int i = 0; i < maxValR; i++)
		{
			posteriorRE[i] = gsl_ran_poisson_pdf(i, lambdaRWork[0]) * pRWork[0];
			posteriorRNotE = gsl_ran_poisson_pdf(i, lambdaRWork[1]) * pRWork[1];

			if ((posteriorRE[i] == 0) && (posteriorRNotE == 0))
			{
				if (i > lambdaRWork[0])
					posteriorRE[i] = 1.0;
			}
			else
				posteriorRE[i] = posteriorRE[i] / (posteriorRE[i] + posteriorRNotE);

			posteriorRE[i] = max(min(posteriorRE[i], 1.0 - DBL_MIN), DBL_MIN);
		}
		
		nFS = 0;
		nFNotS = 0;
		nRE = 0;
		nRNotE = 0;
		FSumS = 0;
		FSumNotS = 0;
		RSumE = 0;
		RSumNotE = 0;

		for (int i = 0; i < sampleSize; i++)
		{
			//cerr << i << endl;
			
			ZF[i] = (gsl_ran_binomial(r, posteriorFS[F[i]], 1) == 1);
			ZR[i] = (gsl_ran_binomial(r, posteriorRE[R[i]], 1) == 1);
			
			if (ZF[i])
			{
				nFS++;
				FSumS += (long)F[i];
			}
			else
			{
				nFNotS++;
				FSumNotS += (long)F[i];
			}
			
			if (ZR[i])
			{
				nRE++;
				RSumE += (long)R[i];
			}
			else
			{
				nRNotE++;
				RSumNotE += (long)R[i];
			}
		}

		// Update parameter estimates
		if (iter >= burnins)
		{
			mult = (double)(iter-burnins) / (double)(iter-burnins+1);
			pF[0] = (mult * pF[0]) + (pFWork[0] / (double)(iter-burnins+1));
			pF[1] = (mult * pF[1]) + (pFWork[1] / (double)(iter-burnins+1));
			pR[0] = (mult * pR[0]) + (pRWork[0] / (double)(iter-burnins+1));
			pR[1] = (mult * pR[1]) + (pRWork[1] / (double)(iter-burnins+1));
			lambdaF[0] = (mult * lambdaF[0]) + (lambdaFWork[0] / (double)(iter-burnins+1));
			lambdaF[1] = (mult * lambdaF[1]) + (lambdaFWork[1] / (double)(iter-burnins+1));
			lambdaR[0] = (mult * lambdaR[0]) + (lambdaRWork[0] / (double)(iter-burnins+1));
			lambdaR[1] = (mult * lambdaR[1]) + (lambdaRWork[1] / (double)(iter-burnins+1));
		}
		
		// Write parameters to file
		(*outf) << pFWork[0] << "\t" << pFWork[1] << "\t" << pRWork[0] << "\t" << pRWork[1] << "\t";
		(*outf) << lambdaFWork[0] << "\t" << lambdaFWork[1] << "\t" << lambdaRWork[0] << "\t" << lambdaRWork[1] << endl;
	}
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete.\r";
	gsl_rng_free(r);
	
	return(1);
}

int readParameters(ifstream* ipff, double* pF, double* pR, double* lambdaF, double* lambdaR,int* lambdaD,
				   int* truncValF, int* truncValR,
				   int* minD,int* maxD, int* estLambdaD,int* estMinD,int* estMaxD,vector<string>* corrv)
{
	string line,id;
	vector<string> data;
	size_t strLen;
 
	getline(*ipff,line); 
	data.clear();
	Tokenize (line,data,"\t ");
	while(!ipff->eof())
    {
		id = data[0];
		strLen = id.length();
		if(id.substr(0,1) == "#") //value
		{
			id = id.substr(1,strLen-1);
			if (strcmp(id.c_str(),"pFpeak")==0)
				pF[0] = atof(data[1].c_str());
			else if (strcmp(id.c_str(),"pFnoise")==0)
				pF[1] = atof(data[1].c_str());
			else if (strcmp(id.c_str(),"pRpeak")==0)
				pR[0] = atof(data[1].c_str());
			else if (strcmp(id.c_str(),"pRnoise")==0)
				pR[1] = atof(data[1].c_str());
			else if (strcmp(id.c_str(),"lambdaFpeak")==0)
				lambdaF[0] = atof(data[1].c_str());
			else if (strcmp(id.c_str(),"lambdaFnoise")==0)
				lambdaF[1] = atof(data[1].c_str());
			else if (strcmp(id.c_str(),"lambdaRpeak")==0)
				lambdaR[0] = atof(data[1].c_str());
			else if (strcmp(id.c_str(),"lambdaRnoise")==0)
				lambdaR[1] = atof(data[1].c_str());
			else if (strcmp(id.c_str(),"lambdaD")==0)
				*lambdaD = atoi(data[1].c_str());
			else if (strcmp(id.c_str(),"truncValF")==0)
				*truncValF = atoi(data[1].c_str());
			else if (strcmp(id.c_str(),"truncValR")==0)
				*truncValR = atoi(data[1].c_str());
			else if (strcmp(id.c_str(),"minD")==0)
				*minD = atoi(data[1].c_str());
			else if (strcmp(id.c_str(),"maxD")==0)
				*maxD = atoi(data[1].c_str());
			else if (strcmp(id.c_str(),"estLambdaD")==0)
				*estLambdaD = atoi(data[1].c_str());
			else if (strcmp(id.c_str(),"estMinD")==0)
				*estMinD = atoi(data[1].c_str());
			else if (strcmp(id.c_str(),"estMaxD")==0)
				*estMaxD = atoi(data[1].c_str());
			else{
				cerr<<"Unknown value parameter specifyer: "<<id<<endl;
				return(-1);
			}
		}
      
		if(id.substr(0,1) == "@") //array
		{
			id = id.substr(1,strLen-1);
			if(strcmp(id.c_str(),"corr") == 0){
				Tokenize(data[1],*corrv,",");
			}else{
				cerr<<"Unknown array parameter specifyer: "<<id<<endl;
				return(-1);
			}
		}
		getline(*ipff,line); 
		data.clear();
		Tokenize (line,data,"\t ");
    } // while
  
	// make sure that the 'minD' and 'maxD' are in synch with the length of 'corr'
	if((int)corrv->size() != (*maxD-*minD +1))
    {
		cerr<<"Number of correlationpoints ("<<corrv->size()<<") does not match [mind,maxd] interval ("<<*maxD-*minD +1<<" points )."<<endl;
		return(-1);
    }
	return(1);
}

int writeParameter(ofstream* opff,string marker,string value)
{
	(*opff) << marker << "\t" << value << endl;
	return(1);
}

string doubleToString(double d)
{
	ostringstream s;
	s << d;
	return(string(s.str()));
}

string intToString(int i)
{
	ostringstream s;
	s << i;
	return(string(s.str()));
}


int writeParameters(ofstream* opff, double* pF, double* pR, double* lambdaF, double* lambdaR, int lambdaD,
					int truncValF, int truncValR,
					int minD,int maxD, int estLambdaD,int estMinD,int estMaxD,double* corr)
{
  
	writeParameter(opff,string("#pFpeak"),doubleToString(pF[0]));
	writeParameter(opff,string("#pFnoise"),doubleToString(pF[1]));
  
	writeParameter(opff,string("#pRpeak"),doubleToString(pR[0]));
	writeParameter(opff,string("#pRnoise"),doubleToString(pR[1]));
  
	writeParameter(opff,string("#lambdaFpeak"),doubleToString(lambdaF[0]));
	writeParameter(opff,string("#lambdaFnoise"),doubleToString(lambdaF[1]));
  
	writeParameter(opff,string("#lambdaRpeak"),doubleToString(lambdaR[0]));
	writeParameter(opff,string("#lambdaRnoise"),doubleToString(lambdaR[1]));
 
	writeParameter(opff,string("#lambdaD"),intToString(lambdaD));

	writeParameter(opff,string("#truncValF"),intToString(truncValF));
	writeParameter(opff,string("#truncValR"),intToString(truncValR));
	
	writeParameter(opff,string("#minD"),intToString(minD));
	writeParameter(opff,string("#maxD"),intToString(maxD));
  
	writeParameter(opff,string("#estLambdaD"),intToString(estLambdaD));
	writeParameter(opff,string("#estMinD"),intToString(estMinD));
	writeParameter(opff,string("#estMaxD"),intToString(estMaxD));
  

	string toWrite(doubleToString(corr[0]));
	for (int i = 1;i<=(maxD-minD);i++)
    {
		toWrite = toWrite + "," + doubleToString(corr[i]); 
    }
	writeParameter(opff,string("@corr"),toWrite);
  
	return(1);
}

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
			   int fCnt, string* fileNames)
{
	/*
	 * 1. Assign the data range variables
	 */

	// Window size of reads
	int w = ceil((double)(maxD - minD) / 2.0);
	// Offsets from center position
	int thetaF1, thetaF2, thetaR1, thetaR2;
	thetaF1 = thetaR2 = ceil((double)(maxD - 1) / 2.0);
	thetaF2 = thetaR1 = thetaF1 - w + 1;
	maxD = thetaF1 + thetaR2 + 1;

	int thetaF0, thetaR0;
	thetaF0 = thetaR0 = thetaF1 + thetaF2;
	
	/*
	 * 2. Slide the files and calculate odds of center positions
	 */
	
	// current window sums.
	long wsumF1 = 0, wsumR1 = 0;
	long wsumF2 = 0, wsumR2 = 0; 

	storageType* tmpF1;
	storageType* tmpR1;
	storageType* tmpF2;
	storageType* tmpR2;

	int read = 0;
	double oddsF;
	double oddsR;

	double oddsF1;
	double oddsR1;
	double oddsF2;
	double oddsR2;
	
	int maxValF1 = 0, maxValF2 = 0, maxValR1 = 0, maxValR2 = 0;
	
	int chunk;
	
	double totalLength = 0;
	for (int i = 0; i < fCnt; i++)
		totalLength += (maxPos[i] - minPos[i] + 1) / 1000.0;
	
	double totalRead = 0.0;
	
	vcerr(2) << "\t* Looking up max values" << endl;
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";

	for (int iter = 0; iter < fCnt; iter++)
	{
		chunk = CHUNK_MULTI*MAXNR;
		// if necessary, reduce to file size.
		chunk = ((maxPos[iter] - minPos[iter] + 1) < chunk ? (maxPos[iter] - minPos[iter] + 1) : chunk);
		
		tmpF1 = new storageType[chunk];
		tmpR1 = new storageType[chunk];
		tmpF2 = new storageType[chunk];
		tmpR2 = new storageType[chunk];

		read = 0;
		
		while (((maxPos[iter] - minPos[iter] + 1) - read >= chunk) && (chunk - (thetaF0+thetaR0+1) >= 0))
		{
			// Read data from F file
			infF1[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minF1[iter])) + sizeof(int), ios::beg);
			infF1[iter].read((char*)tmpF1, sizeof(storageType) * chunk);
			infF2[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minF2[iter])) + sizeof(int), ios::beg);
			infF2[iter].read((char*)tmpF2, sizeof(storageType) * chunk);
			// Read data from R file
			infR1[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minR1[iter])) + sizeof(int), ios::beg);
			infR1[iter].read((char*)tmpR1, sizeof(storageType) * chunk);
			infR2[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minR2[iter])) + sizeof(int), ios::beg);
			infR2[iter].read((char*)tmpR2, sizeof(storageType) * chunk);
			
			// Truncate the data to truncValF and truncValR
			truncData(tmpF1, chunk, truncValF1);
			truncData(tmpR1, chunk, truncValR1);
			truncData(tmpF2, chunk, truncValF2);
			truncData(tmpR2, chunk, truncValR2);

			// loop over the read data and calculate sums over the window size. 
			wsumF1 = 0;
			wsumR1 = 0;
			wsumF2 = 0;
			wsumR2 = 0;
						
			for (int i = 1; i < chunk-thetaR0; i++)
			{
				if (i >= (thetaF0-w+1))
				{
					wsumF1 += (long)tmpF1[i-thetaF2];      // window sum
					wsumR1 += (long)tmpR1[i+w+thetaR1-1];  // window sum
					wsumF2 += (long)tmpF2[i-thetaF2];      // window sum
					wsumR2 += (long)tmpR2[i+w+thetaR1-1];  // window sum
				}
				
				if (i >= thetaF0)
				{
					maxValF1 = max((int)wsumF1, maxValF1);
					maxValR1 = max((int)wsumR1, maxValR1);
					maxValF2 = max((int)wsumF2, maxValF2);
					maxValR2 = max((int)wsumR2, maxValR2);

					// decrease the window sum with the value just before the current window.
					wsumF1 -= (long)tmpF1[i-thetaF1];
					wsumR1 -= (long)tmpR1[i+thetaR1];
					wsumF2 -= (long)tmpF2[i-thetaF1];
					wsumR2 -= (long)tmpR2[i+thetaR1];
				}
			}
			read += chunk-thetaF0-thetaR0;
			
			// Check if we can read in a full chunk of data or simply the rest of the file. 
			if (((maxPos[iter] - minPos[iter] + 1 - read) < chunk) && ((maxPos[iter] - minPos[iter] + 1 - read) >= (thetaF0+thetaR0+1)))
				chunk = (maxPos[iter] - minPos[iter] + 1) - read;
		}
		delete[] tmpF1;
		delete[] tmpR1;
		delete[] tmpF2;
		delete[] tmpR2;
		
		totalRead += (double)read / 1000.0;
		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * totalRead / totalLength << "  \t% complete.\r";

		// Check if we have read past the EOF so the eofbit is set
		if (infF1[iter].eof())
			infF1[iter].clear();
		if (infR1[iter].eof())
			infR1[iter].clear();
		if (infF2[iter].eof())
			infF2[iter].clear();
		if (infR2[iter].eof())
			infR2[iter].clear();
	}
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;

	(*outs) << "Max F window sum (file 1) = " << maxValF1 << endl;
	(*outs) << "Max R window sum (file 1) = " << maxValR1 << endl;
	(*outs) << "Max F window sum (file 2) = " << maxValF2 << endl;
	(*outs) << "Max R window sum (file 2) = " << maxValR2 << endl;

	vcerr(3) << "\t\tMax F window sum (file 1) = " << maxValF1 << endl;
	vcerr(3) << "\t\tMax R window sum (file 1) = " << maxValR1 << endl;
	vcerr(3) << "\t\tMax F window sum (file 2) = " << maxValF2 << endl;
	vcerr(3) << "\t\tMax R window sum (file 2) = " << maxValR2 << endl;
	
	double postFs1[maxValF1];
	double postRs1[maxValR1];
	double postFs2[maxValF2];
	double postRs2[maxValR2];
	double postNotFs1[maxValF1];
	double postNotRs1[maxValR1];
	double postNotFs2[maxValF2];
	double postNotRs2[maxValR2];

	double likeFs1[maxValF1];
	double likeNotFs2[maxValF2];
	double likeRs1[maxValR1];
	double likeNotRs2[maxValR2];

	for (int i = 0; i < maxValF1; i++)
	{
		postFs1[i] = gsl_ran_poisson_pdf(i, lambdaF1[0]) * pF1[0];
		postNotFs1[i] = gsl_ran_poisson_pdf(i, lambdaF1[1]) * pF1[1];
		likeFs1[i] = gsl_ran_poisson_pdf(i, lambdaF1[0]);
	}
	for (int i = 0; i < maxValF2; i++)
	{
		postFs2[i] = gsl_ran_poisson_pdf(i, lambdaF2[0]) * pF2[0];
		postNotFs2[i] = gsl_ran_poisson_pdf(i, lambdaF2[1]) * pF2[1];
		likeNotFs2[i] = gsl_ran_poisson_pdf(i, lambdaF2[1]);
	}
	for (int i = 0; i < maxValR1; i++)
	{
		postRs1[i] = gsl_ran_poisson_pdf(i, lambdaR1[0]) * pR1[0];
		postNotRs1[i] = gsl_ran_poisson_pdf(i, lambdaR1[1]) * pR1[1];
		likeRs1[i] = gsl_ran_poisson_pdf(i, lambdaR1[0]);
	}
	for (int i = 0; i < maxValR2; i++)
	{
		postRs2[i] = gsl_ran_poisson_pdf(i, lambdaR2[0]) * pR2[0];
		postNotRs2[i] = gsl_ran_poisson_pdf(i, lambdaR2[1]) * pR2[1];
		likeNotRs2[i] = gsl_ran_poisson_pdf(i, lambdaR2[1]);
	}

	double priorPosF = pF1[0] * pF2[1];
	double priorPosR = pR1[0] * pR2[1];
	double priorNegF = (pF1[0] * pF2[0]) + (pF1[1] * pF2[1]) + (pF1[1] * pF2[0]);
	double priorNegR = (pR1[0] * pR2[0]) + (pR1[1] * pR2[1]) + (pR1[1] * pR2[0]);
		
	bool FoddsPos = false;
	bool RoddsPos = false;
	
	int pos = thetaF1-1;
	int written = thetaF1;
	bool newchunk = false;
	
	vcerr(2) << "\t* Calculating odds of change" << endl;	
	for (int iter = 0; iter < fCnt; iter++)
	{
		vcerr(3) << "\t  " << fileNames[iter] << endl;
		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";
		
		read = 0;
		chunk = MAXNR;
		chunk = ((maxPos[iter] - minPos[iter] + 1) < chunk ? (maxPos[iter] - minPos[iter] + 1) : chunk);

		tmpF1 = new storageType[chunk];
		tmpR1 = new storageType[chunk];
		tmpF2 = new storageType[chunk];
		tmpR2 = new storageType[chunk];
		
		storageType discreteOdds[chunk-thetaF0-thetaR0];
		short discreteDiffs[chunk-thetaF0-thetaR0];
		storageType discreteOddsStart[chunk-thetaF0-thetaR0];
		storageType discreteOddsEnd[chunk-thetaF0-thetaR0];
		
		// Write min coordinate to odds files:
		outo[iter].write((char*)&minPos[iter], sizeof(int));
		outd[iter].write((char*)&minPos[iter], sizeof(int));

		if (bodds)
		{
			outos[iter].write((char*)&minPos[iter], sizeof(int));
			outoe[iter].write((char*)&minPos[iter], sizeof(int));
		}
		
		// Fill with zeros up to thetaF1 - 1
		storageType zero = 0;
		//for (int i = 0; i < thetaF1; i++)
		for (int i = 0; i < thetaF0; i++)
		{
			outo[iter].write((char*)&zero, sizeof(storageType));
			outd[iter].write((char*)&zero, sizeof(storageType));
			if (bodds)
			{
				outos[iter].write((char*)&zero, sizeof(storageType));
				outoe[iter].write((char*)&zero, sizeof(storageType));
			}
		}
		
		pos = thetaF1-1;
		newchunk = false;
		written = thetaF1;
		
		while (((maxPos[iter] - minPos[iter] + 1) - read >= chunk) && (chunk - (thetaF0+thetaR0+1) >= 0))
		{
			// Read data from F file
			infF1[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minF1[iter])) + sizeof(int), ios::beg);
			infF1[iter].read((char*)tmpF1, sizeof(storageType) * chunk);
			infF2[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minF2[iter])) + sizeof(int), ios::beg);
			infF2[iter].read((char*)tmpF2, sizeof(storageType) * chunk);
			// Read data from R file
			infR1[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minR1[iter])) + sizeof(int), ios::beg);
			infR1[iter].read((char*)tmpR1, sizeof(storageType) * chunk);
			infR2[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minR2[iter])) + sizeof(int), ios::beg);
			infR2[iter].read((char*)tmpR2, sizeof(storageType) * chunk);
			
			// Truncate the data to truncValF and truncValR
			truncData(tmpF1, chunk, truncValF1);
			truncData(tmpR1, chunk, truncValR1);
			truncData(tmpF2, chunk, truncValF2);
			truncData(tmpR2, chunk, truncValR2);

			// loop over the read data and calculate sums over the window size. 
			wsumF1 = 0;
			wsumR1 = 0;
			wsumF2 = 0;
			wsumR2 = 0;

			for (int i = 1; i < chunk-thetaR0; i++)
			{
				if (i >= (thetaF0-w+1))
				{
					wsumF1 += (long)tmpF1[i-thetaF2];      // window sum
					wsumR1 += (long)tmpR1[i+w+thetaR1-1];  // window sum
					wsumF2 += (long)tmpF2[i-thetaF2];      // window sum
					wsumR2 += (long)tmpR2[i+w+thetaR1-1];  // window sum
				}
				
				// if we've read the window size the counts are valid.
				if (i >= thetaF0)
				{
					/*
					oddsF = log(max(postFs1[(int)wsumF1] * postNotFs2[(int)wsumF2], DBL_MIN)) -  // 1, not 2
						log(max((postNotFs1[(int)wsumF1] * postNotFs2[(int)wsumF2]) +            // not 1, not 2
								(postFs1[(int)wsumF1] * postFs2[(int)wsumF2]) +                  // 1, 2
								(postNotFs1[(int)wsumF1] * postFs2[(int)wsumF2]), DBL_MIN));     // not 1, 2

					oddsR = log(max(postRs1[(int)wsumR1] * postNotRs2[(int)wsumR2], DBL_MIN)) -  // 1, not 2
						log(max((postNotRs1[(int)wsumR1] * postNotRs2[(int)wsumR2]) +            // not 1, not 2
								(postRs1[(int)wsumR1] * postRs2[(int)wsumR2]) +                  // 1, 2
								(postNotRs1[(int)wsumR1] * postRs2[(int)wsumR2]), DBL_MIN));     // not 1, 2
					*/

					oddsF = log(max(priorPosF, DBL_MIN)) - log(max(priorNegF, DBL_MIN)) +
						log(max(likeFs1[(int)wsumF1], DBL_MIN)) + log(max(likeNotFs2[(int)wsumF2], DBL_MIN)) -
						(log(max((postFs1[(int)wsumF1] * postFs2[(int)wsumF2]) +
								 (postNotFs1[(int)wsumF1] * postNotFs2[(int)wsumF2]) +
								 (postNotFs1[(int)wsumF1] * postFs2[(int)wsumF2]), DBL_MIN)) -
						 log(max(priorNegF, DBL_MIN)));

					oddsR = log(max(priorPosR, DBL_MIN)) - log(max(priorNegR, DBL_MIN)) +
						log(max(likeRs1[(int)wsumR1], DBL_MIN)) + log(max(likeNotRs2[(int)wsumR2], DBL_MIN)) -
						(log(max((postRs1[(int)wsumR1] * postRs2[(int)wsumR2]) +
								 (postNotRs1[(int)wsumR1] * postNotRs2[(int)wsumR2]) +
								 (postNotRs1[(int)wsumR1] * postRs2[(int)wsumR2]), DBL_MIN)) -
						 log(max(priorNegR, DBL_MIN)));
					
					oddsF1 = log(max(postFs1[(int)wsumF1], DBL_MIN)) - log(max(postNotFs1[(int)wsumF1], DBL_MIN));
					oddsF2 = log(max(postFs2[(int)wsumF2], DBL_MIN)) - log(max(postNotFs2[(int)wsumF2], DBL_MIN));
					oddsR1 = log(max(postRs1[(int)wsumR1], DBL_MIN)) - log(max(postNotRs1[(int)wsumR1], DBL_MIN));
					oddsR2 = log(max(postRs2[(int)wsumR2], DBL_MIN)) - log(max(postNotRs2[(int)wsumR2], DBL_MIN));

					FoddsPos = (oddsF > 0.0) && (oddsF1 > 0.0) && (oddsF2 < 0.0);
					RoddsPos = (oddsR > 0.0) && (oddsR1 > 0.0) && (oddsR2 < 0.0);
					
					if (bodds)
					{
						if (FoddsPos)
							discreteOddsStart[i-thetaF0] = np_round(1.0+oddsF);
						else
							discreteOddsStart[i-thetaF0] = 0;
						if (RoddsPos)
							discreteOddsEnd[i-thetaF0] = np_round(1.0+oddsR);
						else
							discreteOddsEnd[i-thetaF0] = 0;
					}

					discreteDiffs[i-thetaF0] = np_round((oddsF1+oddsR1) - (oddsF2+oddsR2));
					
					if (FoddsPos && RoddsPos)
					{
						discreteOdds[i-thetaF0] = np_round(1.0+oddsF+oddsR);
					}
					else
					{
						discreteOdds[i-thetaF0] = 0;
					}
					pos++;
					
					// decrease the window sum with the value just before the current window.
					wsumF1 -= (long)tmpF1[i-thetaF1];
					wsumR1 -= (long)tmpR1[i+thetaR1];
					wsumF2 -= (long)tmpF2[i-thetaF1];
					wsumR2 -= (long)tmpR2[i+thetaR1];
				}
			}
			
			// Write the chunk of calls to file
			if (chunk == MAXNR)
			{
				outo[iter].write((char*)discreteOdds, sizeof(discreteOdds));
				outd[iter].write((char*)discreteOdds, sizeof(discreteDiffs));
				if (bodds)
				{
					outos[iter].write((char*)discreteOddsStart, sizeof(discreteOddsStart));
					outoe[iter].write((char*)discreteOddsEnd, sizeof(discreteOddsEnd));
				}
			}
			else
			{
				for (int i = 0; i < chunk-thetaF0-thetaR0; i++)
				{
					outo[iter].write((char*)&discreteOdds[i], sizeof(storageType));
					outd[iter].write((char*)&discreteDiffs[i], sizeof(short));
					if (bodds)
					{
						outos[iter].write((char*)&discreteOddsStart[i], sizeof(storageType));
						outoe[iter].write((char*)&discreteOddsEnd[i], sizeof(storageType));
					}
				}
			}
			written += chunk-thetaF0-thetaR0;

			if (!newchunk)
			{
				read += chunk-thetaF0-thetaR0;
			}
			else
				read += chunk;
			
			// Check if we can read in a full chunk of data or simply the rest of the file. 
			if (((maxPos[iter] - minPos[iter] + 1 - read) < chunk) && ((maxPos[iter] - minPos[iter] + 1 - read) >= (thetaF0+thetaR0+1)))
			{
				chunk = (maxPos[iter] - minPos[iter] + 1) - read;
				newchunk = true;
			}

			vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)read / (double)(maxPos[iter] - minPos[iter] + 1) << "  \t% complete.\r";
		}

		delete[] tmpF1;
		delete[] tmpR1;
		delete[] tmpF2;
		delete[] tmpR2;
		
		// Fill odds file with ending zeros
		for (int i = pos+1; i < (maxPos[iter]-minPos[iter]+1); i++)
		{
			outo[iter].write((char*)&zero, sizeof(storageType));
			outd[iter].write((char*)&zero, sizeof(storageType));
			if (bodds)
			{
				outos[iter].write((char*)&zero, sizeof(storageType));
				outoe[iter].write((char*)&zero, sizeof(storageType));
			}
		}

		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;
		
		// Check if we have read past the EOF so the eofbit is set
		if (infF1[iter].eof())
			infF1[iter].clear();
		if (infR1[iter].eof())
			infR1[iter].clear();
		if (infF2[iter].eof())
			infF2[iter].clear();
		if (infR2[iter].eof())
			infR2[iter].clear();
	}
	
	return(1);
}


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
					   int fCnt, string* fileNames)
{	
	/*
	 * 1. Assign the data range variables
	 */

	// Window size of reads
	int w = ceil((double)(maxD - minD) / 2.0);
	// Offsets from center position
	int thetaF1, thetaF2, thetaR1, thetaR2;
	thetaF1 = thetaR2 = ceil((double)(maxD - 1) / 2.0);
	thetaF2 = thetaR1 = thetaF1 - w + 1;
	maxD = thetaF1 + thetaR2 + 1;

	int thetaF0, thetaR0;
	thetaF0 = thetaR0 = thetaF1 + thetaF2;
	
	/*
	 * 2. Slide the files and calculate odds of center positions
	 */
	
	double meanF[2] = {0.0, 0.0};
	double meanR[2] = {0.0, 0.0};

	long cntF[2] = {0,0};
	long cntR[2] = {0,0};
		
	// current window sums.
	long wsumF = 0, wsumR = 0; 
	// current window sums.
	long wFsumF = 0, wFsumR = 0; 

	storageType* tmpF;
	storageType* tmpR;

	int read = 0;
	double oddsF;
	double oddsR;
	int idF;
	int idR;

	double maxSDFR = 0.0;
	double wMeanF = 0.0;
	double wMeanR = 0.0;
	double wVarF = 0.0;
	double wVarR = 0.0;
	int wCntF = 0;
	int wCntR = 0;

	double sd;
	
	*maxValF = 0;
	*maxValR = 0;

	int chunk;
	
	double totalLength = 0;
	for (int i = 0; i < fCnt; i++)
		totalLength += (maxPos[i] - minPos[i] + 1) / 1000.0;
	
	double totalRead = 0.0;
	
	vcerr(2) << "\t* Looking up max value" << endl;
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";

	for (int iter = 0; iter < fCnt; iter++)
	{
		chunk = CHUNK_MULTI*MAXNR;
		// if necessary, reduce to file size.
		chunk = ((maxPos[iter] - minPos[iter] + 1) < chunk ? (maxPos[iter] - minPos[iter] + 1) : chunk);
		
		tmpF = new storageType[chunk];
		tmpR = new storageType[chunk];

		read = 0;
		
		//while (((maxPos[iter] - minPos[iter] + 1) - read >= chunk) && (chunk - maxD >= 0))
		while (((maxPos[iter] - minPos[iter] + 1) - read >= chunk) && (chunk - (thetaF0+thetaR0+1) >= 0))
		{
			// Read data from F file
			infF[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minF[iter])) + sizeof(int), ios::beg);
			infF[iter].read((char*)tmpF, sizeof(storageType) * chunk);
			// Read data from R file
			infR[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minR[iter])) + sizeof(int), ios::beg);
			infR[iter].read((char*)tmpR, sizeof(storageType) * chunk);
			
			// Truncate the data to truncValF and truncValR
			truncData(tmpF, chunk, truncValF);
			truncData(tmpR, chunk, truncValR);

			// loop over the read data and calculate sums over the window size. 
			wsumF = 0;
			wsumR = 0;
			wFsumF = 0;
			wFsumR = 0;
			
			wMeanF = 0.0;
			wMeanR = 0.0;
			
			//for (int i = thetaF1-w+1; i < chunk-thetaR2; i++)
			for (int i = 1; i < chunk-thetaR0; i++)
			{
				//wMeanF = ((wMeanF * (double)(wsumF)) + ((double)(i-thetaF2) * (double)tmpF[i-thetaF2]));
				//wMeanR = ((wMeanR * (double)(wsumR)) + ((double)(i+w+thetaR1-1) * (double)tmpR[i+w+thetaR1-1]));

				//cerr << (int)wMeanF << "\t" << (int)wMeanR << "\t" << wFsumF << "\t" << wFsumR << "\t";

				wMeanF = ((wMeanF * (double)(wFsumF)) + ((double)(i-1) * (double)tmpF[i-1]));
				wMeanR = ((wMeanR * (double)(wFsumR)) + ((double)(i+thetaR0) * (double)tmpR[i+thetaR0]));

				//cerr << "added " << i-1 << " and " << i+thetaR0 << endl;
				
				wFsumF += (long)tmpF[i-1];       // window sum
				wFsumR += (long)tmpR[i+thetaR0]; // window sum
				
				if (i >= (thetaF0-w+1))
				{
					wsumF += (long)tmpF[i-thetaF2];      // window sum
					wsumR += (long)tmpR[i+w+thetaR1-1];  // window sum
				}
				
				if (wFsumF != 0)
					wMeanF /= (double)wFsumF;
				else
					wMeanF = 0.0;
				if (wFsumR != 0)
					wMeanR /= (double)wFsumR;
				else
					wMeanR = 0.0;

				//cerr << (int)wMeanF << "\t" << (int)wMeanR << "\t" << wFsumF << "\t" << wFsumR << "\t";
				
				/*
				if (wsumF != 0)
					wMeanF /= (double)wsumF;
				else
					wMeanF = 0.0;
				if (wsumR != 0)
					wMeanR /= (double)wsumR;
				else
					wMeanR = 0.0;
				*/
				
				//if (i >= thetaF1)
				if (i >= thetaF0)
				{
					*maxValF = max((int)wsumF, (*maxValF));
					*maxValR = max((int)wsumR, (*maxValR));

					//if ((wsumF > 0) || (wsumR > 0))
					if ((wFsumF > 0) || (wFsumR > 0))
					{
						wVarF = 0.0;
						wVarR = 0.0;
						wCntF = 0;
						wCntR = 0;

						/*
						for (int j = i - w + 1; j <= i; j++)
						{
							for (int k = 0; k < (int)tmpF[j-thetaF2]; k++)
							{
								wVarF += pow((double)(j-thetaF2) - wMeanF, 2);
								wCntF++;
							}
							for (int k = 0; k < (int)tmpR[j+w+thetaR1-1]; k++)
							{
								wVarR += pow((double)(j+w+thetaR1-1) - wMeanR, 2);
								wCntR++;
							}
						}
						*/

						for (int j = i - thetaF0; j < i; j++)
						{
							for (int k = 0; k < (int)tmpF[j]; k++)
							{
								wVarF += pow((double)j - wMeanF, 2);
								wCntF++;
							}
							for (int k = 0; k < (int)tmpR[j+thetaF0+1]; k++)
							{
								wVarR += pow((double)(j+thetaF0+1) - wMeanR, 2);
								wCntR++;
							}
						}
						
						if (wCntF != 0)
							wVarF /= (double)wCntF;
						else
							wVarF = 0.0;
						if (wCntR != 0)
							wVarR /= (double)wCntR;
						else
							wVarR = 0.0;

						sd = (sqrt(wVarF) + sqrt(wVarR)) * 0.5;
						
						if (sd > maxSDFR)
							maxSDFR = sd;
					}
					
					//wMeanF = (wMeanF * (double)wsumF) - ((double)tmpF[i-thetaF1] * (double)(i-thetaF1));
					//wMeanR = (wMeanR * (double)wsumR) - ((double)tmpR[i+thetaR1] * (double)(i+thetaR1));

					wMeanF = (wMeanF * (double)wFsumF) - ((double)tmpF[i-thetaF0] * (double)(i-thetaF0));
					wMeanR = (wMeanR * (double)wFsumR) - ((double)tmpR[i+1] * (double)(i+1));

					// decrease the window sum with the value just before the current window.
					wsumF -= (long)tmpF[i-thetaF1];
					// decrease the window sum with the value just before the current window.
					wsumR -= (long)tmpR[i+thetaR1];

					// decrease the window sum with the value just before the current window.
					wFsumF -= (long)tmpF[i-thetaF0];
					// decrease the window sum with the value just before the current window.
					wFsumR -= (long)tmpR[i+1];

					//cerr << "removed " << i-thetaF0 << " and " << i+1 << endl;
					
					/*
					if (wsumF != 0)
						wMeanF /= (double)wsumF;
					else
						wMeanF = 0.0;
					
					if (wsumR != 0)
						wMeanR /= (double)wsumR;
					else
						wMeanR = 0.0;
					*/
					if (wFsumF != 0)
						wMeanF /= (double)wFsumF;
					else
						wMeanF = 0.0;
					
					if (wFsumR != 0)
						wMeanR /= (double)wFsumR;
					else
						wMeanR = 0.0;

					//cerr << (int)wMeanF << "\t" << (int)wMeanR << "\t" << wFsumF << "\t" << wFsumR;
				}
				//cerr << endl;
			}
			//read += chunk-thetaF1-thetaR2;
			read += chunk-thetaF0-thetaR0;
			
			// Check if we can read in a full chunk of data or simply the rest of the file. 
			//if (((maxPos[iter] - minPos[iter] + 1 - read) < chunk) && ((maxPos[iter] - minPos[iter] + 1 - read) >= maxD))
			if (((maxPos[iter] - minPos[iter] + 1 - read) < chunk) && ((maxPos[iter] - minPos[iter] + 1 - read) >= (thetaF0+thetaR0+1)))
				chunk = (maxPos[iter] - minPos[iter] + 1) - read;
		}
		delete[] tmpF;
		delete[] tmpR;
		
		totalRead += (double)read / 1000.0;
		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * totalRead / totalLength << "  \t% complete.\r";

		// Check if we have read past the EOF so the eofbit is set
		if (infF[iter].eof())
			infF[iter].clear();
		if (infR[iter].eof())
			infR[iter].clear();
	}
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;

	(*outs) << "Max Fuzzy score 0.5 * (SD F + SD R) = " << (int)maxSDFR << endl;
	(*outs) << "Max F window sum = " << *maxValF << endl;
	(*outs) << "Max R window sum = " << *maxValR << endl;

	vcerr(3) << "\t\tMax Fuzzy score 0.5 * (SD F + SD R) = " << (int)maxSDFR << endl;
	vcerr(3) << "\t\tMax F window sum = " << *maxValF << endl;
	vcerr(3) << "\t\tMax R window sum = " << *maxValR << endl;
	
	double oddsFs[*maxValF];
	double oddsRs[*maxValR];

	double maxOddsF = 0.0;
	for (int i = 0; i < *maxValF; i++)
	{
		oddsFs[i] = log(max(gsl_ran_poisson_pdf(i, lambdaF[0]) * pF[0], DBL_MIN)) -
			log(max(gsl_ran_poisson_pdf(i, lambdaF[1]) * pF[1], DBL_MIN));
		maxOddsF = (oddsFs[i] > maxOddsF ? oddsFs[i] : maxOddsF);
	}
	double maxOddsR = 0.0;
	for (int i = 0; i < *maxValR; i++)
	{
		oddsRs[i] = log(max(gsl_ran_poisson_pdf(i, lambdaR[0]) * pR[0], DBL_MIN)) -
			log(max(gsl_ran_poisson_pdf(i, lambdaR[1]) * pR[1], DBL_MIN));
		maxOddsR = (oddsRs[i] > maxOddsR ? oddsRs[i] : maxOddsR);
	}

	double maxOdds = maxOddsF + maxOddsR;

	(*outs) << "Max center position odds (odds start + odds end) = " << maxOdds;
	(*outs) << " (" << maxOddsF << " + " << maxOddsR << ")." << endl; // normalized to the [0,100] interval." << endl;

	vcerr(3) << "\t\tMax center position odds (odds start + odds end) = " << maxOdds;
	vcerr(3) << " (" << maxOddsF << " + " << maxOddsR << ")." << endl; // normalized to the [0,100] interval." << endl;
	
	int pos = thetaF1-1;
	int written = thetaF1;
	double tmp;
	bool newchunk = false;
	
	vcerr(2) << "\t* Calculating odds of center positions" << endl;	
	for (int iter = 0; iter < fCnt; iter++)
	{
		vcerr(3) << "\t  " << fileNames[iter] << endl;
		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";
		
		read = 0;
		chunk = MAXNR;
		chunk = ((maxPos[iter] - minPos[iter] + 1) < chunk ? (maxPos[iter] - minPos[iter] + 1) : chunk);

		tmpF = new storageType[chunk];
		tmpR = new storageType[chunk];
		
		//storageType discreteOdds[chunk-thetaF1-thetaR2];
		//storageType discreteFuzzy[chunk-thetaF1-thetaR2];

		storageType discreteOdds[chunk-thetaF0-thetaR0];
		storageType discreteOddsStart[chunk-thetaF0-thetaR0];
		storageType discreteOddsEnd[chunk-thetaF0-thetaR0];
		storageType discreteFuzzy[chunk-thetaF0-thetaR0];
		
		// Write min coordinate to odds files:
		outo[iter].write((char*)&minPos[iter], sizeof(int));
		// Write min coordinate to fuzzy files:
		outf[iter].write((char*)&minPos[iter], sizeof(int));

		if (bodds)
		{
			outos[iter].write((char*)&minPos[iter], sizeof(int));
			outoe[iter].write((char*)&minPos[iter], sizeof(int));
		}
		
		// Fill with zeros up to thetaF1 - 1
		storageType zero = 0;
		//for (int i = 0; i < thetaF1; i++)
		for (int i = 0; i < thetaF0; i++)
		{
			outo[iter].write((char*)&zero, sizeof(storageType));
			outf[iter].write((char*)&zero, sizeof(storageType));
			if (bodds)
			{
				outos[iter].write((char*)&zero, sizeof(storageType));
				outoe[iter].write((char*)&zero, sizeof(storageType));
			}
		}
		
		pos = thetaF1-1;
		newchunk = false;
		written = thetaF1;
		
		//while (((maxPos[iter] - minPos[iter] + 1) - read >= chunk) && (chunk - maxD >= 0))
		while (((maxPos[iter] - minPos[iter] + 1) - read >= chunk) && (chunk - (thetaF0+thetaR0+1) >= 0))
		{
			// Read data from F file
			infF[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minF[iter])) + sizeof(int), ios::beg);
			infF[iter].read((char*)tmpF, sizeof(storageType) * chunk);
			// Read data from R file
			infR[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minR[iter])) + sizeof(int), ios::beg);
			infR[iter].read((char*)tmpR, sizeof(storageType) * chunk);
			
			// Truncate the data to truncValF and truncValR
			truncData(tmpF, chunk, truncValF);
			truncData(tmpR, chunk, truncValR);

			// loop over the read data and calculate sums over the window size. 
			wsumF = 0;
			wsumR = 0;
			wFsumF = 0;
			wFsumR = 0;

			wMeanF = 0.0;
			wMeanR = 0.0;

			//for (int i = thetaF1-w+1; i < chunk-thetaR2; i++)
			for (int i = 1; i < chunk-thetaR0; i++)
			{
				//wMeanF = ((wMeanF * (double)(wsumF)) + ((double)(i-thetaF2) * (double)tmpF[i-thetaF2]));
				//wMeanR = ((wMeanR * (double)(wsumR)) + ((double)(i+w+thetaR1-1) * (double)tmpR[i+w+thetaR1-1]));

				wMeanF = ((wMeanF * (double)(wFsumF)) + ((double)(i-1) * (double)tmpF[i-1]));
				wMeanR = ((wMeanR * (double)(wFsumR)) + ((double)(i+thetaR0) * (double)tmpR[i+thetaR0]));

				wFsumF += (long)tmpF[i-1];       // window sum
				wFsumR += (long)tmpR[i+thetaR0]; // window sum

				if (i >= (thetaF0-w+1))
				{
					wsumF += (long)tmpF[i-thetaF2];      // window sum
					wsumR += (long)tmpR[i+w+thetaR1-1];  // window sum
				}

				if (wFsumF != 0)
					wMeanF /= (double)wFsumF;
				else
					wMeanF = 0.0;
				if (wFsumR != 0)
					wMeanR /= (double)wFsumR;
				else
					wMeanR = 0.0;
				
				/*
				if (wsumF != 0)
					wMeanF /= (double)wsumF;
				else
					wMeanF = 0.0;
				if (wsumR != 0)
					wMeanR /= (double)wsumR;
				else
					wMeanR = 0.0;
				*/
				
				// if we've read the window size the counts are valid.
				//if (i >= thetaF1)
				if (i >= thetaF0)
				{
					oddsF = oddsFs[(int)wsumF];
					oddsR = oddsRs[(int)wsumR];

					if (bodds)
					{
						if (oddsF > 0.0)
							//discreteOddsStart[i-thetaF0] = np_round(1.0 + (99.0 * (oddsF / maxOddsF)));
							discreteOddsStart[i-thetaF0] = np_round(1.0+oddsF);
						else
							discreteOddsStart[i-thetaF0] = 0;
						if (oddsR > 0.0)
							//discreteOddsEnd[i-thetaF0] = np_round(1.0 + (99.0 * (oddsR / maxOddsR)));
							discreteOddsEnd[i-thetaF0] = np_round(1.0+oddsR);
						else
							discreteOddsEnd[i-thetaF0] = 0;
					}
					
					if (oddsF > 0.0 && oddsR > 0.0)
					{
						//tmp = (oddsF+oddsR) / maxOdds;
						//discreteOdds[i-thetaF1] = np_round(1.0 + (99.0 * tmp));
						//discreteOdds[i-thetaF0] = np_round(1.0 + (99.0 * tmp));
						discreteOdds[i-thetaF0] = np_round(1.0+oddsF+oddsR);

						//if ((wsumF > 0) || (wsumR > 0))
						if ((wFsumF > 0) || (wFsumR > 0))
						{
							wVarF = 0.0;
							wVarR = 0.0;
							wCntF = 0;
							wCntR = 0;
							
							for (int j = i - thetaF0; j < i; j++)
							{
								for (int k = 0; k < (int)tmpF[j]; k++)
								{
									wVarF += pow((double)j - wMeanF, 2);
									wCntF++;
								}
								for (int k = 0; k < (int)tmpR[j+thetaF0+1]; k++)
								{
									wVarR += pow((double)(j+thetaF0+1) - wMeanR, 2);
									wCntR++;
								}
							}
							
							/*
							  for (int j = i - w + 1; j <= i; j++)
							  {
							  for (int k = 0; k < (int)tmpF[j-thetaF2]; k++)
							  {
							  wVarF += pow((double)(j-thetaF2) - wMeanF, 2);
							  wCntF++;
							  }
							  for (int k = 0; k < (int)tmpR[j+w+thetaR1-1]; k++)
							  {
							  wVarR += pow((double)(j+w+thetaR1-1) - wMeanR, 2);
							  wCntR++;
							  }
							  }
							*/
							
							if (wCntF != 0)
								wVarF /= (double)wCntF;
							else
								wVarF = 0.0;
							if (wCntR != 0)
								wVarR /= (double)wCntR;
							else
								wVarR = 0.0;
							
							sd = (sqrt(wVarF) + sqrt(wVarR));
							//tmp = (wVarF + wVarR) / maxSDFR;
							
							//if ((wVarF > 0.0) || (wVarR > 0.0))
							//	cerr << wVarF << "\t" << wVarR << "\t" << tmp << "\t" << maxSDFR << "\t" << np_round(1.0 + (99.0 * tmp)) << endl;
							
							//if (tmp > 1.0)
							//	cerr << tmp << "\t" << wVarF << "\t" << wVarR << "\t" << maxSDFR << endl;
							//discreteFuzzy[i-thetaF1] = np_round(1.0 + (99.0 * tmp));
							//discreteFuzzy[i-thetaF0] = np_round(1.0 + (99.0 * tmp));
							discreteFuzzy[i-thetaF0] = np_round(sd);
						}
						else
						{
							//discreteFuzzy[i-thetaF1] = 0;
							discreteFuzzy[i-thetaF0] = 0;
						}
					}
					else
					{
						//discreteOdds[i-thetaF1] = 0;
						discreteOdds[i-thetaF0] = 0;
						discreteFuzzy[i-thetaF0] = 0;
					}
					pos++;
					
					//wMeanF = (wMeanF * (double)wsumF) - ((double)tmpF[i-thetaF1] * (double)(i-thetaF1));
					//wMeanR = (wMeanR * (double)wsumR) - ((double)tmpR[i+thetaR1] * (double)(i+thetaR1));

					wMeanF = (wMeanF * (double)wFsumF) - ((double)tmpF[i-thetaF0] * (double)(i-thetaF0));
					wMeanR = (wMeanR * (double)wFsumR) - ((double)tmpR[i+1] * (double)(i+1));
					
					// Update the means
					idF = (oddsF > 0.0 ? 0 : 1);
					meanF[idF] = (meanF[idF] * ((double)cntF[idF] / (double)(cntF[idF] + 1))) +
						(wsumF / (double)(cntF[idF] + 1));
					cntF[idF]++;
					
					idR = (oddsR > 0.0 ? 0 : 1);
					meanR[idR] = (meanR[idR] * ((double)cntR[idR] / (double)(cntR[idR] + 1))) +
						(wsumR / (double)(cntR[idR] + 1));
					cntR[idR]++;
					
					// decrease the window sum with the value just before the current window.
					wsumF -= (long)tmpF[i-thetaF1];
					// decrease the window sum with the value just before the current window.
					wsumR -= (long)tmpR[i+thetaR1];

					// decrease the window sum with the value just before the current window.
					wFsumF -= (long)tmpF[i-thetaF0];
					// decrease the window sum with the value just before the current window.
					wFsumR -= (long)tmpR[i+1];

					/*
					if (wsumF != 0)
						wMeanF /= (double)wsumF;
					else
						wMeanF = 0.0;
					if (wsumR != 0)
						wMeanR /= (double)wsumR;
					else
						wMeanR = 0.0;
					*/
					if (wFsumF != 0)
						wMeanF /= (double)wFsumF;
					else
						wMeanF = 0.0;
					
					if (wFsumR != 0)
						wMeanR /= (double)wFsumR;
					else
						wMeanR = 0.0;
				}
			}
			
			// Write the chunk of calls to file
			if (chunk == MAXNR)
			{
				outo[iter].write((char*)discreteOdds, sizeof(discreteOdds));
				outf[iter].write((char*)discreteFuzzy, sizeof(discreteFuzzy));
				if (bodds)
				{
					outos[iter].write((char*)discreteOddsStart, sizeof(discreteOddsStart));
					outoe[iter].write((char*)discreteOddsEnd, sizeof(discreteOddsEnd));
				}
			}
			else
			{
				//for (int i = 0; i < chunk-thetaF1-thetaR2; i++)
				for (int i = 0; i < chunk-thetaF0-thetaR0; i++)
				{
					outo[iter].write((char*)&discreteOdds[i], sizeof(storageType));
					outf[iter].write((char*)&discreteFuzzy[i], sizeof(storageType));
					if (bodds)
					{
						outos[iter].write((char*)&discreteOddsStart[i], sizeof(storageType));
						outoe[iter].write((char*)&discreteOddsEnd[i], sizeof(storageType));
					}
				}
			}
			//written += chunk-thetaF1-thetaR2;
			written += chunk-thetaF0-thetaR0;

			if (!newchunk)
			{
				//read += chunk-thetaF1-thetaR2;
				read += chunk-thetaF0-thetaR0;
			}
			else
				read += chunk;
			
			// Check if we can read in a full chunk of data or simply the rest of the file. 
			//if (((maxPos[iter] - minPos[iter] + 1 - read) < chunk) && ((maxPos[iter] - minPos[iter] + 1 - read) >= maxD))
			if (((maxPos[iter] - minPos[iter] + 1 - read) < chunk) && ((maxPos[iter] - minPos[iter] + 1 - read) >= (thetaF0+thetaR0+1)))
			{
				chunk = (maxPos[iter] - minPos[iter] + 1) - read;
				newchunk = true;
			}

			vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)read / (double)(maxPos[iter] - minPos[iter] + 1) << "  \t% complete.\r";
		}

		delete[] tmpF;
		delete[] tmpR;
		
		// Fill odds file with ending zeros
		for (int i = pos+1; i < (maxPos[iter]-minPos[iter]+1); i++)
		{
			outo[iter].write((char*)&zero, sizeof(storageType));
			outf[iter].write((char*)&zero, sizeof(storageType));
			if (bodds)
			{
				outos[iter].write((char*)&zero, sizeof(storageType));
				outoe[iter].write((char*)&zero, sizeof(storageType));
			}
		}

		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;
		
		// Check if we have read past the EOF so the eofbit is set
		if (infF[iter].eof())
			infF[iter].clear();
		if (infR[iter].eof())
			infR[iter].clear();
	}

	// Calculate variance and covariance
	double varF[2] = {0.0, 0.0};
	double varR[2] = {0.0, 0.0};
	double covFR[2] = {0.0, 0.0};
	long cntFR[2] = {0,0};
	int idFR;

	totalRead = 0.0;	
	
	vcerr(2) << "\t* Calculating statistics" << endl;
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";

	cerr << meanF[0] << " " << meanF[1] << endl;
	cerr << cntF[0] << " " << cntF[1] << endl;
	
	cntF[0] = 0;
	cntF[1] = 0;
	cntR[0] = 0;
	cntR[1] = 0;

	for (int iter = 0; iter < fCnt; iter++)
	{
		read = 0;
		chunk = CHUNK_MULTI*MAXNR;
		chunk = ((maxPos[iter] - minPos[iter] + 1) < chunk ? (maxPos[iter] - minPos[iter] + 1) : chunk);

		tmpF = new storageType[chunk];
		tmpR = new storageType[chunk];
	
		while (((maxPos[iter] - minPos[iter] + 1) - read >= chunk) && (chunk - maxD >= 0))
		{
			// Read data from F file
			infF[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minF[iter])) + sizeof(int), ios::beg);
			infF[iter].read((char*)tmpF, sizeof(storageType) * chunk);
			// Read data from R file
			infR[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minR[iter])) + sizeof(int), ios::beg);
			infR[iter].read((char*)tmpR, sizeof(storageType) * chunk);

			// Truncate the data to truncValF and truncValR
			truncData(tmpF, chunk, truncValF);
			truncData(tmpR, chunk, truncValR);
		
			// loop over the read data and calculate sums over the window size. 
			wsumF = 0;
			wsumR = 0;
			for (int i = thetaF1-w+1; i < chunk-thetaR2; i++)
			{
				wsumF += (long)tmpF[i-thetaF2];      // window sum
				wsumR += (long)tmpR[i+w+thetaR1-1];  // window sum
				
				// if we've read the window size the counts are valid.
				if (i >= thetaF1)
				{
					oddsF = oddsFs[(int)wsumF];
					oddsR = oddsRs[(int)wsumR];
					
					// Update the variances and covariances
					idF = (oddsF > 0.0 ? 0 : 1);
					idR = (oddsR > 0.0 ? 0 : 1);
					idFR = (idF == 0 && idR == 0 ? 0 : 1);

					varF[idF] = (varF[idF] * ((double)cntF[idF] / (double)(cntF[idF] + 1))) +
						(pow((double)wsumF - meanF[idF], 2) / (double)(cntF[idF] + 1));
					varR[idR] = (varR[idR] * ((double)cntR[idR] / (double)(cntR[idR] + 1))) +
						(pow((double)wsumR - meanR[idR], 2) / (double)(cntR[idR] + 1));
					covFR[idFR] = (covFR[idFR] * ((double)cntFR[idFR] / (double)(cntFR[idFR] + 1))) +
						(((double)wsumF - meanF[idF]) * ((double)wsumR - meanR[idR]) / (double)(cntFR[idFR] + 1));
				
					cntF[idF]++;
					cntR[idR]++;
					cntFR[idFR]++;
					
					// decrease the window sum with the value just before the current window.
					wsumF -= (long)tmpF[i-thetaF1];
					// decrease the window sum with the value just before the current window.
					wsumR -= (long)tmpR[i+thetaR1];
				}
			}

			read += chunk-thetaF1-thetaR2;
			
			// Check if we can read in a full chunk of data or simply the rest of the file. 
			if (((maxPos[iter] - minPos[iter] + 1 - read) < chunk) && ((maxPos[iter] - minPos[iter] + 1 - read) >= maxD))
				chunk = (maxPos[iter] - minPos[iter] + 1) - read;
		}
		delete[] tmpF;
		delete[] tmpR;

		totalRead += (double)read / 1000.0;
		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * totalRead / totalLength << "  \t% complete.\r";
		
		// Check if we have read past the EOF so the eofbit is set
		if (infF[iter].eof())
			infF[iter].clear();
		if (infR[iter].eof())
			infR[iter].clear();
	}
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;

	vcerr(3) << "\t\tStatistics:" << endl;
	vcerr(3) << "\t\tF: mean = (" <<  meanF[0] << ", " << meanF[1] << ")" << endl;
	vcerr(3) << "\t\tR: mean = (" <<  meanR[0] << ", " << meanR[1] << ")" << endl;
	vcerr(3) << "\t\tF: variance = (" <<  varF[0] << ", " << varF[1] << ")" << endl;
	vcerr(3) << "\t\tR: variance = (" <<  varR[0] << ", " << varR[1] << ")" << endl;
	vcerr(3) << "\t\tCovariance = (" <<  covFR[0] << ", " << covFR[1] << ")" << endl;

	(*outs) << "Statistics:" << endl;
	(*outs) << "F: mean = (" <<  meanF[0] << ", " << meanF[1] << ")" << endl;
	(*outs) << "R: mean = (" <<  meanR[0] << ", " << meanR[1] << ")" << endl;
	(*outs) << "F: variance = (" <<  varF[0] << ", " << varF[1] << ")" << endl;
	(*outs) << "R: variance = (" <<  varR[0] << ", " << varR[1] << ")" << endl;
	(*outs) << "Covariance = (" <<  covFR[0] << ", " << covFR[1] << ")" << endl;
	
	return(1);
}

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
						  int fCnt, string* fileNames)
{	
	/*
	 * 1. Assign the data range variables
	 */

	// Window size of reads
	int w = ceil((double)(maxD - minD) / 2.0);
	// Offsets from center position
	int thetaF1, thetaF2, thetaR1, thetaR2;
	thetaF1 = thetaR2 = ceil((double)(maxD - 1) / 2.0);
	thetaF2 = thetaR1 = thetaF1 - w + 1;
	maxD = thetaF1 + thetaR2 + 1;

	int thetaF0, thetaR0;
	thetaF0 = thetaR0 = thetaF1 + thetaF2;
	
	/*
	 * 2. Slide the files and calculate odds of center positions
	 */
	
	double meanF[2] = {0.0, 0.0};
	double meanR[2] = {0.0, 0.0};

	long cntF[2] = {0,0};
	long cntR[2] = {0,0};
		
	// current window sums.
	long wsumF = 0, wsumR = 0; 
	// current window sums.
	long wFsumF = 0, wFsumR = 0; 

	storageType* tmpF;
	storageType* tmpR;

	int read = 0;
	double oddsF;
	double oddsR;
	int idF;
	int idR;

	double maxSDFR = 0.0;
	double wMeanF = 0.0;
	double wMeanR = 0.0;
	double wVarF = 0.0;
	double wVarR = 0.0;
	int wCntF = 0;
	int wCntR = 0;

	double sd;
	
	*maxValF = 0;
	*maxValR = 0;

	int chunk;
	
	double totalLength = 0;
	for (int i = 0; i < fCnt; i++)
		totalLength += (maxPos[i] - minPos[i] + 1) / 1000.0;
	
	double totalRead = 0.0;
	
	vcerr(2) << "\t* Looking up max value" << endl;
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";

	for (int iter = 0; iter < fCnt; iter++)
	{
		chunk = CHUNK_MULTI*MAXNR;
		// if necessary, reduce to file size.
		chunk = ((maxPos[iter] - minPos[iter] + 1) < chunk ? (maxPos[iter] - minPos[iter] + 1) : chunk);
		
		tmpF = new storageType[chunk];
		tmpR = new storageType[chunk];

		read = 0;
		
		//while (((maxPos[iter] - minPos[iter] + 1) - read >= chunk) && (chunk - maxD >= 0))
		while (((maxPos[iter] - minPos[iter] + 1) - read >= chunk) && (chunk - (thetaF0+thetaR0+1) >= 0))
		{
			// Read data from F file
			infF[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minF[iter])) + sizeof(int), ios::beg);
			infF[iter].read((char*)tmpF, sizeof(storageType) * chunk);
			// Read data from R file
			infR[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minR[iter])) + sizeof(int), ios::beg);
			infR[iter].read((char*)tmpR, sizeof(storageType) * chunk);
			
			// Truncate the data to truncValF and truncValR
			truncData(tmpF, chunk, truncValF);
			truncData(tmpR, chunk, truncValR);

			// loop over the read data and calculate sums over the window size. 
			wsumF = 0;
			wsumR = 0;
			wFsumF = 0;
			wFsumR = 0;
			
			wMeanF = 0.0;
			wMeanR = 0.0;
			
			//for (int i = thetaF1-w+1; i < chunk-thetaR2; i++)
			for (int i = 1; i < chunk-thetaR0; i++)
			{
				//wMeanF = ((wMeanF * (double)(wsumF)) + ((double)(i-thetaF2) * (double)tmpF[i-thetaF2]));
				//wMeanR = ((wMeanR * (double)(wsumR)) + ((double)(i+w+thetaR1-1) * (double)tmpR[i+w+thetaR1-1]));

				//cerr << (int)wMeanF << "\t" << (int)wMeanR << "\t" << wFsumF << "\t" << wFsumR << "\t";

				wMeanF = ((wMeanF * (double)(wFsumF)) + ((double)(i-1) * (double)tmpF[i-1]));
				wMeanR = ((wMeanR * (double)(wFsumR)) + ((double)(i+thetaR0) * (double)tmpR[i+thetaR0]));

				//cerr << "added " << i-1 << " and " << i+thetaR0 << endl;
				
				wFsumF += (long)tmpF[i-1];       // window sum
				wFsumR += (long)tmpR[i+thetaR0]; // window sum
				
				if (i >= (thetaF0-w+1))
				{
					wsumF += (long)tmpF[i-thetaF2];      // window sum
					wsumR += (long)tmpR[i+w+thetaR1-1];  // window sum
				}
				
				if (wFsumF != 0)
					wMeanF /= (double)wFsumF;
				else
					wMeanF = 0.0;
				if (wFsumR != 0)
					wMeanR /= (double)wFsumR;
				else
					wMeanR = 0.0;

				//cerr << (int)wMeanF << "\t" << (int)wMeanR << "\t" << wFsumF << "\t" << wFsumR << "\t";
				
				/*
				if (wsumF != 0)
					wMeanF /= (double)wsumF;
				else
					wMeanF = 0.0;
				if (wsumR != 0)
					wMeanR /= (double)wsumR;
				else
					wMeanR = 0.0;
				*/
				
				//if (i >= thetaF1)
				if (i >= thetaF0)
				{
					*maxValF = max((int)wsumF, (*maxValF));
					*maxValR = max((int)wsumR, (*maxValR));

					//if ((wsumF > 0) || (wsumR > 0))
					if ((wFsumF > 0) || (wFsumR > 0))
					{
						wVarF = 0.0;
						wVarR = 0.0;
						wCntF = 0;
						wCntR = 0;

						/*
						for (int j = i - w + 1; j <= i; j++)
						{
							for (int k = 0; k < (int)tmpF[j-thetaF2]; k++)
							{
								wVarF += pow((double)(j-thetaF2) - wMeanF, 2);
								wCntF++;
							}
							for (int k = 0; k < (int)tmpR[j+w+thetaR1-1]; k++)
							{
								wVarR += pow((double)(j+w+thetaR1-1) - wMeanR, 2);
								wCntR++;
							}
						}
						*/

						for (int j = i - thetaF0; j < i; j++)
						{
							for (int k = 0; k < (int)tmpF[j]; k++)
							{
								wVarF += pow((double)j - wMeanF, 2);
								wCntF++;
							}
							for (int k = 0; k < (int)tmpR[j+thetaF0+1]; k++)
							{
								wVarR += pow((double)(j+thetaF0+1) - wMeanR, 2);
								wCntR++;
							}
						}
						
						if (wCntF != 0)
							wVarF /= (double)wCntF;
						else
							wVarF = 0.0;
						if (wCntR != 0)
							wVarR /= (double)wCntR;
						else
							wVarR = 0.0;

						sd = (sqrt(wVarF) + sqrt(wVarR)) * 0.5;
						
						if (sd > maxSDFR)
							maxSDFR = sd;
					}
					
					//wMeanF = (wMeanF * (double)wsumF) - ((double)tmpF[i-thetaF1] * (double)(i-thetaF1));
					//wMeanR = (wMeanR * (double)wsumR) - ((double)tmpR[i+thetaR1] * (double)(i+thetaR1));

					wMeanF = (wMeanF * (double)wFsumF) - ((double)tmpF[i-thetaF0] * (double)(i-thetaF0));
					wMeanR = (wMeanR * (double)wFsumR) - ((double)tmpR[i+1] * (double)(i+1));

					// decrease the window sum with the value just before the current window.
					wsumF -= (long)tmpF[i-thetaF1];
					// decrease the window sum with the value just before the current window.
					wsumR -= (long)tmpR[i+thetaR1];

					// decrease the window sum with the value just before the current window.
					wFsumF -= (long)tmpF[i-thetaF0];
					// decrease the window sum with the value just before the current window.
					wFsumR -= (long)tmpR[i+1];

					//cerr << "removed " << i-thetaF0 << " and " << i+1 << endl;
					
					/*
					if (wsumF != 0)
						wMeanF /= (double)wsumF;
					else
						wMeanF = 0.0;
					
					if (wsumR != 0)
						wMeanR /= (double)wsumR;
					else
						wMeanR = 0.0;
					*/
					if (wFsumF != 0)
						wMeanF /= (double)wFsumF;
					else
						wMeanF = 0.0;
					
					if (wFsumR != 0)
						wMeanR /= (double)wFsumR;
					else
						wMeanR = 0.0;

					//cerr << (int)wMeanF << "\t" << (int)wMeanR << "\t" << wFsumF << "\t" << wFsumR;
				}
				//cerr << endl;
			}
			//read += chunk-thetaF1-thetaR2;
			read += chunk-thetaF0-thetaR0;
			
			// Check if we can read in a full chunk of data or simply the rest of the file. 
			//if (((maxPos[iter] - minPos[iter] + 1 - read) < chunk) && ((maxPos[iter] - minPos[iter] + 1 - read) >= maxD))
			if (((maxPos[iter] - minPos[iter] + 1 - read) < chunk) && ((maxPos[iter] - minPos[iter] + 1 - read) >= (thetaF0+thetaR0+1)))
				chunk = (maxPos[iter] - minPos[iter] + 1) - read;
		}
		delete[] tmpF;
		delete[] tmpR;
		
		totalRead += (double)read / 1000.0;
		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * totalRead / totalLength << "  \t% complete.\r";

		// Check if we have read past the EOF so the eofbit is set
		if (infF[iter].eof())
			infF[iter].clear();
		if (infR[iter].eof())
			infR[iter].clear();
	}
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;

	(*outs) << "Max Fuzzy score 0.5 * (SD F + SD R) = " << (int)maxSDFR << endl;
	(*outs) << "Max F window sum = " << *maxValF << endl;
	(*outs) << "Max R window sum = " << *maxValR << endl;

	vcerr(3) << "\t\tMax Fuzzy score 0.5 * (SD F + SD R) = " << (int)maxSDFR << endl;
	vcerr(3) << "\t\tMax F window sum = " << *maxValF << endl;
	vcerr(3) << "\t\tMax R window sum = " << *maxValR << endl;
	
	double oddsFs[*maxValF];
	double oddsRs[*maxValR];

	double maxOddsF = 0.0;
	for (int i = 0; i < *maxValF; i++)
	{
		oddsFs[i] = log(max(gsl_ran_poisson_pdf(i, lambdaF[0]) * pF[0], DBL_MIN)) -
			log(max(gsl_ran_poisson_pdf(i, lambdaF[1]) * pF[1], DBL_MIN));
		maxOddsF = (oddsFs[i] > maxOddsF ? oddsFs[i] : maxOddsF);
	}
	double maxOddsR = 0.0;
	for (int i = 0; i < *maxValR; i++)
	{
		oddsRs[i] = log(max(gsl_ran_poisson_pdf(i, lambdaR[0]) * pR[0], DBL_MIN)) -
			log(max(gsl_ran_poisson_pdf(i, lambdaR[1]) * pR[1], DBL_MIN));
		maxOddsR = (oddsRs[i] > maxOddsR ? oddsRs[i] : maxOddsR);
	}

	double maxOdds = maxOddsF + maxOddsR;

	(*outs) << "Max center position odds (odds start + odds end) = " << maxOdds;
	(*outs) << " (" << maxOddsF << " + " << maxOddsR << ")." << endl; // normalized to the [0,100] interval." << endl;

	vcerr(3) << "\t\tMax center position odds (odds start + odds end) = " << maxOdds;
	vcerr(3) << " (" << maxOddsF << " + " << maxOddsR << ")." << endl; // normalized to the [0,100] interval." << endl;
	
	int pos = thetaF1-1;
	int written = thetaF1;
	double tmp;
	bool newchunk = false;
	
	vcerr(2) << "\t* Calculating odds of center positions" << endl;	
	for (int iter = 0; iter < fCnt; iter++)
	{
		vcerr(3) << "\t  " << fileNames[iter] << endl;
		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";
		
		read = 0;
		chunk = MAXNR;
		chunk = ((maxPos[iter] - minPos[iter] + 1) < chunk ? (maxPos[iter] - minPos[iter] + 1) : chunk);

		tmpF = new storageType[chunk];
		tmpR = new storageType[chunk];
		
		//storageType discreteOdds[chunk-thetaF1-thetaR2];
		//storageType discreteFuzzy[chunk-thetaF1-thetaR2];

		short discreteOdds[chunk-thetaF0-thetaR0];
		short discreteOddsStart[chunk-thetaF0-thetaR0];
		short discreteOddsEnd[chunk-thetaF0-thetaR0];
		storageType discreteFuzzy[chunk-thetaF0-thetaR0];
		
		// Write min coordinate to odds files:
		outo[iter].write((char*)&minPos[iter], sizeof(int));
		// Write min coordinate to fuzzy files:
		outf[iter].write((char*)&minPos[iter], sizeof(int));

		if (bodds)
		{
			outos[iter].write((char*)&minPos[iter], sizeof(int));
			outoe[iter].write((char*)&minPos[iter], sizeof(int));
		}
		
		// Fill with zeros up to thetaF1 - 1
		storageType zero = 0;
		//for (int i = 0; i < thetaF1; i++)
		for (int i = 0; i < thetaF0; i++)
		{
			outo[iter].write((char*)&zero, sizeof(short));
			outf[iter].write((char*)&zero, sizeof(storageType));
			if (bodds)
			{
				outos[iter].write((char*)&zero, sizeof(short));
				outoe[iter].write((char*)&zero, sizeof(short));
			}
		}
		
		pos = thetaF1-1;
		newchunk = false;
		written = thetaF1;
		
		//while (((maxPos[iter] - minPos[iter] + 1) - read >= chunk) && (chunk - maxD >= 0))
		while (((maxPos[iter] - minPos[iter] + 1) - read >= chunk) && (chunk - (thetaF0+thetaR0+1) >= 0))
		{
			// Read data from F file
			infF[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minF[iter])) + sizeof(int), ios::beg);
			infF[iter].read((char*)tmpF, sizeof(storageType) * chunk);
			// Read data from R file
			infR[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minR[iter])) + sizeof(int), ios::beg);
			infR[iter].read((char*)tmpR, sizeof(storageType) * chunk);
			
			// Truncate the data to truncValF and truncValR
			truncData(tmpF, chunk, truncValF);
			truncData(tmpR, chunk, truncValR);

			// loop over the read data and calculate sums over the window size. 
			wsumF = 0;
			wsumR = 0;
			wFsumF = 0;
			wFsumR = 0;

			wMeanF = 0.0;
			wMeanR = 0.0;

			//for (int i = thetaF1-w+1; i < chunk-thetaR2; i++)
			for (int i = 1; i < chunk-thetaR0; i++)
			{
				//wMeanF = ((wMeanF * (double)(wsumF)) + ((double)(i-thetaF2) * (double)tmpF[i-thetaF2]));
				//wMeanR = ((wMeanR * (double)(wsumR)) + ((double)(i+w+thetaR1-1) * (double)tmpR[i+w+thetaR1-1]));

				wMeanF = ((wMeanF * (double)(wFsumF)) + ((double)(i-1) * (double)tmpF[i-1]));
				wMeanR = ((wMeanR * (double)(wFsumR)) + ((double)(i+thetaR0) * (double)tmpR[i+thetaR0]));

				wFsumF += (long)tmpF[i-1];       // window sum
				wFsumR += (long)tmpR[i+thetaR0]; // window sum

				if (i >= (thetaF0-w+1))
				{
					wsumF += (long)tmpF[i-thetaF2];      // window sum
					wsumR += (long)tmpR[i+w+thetaR1-1];  // window sum
				}

				if (wFsumF != 0)
					wMeanF /= (double)wFsumF;
				else
					wMeanF = 0.0;
				if (wFsumR != 0)
					wMeanR /= (double)wFsumR;
				else
					wMeanR = 0.0;
				
				/*
				if (wsumF != 0)
					wMeanF /= (double)wsumF;
				else
					wMeanF = 0.0;
				if (wsumR != 0)
					wMeanR /= (double)wsumR;
				else
					wMeanR = 0.0;
				*/
				
				// if we've read the window size the counts are valid.
				//if (i >= thetaF1)
				if (i >= thetaF0)
				{
					oddsF = oddsFs[(int)wsumF];
					oddsR = oddsRs[(int)wsumR];

					if (bodds)
					{
						//if (oddsF > 0.0)
							//discreteOddsStart[i-thetaF0] = np_round(1.0 + (99.0 * (oddsF / maxOddsF)));
							//discreteOddsStart[i-thetaF0] = np_round(1.0+oddsF);
						discreteOddsStart[i-thetaF0] = np_round(oddsF);
						//else
						//	discreteOddsStart[i-thetaF0] = 0;
						//if (oddsR > 0.0)
							//discreteOddsEnd[i-thetaF0] = np_round(1.0 + (99.0 * (oddsR / maxOddsR)));
						//	discreteOddsEnd[i-thetaF0] = np_round(1.0+oddsR);
						discreteOddsEnd[i-thetaF0] = np_round(oddsR);
						//else
						//	discreteOddsEnd[i-thetaF0] = 0;
					}

					discreteOdds[i-thetaF0] = np_round(oddsF+oddsR);
					if (oddsF > 0.0 && oddsR > 0.0)
					{
						//tmp = (oddsF+oddsR) / maxOdds;
						//discreteOdds[i-thetaF1] = np_round(1.0 + (99.0 * tmp));
						//discreteOdds[i-thetaF0] = np_round(1.0 + (99.0 * tmp));
						//discreteOdds[i-thetaF0] = np_round(1.0+oddsF+oddsR);

						//if ((wsumF > 0) || (wsumR > 0))
						if ((wFsumF > 0) || (wFsumR > 0))
						{
							wVarF = 0.0;
							wVarR = 0.0;
							wCntF = 0;
							wCntR = 0;
							
							for (int j = i - thetaF0; j < i; j++)
							{
								for (int k = 0; k < (int)tmpF[j]; k++)
								{
									wVarF += pow((double)j - wMeanF, 2);
									wCntF++;
								}
								for (int k = 0; k < (int)tmpR[j+thetaF0+1]; k++)
								{
									wVarR += pow((double)(j+thetaF0+1) - wMeanR, 2);
									wCntR++;
								}
							}
							
							/*
							  for (int j = i - w + 1; j <= i; j++)
							  {
							  for (int k = 0; k < (int)tmpF[j-thetaF2]; k++)
							  {
							  wVarF += pow((double)(j-thetaF2) - wMeanF, 2);
							  wCntF++;
							  }
							  for (int k = 0; k < (int)tmpR[j+w+thetaR1-1]; k++)
							  {
							  wVarR += pow((double)(j+w+thetaR1-1) - wMeanR, 2);
							  wCntR++;
							  }
							  }
							*/
							
							if (wCntF != 0)
								wVarF /= (double)wCntF;
							else
								wVarF = 0.0;
							if (wCntR != 0)
								wVarR /= (double)wCntR;
							else
								wVarR = 0.0;
							
							sd = (sqrt(wVarF) + sqrt(wVarR));
							//tmp = (wVarF + wVarR) / maxSDFR;
							
							//if ((wVarF > 0.0) || (wVarR > 0.0))
							//	cerr << wVarF << "\t" << wVarR << "\t" << tmp << "\t" << maxSDFR << "\t" << np_round(1.0 + (99.0 * tmp)) << endl;
							
							//if (tmp > 1.0)
							//	cerr << tmp << "\t" << wVarF << "\t" << wVarR << "\t" << maxSDFR << endl;
							//discreteFuzzy[i-thetaF1] = np_round(1.0 + (99.0 * tmp));
							//discreteFuzzy[i-thetaF0] = np_round(1.0 + (99.0 * tmp));
							discreteFuzzy[i-thetaF0] = np_round(sd);
						}
						else
						{
							//discreteFuzzy[i-thetaF1] = 0;
							discreteFuzzy[i-thetaF0] = 0;
						}
					}
					else
					{
						//discreteOdds[i-thetaF1] = 0;
						//discreteOdds[i-thetaF0] = 0;
						discreteFuzzy[i-thetaF0] = 0;
					}
					pos++;
					
					//wMeanF = (wMeanF * (double)wsumF) - ((double)tmpF[i-thetaF1] * (double)(i-thetaF1));
					//wMeanR = (wMeanR * (double)wsumR) - ((double)tmpR[i+thetaR1] * (double)(i+thetaR1));

					wMeanF = (wMeanF * (double)wFsumF) - ((double)tmpF[i-thetaF0] * (double)(i-thetaF0));
					wMeanR = (wMeanR * (double)wFsumR) - ((double)tmpR[i+1] * (double)(i+1));
					
					// Update the means
					idF = (oddsF > 0.0 ? 0 : 1);
					meanF[idF] = (meanF[idF] * ((double)cntF[idF] / (double)(cntF[idF] + 1))) +
						(wsumF / (double)(cntF[idF] + 1));
					cntF[idF]++;
					
					idR = (oddsR > 0.0 ? 0 : 1);
					meanR[idR] = (meanR[idR] * ((double)cntR[idR] / (double)(cntR[idR] + 1))) +
						(wsumR / (double)(cntR[idR] + 1));
					cntR[idR]++;
					
					// decrease the window sum with the value just before the current window.
					wsumF -= (long)tmpF[i-thetaF1];
					// decrease the window sum with the value just before the current window.
					wsumR -= (long)tmpR[i+thetaR1];

					// decrease the window sum with the value just before the current window.
					wFsumF -= (long)tmpF[i-thetaF0];
					// decrease the window sum with the value just before the current window.
					wFsumR -= (long)tmpR[i+1];

					/*
					if (wsumF != 0)
						wMeanF /= (double)wsumF;
					else
						wMeanF = 0.0;
					if (wsumR != 0)
						wMeanR /= (double)wsumR;
					else
						wMeanR = 0.0;
					*/
					if (wFsumF != 0)
						wMeanF /= (double)wFsumF;
					else
						wMeanF = 0.0;
					
					if (wFsumR != 0)
						wMeanR /= (double)wFsumR;
					else
						wMeanR = 0.0;
				}
			}
			
			// Write the chunk of calls to file
			if (chunk == MAXNR)
			{
				outo[iter].write((char*)discreteOdds, sizeof(discreteOdds));
				outf[iter].write((char*)discreteFuzzy, sizeof(discreteFuzzy));
				if (bodds)
				{
					outos[iter].write((char*)discreteOddsStart, sizeof(discreteOddsStart));
					outoe[iter].write((char*)discreteOddsEnd, sizeof(discreteOddsEnd));
				}
			}
			else
			{
				//for (int i = 0; i < chunk-thetaF1-thetaR2; i++)
				for (int i = 0; i < chunk-thetaF0-thetaR0; i++)
				{
					outo[iter].write((char*)&discreteOdds[i], sizeof(short));
					outf[iter].write((char*)&discreteFuzzy[i], sizeof(storageType));
					if (bodds)
					{
						outos[iter].write((char*)&discreteOddsStart[i], sizeof(short));
						outoe[iter].write((char*)&discreteOddsEnd[i], sizeof(short));
					}
				}
			}
			//written += chunk-thetaF1-thetaR2;
			written += chunk-thetaF0-thetaR0;

			if (!newchunk)
			{
				//read += chunk-thetaF1-thetaR2;
				read += chunk-thetaF0-thetaR0;
			}
			else
				read += chunk;
			
			// Check if we can read in a full chunk of data or simply the rest of the file. 
			//if (((maxPos[iter] - minPos[iter] + 1 - read) < chunk) && ((maxPos[iter] - minPos[iter] + 1 - read) >= maxD))
			if (((maxPos[iter] - minPos[iter] + 1 - read) < chunk) && ((maxPos[iter] - minPos[iter] + 1 - read) >= (thetaF0+thetaR0+1)))
			{
				chunk = (maxPos[iter] - minPos[iter] + 1) - read;
				newchunk = true;
			}

			vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)read / (double)(maxPos[iter] - minPos[iter] + 1) << "  \t% complete.\r";
		}

		delete[] tmpF;
		delete[] tmpR;
		
		// Fill odds file with ending zeros
		for (int i = pos+1; i < (maxPos[iter]-minPos[iter]+1); i++)
		{
			outo[iter].write((char*)&zero, sizeof(short));
			outf[iter].write((char*)&zero, sizeof(storageType));
			if (bodds)
			{
				outos[iter].write((char*)&zero, sizeof(short));
				outoe[iter].write((char*)&zero, sizeof(short));
			}
		}

		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;
		
		// Check if we have read past the EOF so the eofbit is set
		if (infF[iter].eof())
			infF[iter].clear();
		if (infR[iter].eof())
			infR[iter].clear();
	}

	// Calculate variance and covariance
	double varF[2] = {0.0, 0.0};
	double varR[2] = {0.0, 0.0};
	double covFR[2] = {0.0, 0.0};
	long cntFR[2] = {0,0};
	int idFR;

	totalRead = 0.0;	
	
	vcerr(2) << "\t* Calculating statistics" << endl;
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";

	cerr << meanF[0] << " " << meanF[1] << endl;
	cerr << cntF[0] << " " << cntF[1] << endl;
	
	cntF[0] = 0;
	cntF[1] = 0;
	cntR[0] = 0;
	cntR[1] = 0;

	for (int iter = 0; iter < fCnt; iter++)
	{
		read = 0;
		chunk = CHUNK_MULTI*MAXNR;
		chunk = ((maxPos[iter] - minPos[iter] + 1) < chunk ? (maxPos[iter] - minPos[iter] + 1) : chunk);

		tmpF = new storageType[chunk];
		tmpR = new storageType[chunk];
	
		while (((maxPos[iter] - minPos[iter] + 1) - read >= chunk) && (chunk - maxD >= 0))
		{
			// Read data from F file
			infF[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minF[iter])) + sizeof(int), ios::beg);
			infF[iter].read((char*)tmpF, sizeof(storageType) * chunk);
			// Read data from R file
			infR[iter].seekg(sizeof(storageType) * (read + (minPos[iter]-minR[iter])) + sizeof(int), ios::beg);
			infR[iter].read((char*)tmpR, sizeof(storageType) * chunk);

			// Truncate the data to truncValF and truncValR
			truncData(tmpF, chunk, truncValF);
			truncData(tmpR, chunk, truncValR);
		
			// loop over the read data and calculate sums over the window size. 
			wsumF = 0;
			wsumR = 0;
			for (int i = thetaF1-w+1; i < chunk-thetaR2; i++)
			{
				wsumF += (long)tmpF[i-thetaF2];      // window sum
				wsumR += (long)tmpR[i+w+thetaR1-1];  // window sum
				
				// if we've read the window size the counts are valid.
				if (i >= thetaF1)
				{
					oddsF = oddsFs[(int)wsumF];
					oddsR = oddsRs[(int)wsumR];
					
					// Update the variances and covariances
					idF = (oddsF > 0.0 ? 0 : 1);
					idR = (oddsR > 0.0 ? 0 : 1);
					idFR = (idF == 0 && idR == 0 ? 0 : 1);

					varF[idF] = (varF[idF] * ((double)cntF[idF] / (double)(cntF[idF] + 1))) +
						(pow((double)wsumF - meanF[idF], 2) / (double)(cntF[idF] + 1));
					varR[idR] = (varR[idR] * ((double)cntR[idR] / (double)(cntR[idR] + 1))) +
						(pow((double)wsumR - meanR[idR], 2) / (double)(cntR[idR] + 1));
					covFR[idFR] = (covFR[idFR] * ((double)cntFR[idFR] / (double)(cntFR[idFR] + 1))) +
						(((double)wsumF - meanF[idF]) * ((double)wsumR - meanR[idR]) / (double)(cntFR[idFR] + 1));
				
					cntF[idF]++;
					cntR[idR]++;
					cntFR[idFR]++;
					
					// decrease the window sum with the value just before the current window.
					wsumF -= (long)tmpF[i-thetaF1];
					// decrease the window sum with the value just before the current window.
					wsumR -= (long)tmpR[i+thetaR1];
				}
			}

			read += chunk-thetaF1-thetaR2;
			
			// Check if we can read in a full chunk of data or simply the rest of the file. 
			if (((maxPos[iter] - minPos[iter] + 1 - read) < chunk) && ((maxPos[iter] - minPos[iter] + 1 - read) >= maxD))
				chunk = (maxPos[iter] - minPos[iter] + 1) - read;
		}
		delete[] tmpF;
		delete[] tmpR;

		totalRead += (double)read / 1000.0;
		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * totalRead / totalLength << "  \t% complete.\r";
		
		// Check if we have read past the EOF so the eofbit is set
		if (infF[iter].eof())
			infF[iter].clear();
		if (infR[iter].eof())
			infR[iter].clear();
	}
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;

	vcerr(3) << "\t\tStatistics:" << endl;
	vcerr(3) << "\t\tF: mean = (" <<  meanF[0] << ", " << meanF[1] << ")" << endl;
	vcerr(3) << "\t\tR: mean = (" <<  meanR[0] << ", " << meanR[1] << ")" << endl;
	vcerr(3) << "\t\tF: variance = (" <<  varF[0] << ", " << varF[1] << ")" << endl;
	vcerr(3) << "\t\tR: variance = (" <<  varR[0] << ", " << varR[1] << ")" << endl;
	vcerr(3) << "\t\tCovariance = (" <<  covFR[0] << ", " << covFR[1] << ")" << endl;

	(*outs) << "Statistics:" << endl;
	(*outs) << "F: mean = (" <<  meanF[0] << ", " << meanF[1] << ")" << endl;
	(*outs) << "R: mean = (" <<  meanR[0] << ", " << meanR[1] << ")" << endl;
	(*outs) << "F: variance = (" <<  varF[0] << ", " << varF[1] << ")" << endl;
	(*outs) << "R: variance = (" <<  varR[0] << ", " << varR[1] << ")" << endl;
	(*outs) << "Covariance = (" <<  covFR[0] << ", " << covFR[1] << ")" << endl;
	
	return(1);
}


