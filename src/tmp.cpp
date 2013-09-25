int checkFileLengths(ifstream* infF, ifstream* infR,
					 int* minPos, int* maxPos,
					 int* minR, int* minF,
					 int* maxF, int* maxR)
{
	vcerr(2) << "\t* Checking file lengths" << endl;
	
	getMinMaxF(infF, minF, maxF);
	getMinMaxF(infR, minR, maxR);

	if ((*minF != *minR) || (*maxF != *maxR)) // Input files have not the same coordinates
    {
 		vcerr(3) << "\t\tFiles start at different positions. Largest common region: ";
		
		*minPos = (*minF > *minR ? *minF : *minR);
		*maxPos = (*maxF < *maxR ? *maxF : *maxR);
		
		vcerr(3) << "[ " << *minPos << ", " << *maxPos << " ]." << endl;
    }
	else // input files have the exact same coordinates
    {
		*minPos = *minF;
		*maxPos = *maxF;
    }

	return(1);
}

int estimateDistance(ifstream* infF, ifstream* infR, ofstream* outf,
					 int* minD, int* maxD, int* lambdaD,
					 int minLag, int maxLag,
					 int minPos, int maxPos,
					 int minR, int minF,
					 int maxF, int maxR,
		     int truncValF, int truncValR, double* corr)
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
	double sdF = 0.0;
	double sdR = 0.0;

	int read = 0;

	//double corr[maxLag-minLag+1];
	double tmpcorr[maxLag-minLag+1];

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

	//minPos = 121150000;
	//maxPos = 121200000;
	
	chunk = CHUNK_MULTI*MAXNR;
	// if necessary, reduce to file size.
	chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);

	read = 0;

	for (int lag = minLag; lag <= maxLag; lag++)
		corr[lag-minLag]  = 0.0;
	
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";
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
		
		for (int lag = minLag; lag <= maxLag; lag++)
			tmpcorr[lag-minLag]  = 0.0;
		
		//Loop all F positions
		for (int i = 0; i <= chunk - maxLag; i++)
		{
			//Traverse from minLag to maxLag
			for (int lag = minLag; lag <= maxLag; lag++)
				tmpcorr[lag-minLag] += ((double)tmpF[i] - meanF) * ((double)tmpR[i+lag-1] - meanR);
			
			sdF += pow((double)tmpF[i] - meanF,2);
			sdR += pow((double)tmpR[i] - meanR,2);
		}

		for (int lag = minLag; lag <= maxLag; lag++)
			corr[lag-minLag] = (corr[lag-minLag] * ((double)read / (double)(read+chunk-maxLag+1))) +
				(tmpcorr[lag-minLag] / (double)(read+chunk-maxLag+1));

		//cerr << endl << corr[47-minLag] << "\t" << corr[150-minLag] << endl;
		//cerr << tmpcorr[47-minLag] / (double)(chunk-maxLag+1) << "\t" << tmpcorr[150-minLag] / (double)(chunk-maxLag+1) << endl;
		
		read += chunk-maxLag+1;
		
		// Check if we can read in a full chunk of data or simply the rest of the file. 
		if (((maxPos - minPos + 1 - read) < chunk) && ((maxPos - minPos + 1 - read) >= maxLag))
			chunk = (maxPos - minPos + 1) - read;
		
		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)read / (double)(maxPos - minPos + 1) << "  \t% complete.\r";
	}
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;

	// Check if we have read past the EOF so the eofbit is set
	if (infF->eof())
		infF->clear();
	if (infR->eof())
		infR->clear();
	
	sdF = sqrt(sdF / (double)(maxPos - minPos + 1 - maxLag));
	sdR = sqrt(sdR / (double)(maxPos - minPos + 1 - maxLag));

	double sum = 0.0;
	for (int lag = minLag; lag <= maxLag; lag++)
	{
		corr[lag-minLag] /= sdF * sdR;
		sum += corr[lag-minLag];
		*(outf) << lag << "\t" << corr[lag-minLag] << endl;
	}

	double maxCorr = 0.0;
	for (int lag = minLag; lag <= maxLag; lag++)
	{
		corr[lag-minLag] /= sum;
		maxCorr = max(maxCorr, corr[lag-minLag]);
	}

	double steps[1000];
	for (int i = 0; i < 1000; i++)
		steps[i] = maxCorr - ((double)(i-1) * maxCorr) / 1000.0;

	double peak = 0.0;
	int peakId = -1;
	double sum = 0.0;
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
			if (corr[lag-minLag] >= val)
			{
				if (sum == 0.0)
				{
					sum = corr[lag-minLag];
					tmpPeak = corr[lag-minLag];
                    tmpPeakId = lag;
                    tmpStart = lag;
                    tmpEnd = lag;
                }
                else
                {
                    sum = sum + corr[lag-minLag];
                    tmpEnd = lag;
                    if (corr[lag-minLag] > tmpPeak)
                    {
                        tmpPeak = corr[lag-minLag];
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
        if (peak != -1)
            break;
    }

	*minD = start;
	*maxD = end;
	*lambdaD = peakId;
	
	delete[] tmpF;
	delete[] tmpR;
	
	return(1);
}
