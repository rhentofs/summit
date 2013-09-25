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
	double sdF = 0.0;
	double sdR = 0.0;

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

	truncValF = 200;
	
	int maxDist = 2*((maxLag-minLag)+np_round((double)minLag/2.0))+1;
	long cnt[maxDist][truncValF];

	for (int dist = 0; dist < maxDist; dist++)
		for (int i = 0; i < truncValF; i++)
			cnt[dist][i] = 0;
	
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";

	bool posPos[maxDist];

	double avg = 0.0;
	long avgcnt = 0;
	
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
							avg = ((avg * (double)avgcnt) + ((double)tmp * (double)tmpF[lowPos+j])) /
								(avgcnt + tmpF[lowPos+j]);
							avgcnt += tmpF[lowPos+j];
							//cnt[tmp] += tmpF[lowPos+j];
							cnt[tmp][tmpF[lowPos+j]] += 1;
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
					avg = ((avg * (double)avgcnt) + ((double)tmp * (double)tmpF[lowPos+j])) /
						(avgcnt + tmpF[lowPos+j]);
					avgcnt += tmpF[lowPos+j];
					//cnt[tmp] += tmpF[lowPos+j];
					cnt[tmp][tmpF[lowPos+j]] += 1;
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
	
	// Check if we have read past the EOF so the eofbit is set
	if (infF->eof())
		infF->clear();
	if (infR->eof())
		infR->clear();
	
	double sum = 0.0;
	//*(outf) << "lag\t" << "cnt\t" << endl;
	
	for (int i = 0; i < maxDist; i++)
	{
		//*(outf) << i << "\t" << cnt[i] << endl;
		*(outf) << i << "\t";
		for (int j = 0; j < truncValF-1; j++)
			*(outf) << cnt[i][j] << "\t";
		*(outf) << cnt[i][truncValF-1] << endl;
	}

	/*
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
