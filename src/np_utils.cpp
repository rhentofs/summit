#include "np_utils.h"

int checkFileLengths(ifstream* infF, ifstream* infR,
		     int* minPos, int* maxPos,
		     int* minR, int* minF,
		     int* maxF, int* maxR)
{
  vcerr(2) << "*** Checking file lengths ***" << endl;
  
  getMinMaxF(infF, minF, maxF);
  getMinMaxF(infR, minR, maxR);
  
  if ((*minF != *minR) || (*maxF != *maxR)) // Input files have not the same coordinates
    {
 		vcerr(3) << "\tFiles start at different positions. Largest common region: ";
		
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

/*
 * Top level function for doing MCMC estimation of parameters
 * infF:        file with forward read counts
 * infR:        file with forward read counts
 * outf:        file where results are stored
 * minD:        min allowed distance between forward and reverse reads to define a nucleosome center
 * maxD:        max allowed distance between forward and reverse reads to define a nucleosome center
 * sampleSize:  number of observations from the window sum data to base inference on
 * burnins:     number of burnin iterations
 * iterations:  number of iterations to base estimation on
 * seed:        use seed for random number generation (-1: base seed on CPU time)
 */

int estimateParameters(ifstream* infF, ifstream* infR, ofstream* outf,
					   int minD, int maxD, int lambdaD,
					   double* pF, double* pR, double* lambdaF, double* lambdaR,
					   int sampleSize, int burnins, int iterations,
					   int seed, int minPos, int maxPos,
					   int minR, int minF,
					   int maxF, int maxR,
					   int truncValF, int truncValR)
{
	/*
	 * 2. Assign the data range variables
	 */

	// Window size of reads
	int w = ceil((double)(maxD - minD) / 2.0);
	// Data range
	int dataRange = maxPos - minPos + 1;
	// Offsets from center position
	int thetaF1, thetaF2, thetaR1, thetaR2;
	thetaF1 = thetaR2 = ceil((double)(maxD - 1) / 2.0);
	thetaF2 = thetaR1 = thetaF1 - w + 1;
	maxD = thetaF1 + thetaR2 + 1;
	// Maximum number of non-overlapping maxDs
	int maxDs = dataRange / (maxD + 1);
	// Valid sample size?
	bool wholeSample = (sampleSize < 1);
	if (wholeSample)
		sampleSize = maxDs;
	bool validSize = (sampleSize <= maxDs);
	if (!validSize)
	{
		sampleSize = maxDs;
		wholeSample = true;
	}

	vcerr(3) << "\t* Settings" << endl;
	vcerr(3) << setprecision(3);
	vcerr(3) << "\t\tRead count window length: " << w << endl;
	vcerr(3) << "\t\tNumber of observations: " << dataRange << endl;
	vcerr(3) << "\t\tMaximum number of non-overlapping regions: " << maxDs << endl;
	if (!validSize)
		vcerr(3) << "\t\tSample size too big, reset to whole sample." << endl;
	vcerr(3) << "\t\tNumber of non-overlapping regions in sample: " << sampleSize;
	vcerr(3) << " " << (wholeSample ? "(whole sample)" : "(user specified)") << endl;
	vcerr(3) << "\t\tOffset forward reads: [" << -thetaF1 << "," << -thetaF2 << "]" << endl;
	vcerr(3) << "\t\tOffset reverse reads: [" << thetaR1 << "," << thetaR2 << "]" << endl;

	/*
	 * 3. Sample the data used for estimation of parameters and calculate sample statistics
	 */
	
	vcerr(2) << "\t* Sampling data" << endl;
	
	// Set seed to random int if not specified (-1)
	if (seed == -1)
	{
		srand(time(NULL));
		seed = rand();
	}

	// Allocate vector of positions on which to base parameter estimation
	int* centers = new int[sampleSize];	
	if (!wholeSample)
	{
		vcerr(3) << "\t\tSampling center positions" << endl;
		np_sample(seed, minPos+thetaF1, maxPos-thetaR2, maxD+1, sampleSize, centers);
	}
	else
	{
		for (int i=0; i < sampleSize; i++)
			centers[i] = i * (maxD+1) + minPos + thetaF1;
	}
	
	vcerr(3) << "\t\tCalculating window sums" << endl;
	int* F = new int[sampleSize];
	int* R = new int[sampleSize];
	double meanF, meanR, varF, varR, covFR;
	int maxValF, maxValR;

	if (windowSums(infF, infR, centers, sampleSize, maxPos, minPos,
				   minF, minR, minD, maxD,
				   &maxValF, &maxValR,
				   thetaF1, thetaF2, thetaR1, thetaR2,
				   F, R, &meanF, &meanR, &varF, &varR, &covFR,
				   truncValF, truncValR) < 0)
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
	 * 4. Estimate parameters through MCMC sampling
	 */

	vcerr(2) << "\t* Estimating parameters" << endl;
		
	// Allocations
	bool ZF[sampleSize];
	bool ZR[sampleSize];

	// Parameters
	double alphaF = pow(meanF,2) / (varF - pow(meanF, 2));
	double alphaR = pow(meanR,2) / (varR - pow(meanR, 2));
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

	// Gibbs sampling
	if (gibbs(F, R, ZF, ZR,
			  pF, pR,
			  lambdaF, lambdaR,
			  deltaF, deltaR,
			  alphaF, alphaR, betaF, betaR,
			  maxValF, maxValR,
			  sampleSize, iterations, burnins, seed,
			  outf) < 0)
	{
		delete[] F;
		delete[] R;

		cerr << endl << "ERROR: gibbs failed, aborting." << endl;
		return(-1);
	}

	delete[] F;
	delete[] R;
	
	return(1);
}

int truncData(storageType* data, int size, int truncVal)
{
	storageType tv = (storageType)truncVal;
	for (int i = 0; i < size; i++)
		data[i] = (data[i] > tv ? tv : data[i]);
	return(1);
}


int calculateTruncVal(ifstream* infF, ifstream* infR, ofstream* outf,
					  int* truncValF, int* truncValR, double truncLimit,
					  int minPos, int maxPos,
					  int minR, int minF,
					  int maxF, int maxR)
{
	int chunk = CHUNK_MULTI*MAXNR;
	// if necessary, reduce to file size.
	chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);

	storageType* tmpF = new storageType[chunk];
	storageType* tmpR = new storageType[chunk];

	int maxValF = 0;
	int maxValR = 0;

	int read = 0;
	
	// Identify max values in F and R files
	vcerr(3) << "\t\t" << "Looking up max values" << endl;
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";
	while ((maxPos - minPos + 1 - read) >= chunk && chunk > 0)
	{
		// Read data from F file
		infF->seekg(sizeof(storageType) * (read + (minPos-minF)) + sizeof(int), ios::beg);
		infF->read((char*)tmpF, sizeof(storageType) * chunk);
		// Read data from R file
		infR->seekg(sizeof(storageType) * (read + (minPos-minR)) + sizeof(int), ios::beg);
		infR->read((char*)tmpR, sizeof(storageType) * chunk);
		
		//Loop all positions
		for (int i = 0; i < chunk; i++)
		{
			maxValF = ((int)tmpF[i] > maxValF ? (int)tmpF[i] : maxValF);
			maxValR = ((int)tmpR[i] > maxValR ? (int)tmpR[i] : maxValR);
		}
		
		read += chunk;
		
		// Check if we can read in a full chunk of data or simply the rest of the file. 
		if (((maxPos - minPos + 1 - read) < chunk))
			chunk = (maxPos - minPos + 1) - read;
		
		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)read / (double)(maxPos - minPos + 1) << "  \t% complete.\r";
	}
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;

	// Check if we have read past the EOF so the eofbit is set
	if (infF->eof())
		infF->clear();
	if (infR->eof())
		infR->clear();

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
	
	chunk = CHUNK_MULTI*MAXNR;
	// if necessary, reduce to file size.
	chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);

	read = 0;
	
	// Populate count vectors
	vcerr(3) << "\t\t" << "Counting value ocurrences" << endl;
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";
	while ((maxPos - minPos + 1) - read >= chunk && chunk > 0)
	{
		// Read data from F file
		infF->seekg(sizeof(storageType) * (read + (minPos-minF)) + sizeof(int), ios::beg);
		infF->read((char*)tmpF, sizeof(storageType) * chunk);
		// Read data from R file
		infR->seekg(sizeof(storageType) * (read + (minPos-minR)) + sizeof(int), ios::beg);
		infR->read((char*)tmpR, sizeof(storageType) * chunk);

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
		if (((maxPos - minPos + 1 - read) < chunk))
			chunk = (maxPos - minPos + 1) - read;
		
		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)read / (double)(maxPos - minPos + 1) << "  \t% complete.\r";
	}
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;

	// Check if we have read past the EOF so the eofbit is set
	if (infF->eof())
		infF->clear();
	if (infR->eof())
		infR->clear();

	if(truncLimit == 1)
	  {
	    *truncValF = maxValF;
	    *truncValR = maxValR;
	  }else{
	  // Identify truncate values
	  vcerr(3) << "\t\t" << "Identifying truncate values" << endl;
	  
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
	}
	vcerr(3) << "\t\t" << "truncValF: " << (*truncValF) << " truncValR: " << (*truncValR) << endl;
	
	// Write counts to outf
	for (int i = 1; i <= max(maxValF, maxValR); i++)
		(*outf) << i << "\t" << (i > maxValF ? 0 : countsF[i]) << "\t" << (i > maxValR ? 0 : countsR[i]) << endl;

	return(1);
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

	//minPos = 121150000;
	//maxPos = 121200000;
	
	chunk = CHUNK_MULTI*MAXNR;
	// if necessary, reduce to file size.
	chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);

	read = 0;

	//ofstream ofmin;
	//ofstream ofone;
	//ofstream off;
	//ofstream ofr;
	//ofmin.open("oMin.txt",ios::trunc);
	//ofone.open("oOne.txt",ios::trunc);
	//off.open("oF.txt",ios::trunc);
	//ofr.open("oR.txt",ios::trunc);

	double corr[maxLag-minLag+1];
	double tmpcorr[maxLag-minLag+1];
	long testmin[maxLag-minLag+1];
	long testone[maxLag-minLag+1];
	long testf[maxLag-minLag+1];
	long testr[maxLag-minLag+1];

	// calculate auto-correlations while we're at it
	double autoCorrF[maxLag-minLag+1];
	double autoTmpCorrF[maxLag-minLag+1];
	long autoPairF[maxLag-minLag+1];

	double autoCorrR[maxLag-minLag+1];
	double autoTmpCorrR[maxLag-minLag+1];
	long autoPairR[maxLag-minLag+1];
	
	
	for (int lag = minLag; lag <= maxLag; lag++)
	{
		corr[lag-minLag]  = 0.0;
		distDens[lag-minLag]  = 0.0;
		testmin[lag-minLag]  = 0;
		testone[lag-minLag]  = 0;
		testf[lag-minLag]  = 0;
		testr[lag-minLag]  = 0;

		autoCorrF[lag - minLag] = 0;
		autoTmpCorrF[lag - minLag] = 0;
		autoPairF[lag - minLag] = 0;
		autoCorrR[lag - minLag] = 0;
		autoTmpCorrR[lag - minLag] = 0;
		autoPairR[lag - minLag] = 0;
	}
	
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
		  {
		    tmpcorr[lag-minLag]  = 0.0;
		    autoTmpCorrF[lag-minLag]  = 0.0;
		    autoTmpCorrR[lag-minLag]  = 0.0;
		    
		  }
		//Loop all F positions
		bool fabove = false;
		bool rabove = false;
		for (int i = 0; i <= chunk - maxLag; i++)
		{
			fabove = (int)tmpF[i] > 0;
			//Traverse from minLag to maxLag
			for (int lag = minLag; lag <= maxLag; lag++)
			{
				rabove = (int)tmpR[i+lag-1] > 0;
				tmpcorr[lag-minLag] += ((double)tmpF[i] - meanF) * ((double)tmpR[i+lag-1] - meanR);
				autoTmpCorrF[lag-minLag] += ((double)tmpF[i] - meanF) * ((double)tmpF[i+lag-1] - meanF);
				autoTmpCorrR[lag-minLag] += ((double)tmpR[i] - meanR) * ((double)tmpR[i+lag-1] - meanR);
				if(tmpF[i+lag-1]>0 && tmpF[i] > 0)
				  autoPairF[lag-minLag] += 1;
				if(tmpR[i+lag -1] > 0 && tmpR[i]>0)
				  autoPairR[lag-minLag] += 1;
				if (fabove && rabove)
				{
				  testmin[lag-minLag] += min(tmpF[i], tmpR[i+lag-1]);
				  testone[lag-minLag] += 1;
				  testf[lag-minLag] += tmpF[i];
				  testr[lag-minLag] += tmpR[i+lag-1];
				}
			}
			
			sdF += pow((double)tmpF[i] - meanF,2);
			sdR += pow((double)tmpR[i] - meanR,2);
		}

		for (int lag = minLag; lag <= maxLag; lag++)
		  {
		    corr[lag-minLag] = (corr[lag-minLag] * ((double)read / (double)(read+chunk-maxLag+1))) +
		      (tmpcorr[lag-minLag] / (double)(read+chunk-maxLag+1));
		    autoCorrF[lag-minLag] = (autoCorrF[lag-minLag] * ((double)read / (double)(read+chunk-maxLag+1))) +
		      (autoTmpCorrF[lag-minLag] / (double)(read+chunk-maxLag+1));
		    autoCorrR[lag-minLag] = (autoCorrR[lag-minLag] * ((double)read / (double)(read+chunk-maxLag+1))) +
		      (autoTmpCorrR[lag-minLag] / (double)(read+chunk-maxLag+1));
		  }


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
	*(outf) << "lag\t" << "pearson\t" << "pearson_F\t" << "pearson_R\t";
	*(outf) << "pair_one\t" << "pair_one_F\t" << "pair_one_R\t";
	*(outf) << "pair_min_F_R" << "\t" << "pair_F" << "\t" << "pair_R" <<endl;
	for (int lag = minLag; lag <= maxLag; lag++)
	  {
		corr[lag-minLag] /= sdF * sdR;
		autoCorrF[lag-minLag] /= sdF * sdF;
		autoCorrR[lag-minLag] /= sdR * sdR;
		sum += testone[lag-minLag];
		//*(outf) << lag << "\t" << corr[lag-minLag] << endl;
		*(outf) << lag << "\t" << corr[lag-minLag] <<"\t" << autoCorrF[lag-minLag] << "\t";
		*(outf) << autoCorrR[lag-minLag] <<"\t"<<testone[lag-minLag] << "\t";
		*(outf) << autoPairF[lag-minLag] <<"\t"<<autoPairR[lag-minLag]<< "\t";
		*(outf) << testmin[lag-minLag] << "\t" << testf[lag-minLag] << "\t" << testr[lag-minLag] << endl;
	  }
	
	//ofmin.close();	
	//ofone.close();
	//off.close();
	//ofr.close();

	double maxVal = 0.0;
	for (int lag = minLag; lag <= maxLag; lag++)
	{
		distDens[lag-minLag] = testone[lag-minLag] / sum;
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
	
	delete[] tmpF;
	delete[] tmpR;
	
	return(1);
}

int windowSums(ifstream* infF, ifstream* infR, int* centers,
			   int sampleSize, int maxPos, int minPos,
			   int minF, int minR, int minD, int maxD,
			   int* maxValF, int* maxValR,
			   int thetaF1, int thetaF2, int thetaR1, int thetaR2,
			   int* F, int* R,
			   double* meanF, double* meanR,
			   double* varF, double* varR,
			   double* covFR,
			   int truncValF, int truncValR)
{
	int chunk = CHUNK_MULTI*MAXNR;
	// if necessary, reduce to file size.
	chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);

	storageType* tmpF = new storageType[chunk];
	storageType* tmpR = new storageType[chunk];
	
	*meanF = 0.0;
	*meanR = 0.0;
	
	int sampleCount = 0;
	int read = 0;

	int w = thetaF1 - thetaF2 + 1;

	*maxValF = 0;
	*maxValR = 0;
	
	// Assign Fi and Ri
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";
	while (((maxPos - minPos + 1) - read >= chunk) && (chunk - maxD >= 0) && (sampleSize > sampleCount))
	{
		if (centers[sampleCount] + thetaR2 <= read + chunk)
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

			while ((centers[sampleCount] < read + chunk) && (sampleSize > sampleCount) && (centers[sampleCount] + thetaR2 <= read + chunk))
			{
				F[sampleCount] = 0;
				R[sampleCount] = 0;

				for (int j = 0; j < w; j++)
				{
					F[sampleCount] += tmpF[centers[sampleCount] - thetaF1 + j - read];
					R[sampleCount] += tmpR[centers[sampleCount] + thetaR1 + j - read];
				}

				if (F[sampleCount] > (*maxValF))
					*maxValF = F[sampleCount];

				if (R[sampleCount] > (*maxValR))
					*maxValR = R[sampleCount];
				
				// Update the mean
				*meanF = ((*meanF) * ((double)sampleCount / (double)(sampleCount+1))) +
					(F[sampleCount] / (double)(sampleCount+1));
				*meanR = ((*meanR) * ((double)sampleCount / (double)(sampleCount+1))) +
					(R[sampleCount] / (double)(sampleCount+1));

				sampleCount++;
			}
		}

		read += chunk-thetaF1-1;

		// Check if we can read in a full chunk of data or simply the rest of the file. 
		if (((maxPos - minPos + 1 - read) < chunk) && ((maxPos - minPos + 1 - read) >= maxD))
			chunk = (maxPos - minPos + 1) - read;

		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)sampleCount / (double)sampleSize << "  \t% complete.\r";
	}
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;
	delete[] tmpF;
	delete[] tmpR;

	// Calculate variance and covariance
	*varF = 0.0;
	*varR = 0.0;
	*covFR = 0.0;
	double mult;
	for (int i = 0; i < sampleSize; i++)
	{
		mult = ((double)i / (double)(i+1));
		*varF = ((*varF) * mult) + (pow((double)F[i] - (*meanF), 2) / (double)(i+1));
		*varR = ((*varR) * mult) + (pow((double)R[i] - (*meanF), 2) / (double)(i+1));
		*covFR = ((*covFR) * mult) + (((double)F[i] - (*meanF)) * ((double)R[i] - (*meanR)) / (double)(i+1));
	}

	// Check if we have read past the EOF so the eofbit is set
	if (infF->eof())
		infF->clear();
	if (infR->eof())
		infR->clear();
	
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
		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)iter / (double)(burnins + iterations) << "  \t% complete.\r";
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
		}

		criteriaMet = false;
		// TODO: Make sure that this loop don't get stuck somehow
		while (!criteriaMet)
		{
			lambdaRWork[0] = gsl_ran_gamma(r, alphaR + (double)RSumE, 1.0 / (betaR + (double)nRE));
			lambdaRWork[1] = gsl_ran_gamma(r, alphaR + (double)RSumNotE, 1.0 / (betaR + (double)nRNotE));

			criteriaMet = (lambdaRWork[0] > lambdaRWork[1]);
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

int readParameters(ifstream* ipff, double* pF, double* pR, double* lambdaF, double* lambdaR,int* lambdaD,int* minD,int* maxD,
				   int* estLambdaD,int* estMinD,int* estMaxD,vector<string>* corrv)
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


int writeParameters(ofstream* opff, double* pF, double* pR, double* lambdaF, double* lambdaR, int lambdaD,int minD,int maxD,
		    int estLambdaD,int estMinD,int estMaxD,double* corr)
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

int centerPositionOdds(ifstream* infF, ifstream* infR,
		       ofstream* outf,
		       double* pF, double* pR, double* lambdaF, double* lambdaR,
		       int minD, int maxD,
		       int* maxValF, int* maxValR,
		       int startCrd, int endCrd,
		       int minPos, int maxPos,
		       int minR, int minF,
		       int maxF, int maxR,
		       int truncValF, int truncValR)
{	
	//if(startCrd >=0)
	//  {
	//    minPos = max(minPos,startCrd);
	//    vcerr(2) << "\tLimiting analysis (possibly restricted by actual genome coordinates). Staring from: "<<minPos<<endl;
	//  }
	
	//	if(endCrd >=0)
	//  {
	//    maxPos = min(maxPos,endCrd);
	//    vcerr(2) << "\tLimiting analysis (possibly restricted by actual genome coordinates). Ending at: "<<maxPos<<endl;
	//  }



	/*
	 * 2. Assign the data range variables
	 */

	// Window size of reads
	int w = ceil((double)(maxD - minD) / 2.0);
	// Offsets from center position
	int thetaF1, thetaF2, thetaR1, thetaR2;
	thetaF1 = thetaR2 = ceil((double)(maxD - 1) / 2.0);
	thetaF2 = thetaR1 = thetaF1 - w + 1;
	maxD = thetaF1 + thetaR2 + 1;

	/*
	 * 3. Slide the files and calculate odds of center positions
	 */
	
	double meanF[2] = {0.0, 0.0};
	double meanR[2] = {0.0, 0.0};

	long cntF[2] = {0,0};
	long cntR[2] = {0,0};
		
	// current window sums.
	long wsumF = 0, wsumR = 0; 

	int chunk = CHUNK_MULTI*MAXNR;
	// if necessary, reduce to file size.
	chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);

	storageType* tmpF = new storageType[chunk];
	storageType* tmpR = new storageType[chunk];

	int read = 0;
	double oddsF;
	double oddsR;
	int idF;
	int idR;

	*maxValF = 0;
	*maxValR = 0;

	vcerr(2) << "\t* Looking up max value" << endl;
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";
	while (((maxPos - minPos + 1) - read >= chunk) && (chunk - maxD >= 0))
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

		// loop over the read data and calculate sums over the window size. 
		wsumF = 0;
		wsumR = 0;
		for (int i = thetaF1-w+1; i < chunk-thetaR2; i++)
		{
			wsumF += (long)tmpF[i-thetaF2];      // window sum
			wsumR += (long)tmpR[i+w+thetaR1-1];  // window sum
			
			if (i >= thetaF1)
			{
				*maxValF = max((int)wsumF, (*maxValF));
				*maxValR = max((int)wsumR, (*maxValR));
				
				// decrease the window sum with the value just before the current window.
				wsumF -= (long)tmpF[i-thetaF1];
				// decrease the window sum with the value just before the current window.
				wsumR -= (long)tmpR[i+thetaR1];
			}
		}
		read += chunk-thetaF1-thetaR2;
		
		// Check if we can read in a full chunk of data or simply the rest of the file. 
		if (((maxPos - minPos + 1 - read) < chunk) && ((maxPos - minPos + 1 - read) >= maxD))
			chunk = (maxPos - minPos + 1) - read;

		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)read / (double)(maxPos - minPos + 1) << "  \t% complete.\r";
	}
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;

	// Check if we have read past the EOF so the eofbit is set
	if (infF->eof())
		infF->clear();
	if (infR->eof())
		infR->clear();
	
	double oddsFs[*maxValF];
	double oddsRs[*maxValR];

	vcerr(3) << "\t\tF max = " << *maxValF << " R max = " << *maxValR << endl;

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
	read = 0;
	chunk = MAXNR;
	chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);

	// Write min coordinate to odds files:
	outf->write((char*)&minPos, sizeof(int));

	// Fill with zeros up to thetaF1 - 1
	storageType zero = 0;
	for (int i = 0; i < thetaF1; i++)
		outf->write((char*)&zero, sizeof(storageType));
	
	vcerr(2) << "\t* Calculating odds of center positions" << endl;
	storageType discreteOdds[chunk-thetaF1-thetaR2];
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";
	int pos = thetaF1-1;
	bool newchunk = false;
	int written = thetaF1;
	double tmp;
	while (((maxPos - minPos + 1) - read >= chunk) && (chunk - maxD >= 0))
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

				//discreteOdds[i-thetaF1] = (oddsF > 0.0 && oddsR > 0.0 ? 1 : 0);
				if (oddsF > 0.0 && oddsR > 0.0)
				{
					tmp = (oddsF+oddsR) / maxOdds;
					discreteOdds[i-thetaF1] = np_round(1.0 + (99.0 * tmp));
				}
				else
					discreteOdds[i-thetaF1] = 0;
				pos++;
				
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
			}
		}

		// Write the chunk of calls to file
		if (chunk == MAXNR)
		{
			outf->write((char*)discreteOdds, sizeof(discreteOdds));
		}
		else
		{
			for (int i = 0; i < chunk-thetaF1-thetaR2; i++)
				outf->write((char*)&discreteOdds[i], sizeof(storageType));
		}
		written += chunk-thetaF1-thetaR2;

		if (!newchunk)
			read += chunk-thetaF1-thetaR2;
		else
			read += chunk;

		// Check if we can read in a full chunk of data or simply the rest of the file. 
		if (((maxPos - minPos + 1 - read) < chunk) && ((maxPos - minPos + 1 - read) >= maxD))
		{
			chunk = (maxPos - minPos + 1) - read;
			newchunk = true;
		}

		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)read / (double)(maxPos - minPos + 1) << "  \t% complete.\r";
	}
	
	// Fill odds file with ending zeros
	for (int i = pos+1; i < (maxPos-minPos+1); i++)
		outf->write((char*)&zero, sizeof(storageType));
	
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;

	// Check if we have read past the EOF so the eofbit is set
	if (infF->eof())
		infF->clear();
	if (infR->eof())
		infR->clear();

	// Calculate variance and covariance
	read = 0;
	chunk = CHUNK_MULTI*MAXNR;
	chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);

	double varF[2] = {0.0, 0.0};
	double varR[2] = {0.0, 0.0};
	double covFR[2] = {0.0, 0.0};
	long cntFR[2] = {0,0};
	int idFR;
	cntF[0] = 0;
	cntF[1] = 0;
	cntR[0] = 0;
	cntR[1] = 0;
	
	vcerr(2) << "\t* Calculating statistics" << endl;
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";
	while (((maxPos - minPos + 1) - read >= chunk) && (chunk - maxD >= 0))
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
		if (((maxPos - minPos + 1 - read) < chunk) && ((maxPos - minPos + 1 - read) >= maxD))
			chunk = (maxPos - minPos + 1) - read;

		vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)read / (double)(maxPos - minPos + 1) << "  \t% complete.\r";
	}
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;

	vcerr(3) << "\t\tStatistics:" << endl;
	vcerr(3) << "\t\tF: mean = (" <<  meanF[0] << ", " << meanF[1] << ")" << endl;
	vcerr(3) << "\t\tR: mean = (" <<  meanR[0] << ", " << meanR[1] << ")" << endl;
	vcerr(3) << "\t\tF: variance = (" <<  varF[0] << ", " << varF[1] << ")" << endl;
	vcerr(3) << "\t\tR: variance = (" <<  varR[0] << ", " << varR[1] << ")" << endl;
	vcerr(3) << "\t\tCovariance = (" <<  covFR[0] << ", " << covFR[1] << ")" << endl;
	
	delete[] tmpF;
	delete[] tmpR;
	
	// Check if we have read past the EOF so the eofbit is set
	if (infF->eof())
		infF->clear();
	if (infR->eof())
		infR->clear();

	return(1);
}


int predictBoundaries(ifstream* infF, ifstream* infR,
					  ifstream* infO, ofstream* outf,
					  double* pF, double* pR, double* lambdaF, double* lambdaR,
					  int lambdaD, int minD, int maxD,
					  int minLag, int maxLag,
					  int maxValF, int maxValR,
					  int startCrd, int endCrd,
					  int minPos, int maxPos,
					  int minR, int minF,
					  int maxF, int maxR,
					  int truncValF, int truncValR,
					  double* distDens)
{
	// Window size of reads
	int w = ceil((double)(maxD - minD) / 2.0);
	// Offsets from center position
	int thetaF1, thetaF2, thetaR1, thetaR2;
	thetaF1 = thetaR2 = ceil((double)(maxD - 1) / 2.0);
	thetaF2 = thetaR1 = thetaF1 - w + 1;
	maxD = thetaF1 + thetaR2 + 1;
		
	int file_minPos = 0;
	file_minPos = minPos;

	//cerr << minF << " " << maxF << endl;
	//cerr << minR << " " << maxR << endl;
	//cerr << minPos << " " << maxPos << endl;
	
	if(startCrd >0)
	  {
	    minPos = max(minPos,startCrd);
	    vcerr(2) << "\tLimiting analysis (possibly restricted by actual genome coordinates). Staring from: "<<minPos<<endl;
	  }

	if(endCrd != -1)
	  {
	    maxPos = min(maxPos,endCrd);
	    vcerr(2) << "\tLimiting analysis (possibly restricted by actual genome coordinates). Ending at: "<<maxPos<<endl;
	  }

	double oddsFs[maxValF];
	double oddsRs[maxValR];

	for (int i = 0; i < maxValF; i++)
		oddsFs[i] = log(max(gsl_ran_poisson_pdf(i, lambdaF[0]) * pF[0], DBL_MIN)) -
			log(max(gsl_ran_poisson_pdf(i, lambdaF[1]) * pF[1], DBL_MIN));
	for (int i = 0; i < maxValR; i++)
		oddsRs[i] = log(max(gsl_ran_poisson_pdf(i, lambdaR[0]) * pR[0], DBL_MIN)) -
			log(max(gsl_ran_poisson_pdf(i, lambdaR[1]) * pR[1], DBL_MIN));

	NPPredictor* pred =
		new NPPredictor(outf, thetaF1,
						minD, maxD, minLag, maxLag,
						infF, infR,
						file_minPos, maxPos, thetaF2, distDens,
						pF, pR, lambdaF, lambdaR);

	int chunk = CHUNK_MULTI*MAXNR;
	// if necessary, reduce to file size.
	chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);

	storageType* tmpF = new storageType[chunk];
	storageType* tmpR = new storageType[chunk];
	storageType* tmpO = new storageType[chunk];
	
	int read = 0;
	// current window sums.
	long wsumF = 0, wsumR = 0; 
	double oddsF;
	double oddsR;

	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.\r";
	while (((maxPos - minPos + 1) - read >= chunk) && (chunk - maxD >= 0))
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
	    
	    // Read data from O file
	    infO->seekg(sizeof(storageType) * (read + (minPos - file_minPos)) + sizeof(int), ios::beg);
	    infO->read((char*)tmpO, sizeof(storageType) * chunk);

	    // reset the wsums since the 'read += chunk-thetaF1-thetaR2;' below accounts for the "lost" values. 
	    wsumF = 0;
	    wsumR = 0;
	    for (int i = thetaF1-w+1; i < chunk-thetaR2; i++)
	      {
		wsumF += (long)tmpF[i-thetaF2];      // window sum
		wsumR += (long)tmpR[i+w+thetaR1-1];  // window sum
		
		// if we've read the window size the counts are valid.
		if (i >= thetaF1)
		  {
		    //if (tmpO[i] == 1)
		    if (tmpO[i] > 0)
		      {
			oddsF = oddsFs[(int)wsumF];
			oddsR = oddsRs[(int)wsumR];
			
			pred->addPrediction((double)(i-thetaF1+read+minPos),
					    (double)(i+thetaR2+read+minPos),
					    oddsF+oddsR, i+read+minPos,0);
		      }
		    
		    // decrease the window sum with the value just before the current window.
		    wsumF -= (long)tmpF[i-thetaF1];
		    // decrease the window sum with the value just before the current window.
		    wsumR -= (long)tmpR[i+thetaR1];
		    
		    if ((i + read) % 1000 == 1)
		      vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)(i + read) / (double)(maxPos - minPos + 1) << "  \t% complete.\r";
		  }
	      }
	    read += chunk-thetaF1-thetaR2;
	    
	    // Old variant:
	    // 		for (int i = thetaF1; i < chunk-thetaR2; i++)
	    // 	    {
	    // 			if (tmpO[i] == 1)
	    // 			{
	    
	    // 				pred->addPrediction((double)(i-thetaF1+read+minPos),
	    // 									(double)(i+thetaR2+read+minPos),
	    // 									1.0, 1.0, 1.0, i+read+minPos,0);
	    // 			}
	    
	    // 			if ((i + read) % 1000 == 1)
	    // 				vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)(i + read) / (double)(maxPos - minPos + 1) << "  \t% complete.\r";
	    // 	    }
	    // 		read += chunk-thetaF1-thetaR2;
	    
	    // Check if we can read in a full chunk of data or simply the rest of the file. 
	    if (((maxPos - minPos + 1 - read) < chunk) && ((maxPos - minPos + 1 - read) >= maxD))
	      chunk = (maxPos - minPos + 1) - read;
	  
	    vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)read / (double)(maxPos - minPos + 1) << "  \t% complete.\r";
	  }
	vcerr(3) << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;
	
	delete pred;
	delete[] tmpF;
	delete[] tmpR;
	delete[] tmpO;
	
	return(1);
}

