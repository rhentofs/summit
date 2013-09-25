int centerPositionOdds(ifstream* infF, ifstream* infR,
					   ofstream* outo, ofstream* outf,
					   ofstream* outs,
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

	/*
	 * 2. Slide the files and calculate odds of center positions
	 */
	
	double meanF[2] = {0.0, 0.0};
	double meanR[2] = {0.0, 0.0};

	long cntF[2] = {0,0};
	long cntR[2] = {0,0};
		
	// current window sums.
	long wsumF = 0, wsumR = 0; 

	storageType* tmpF;
	storageType* tmpR;

	int read = 0;
	double oddsF;
	double oddsR;
	int idF;
	int idR;

	double maxVarFR = 0.0;
	double wMeanF = 0.0;
	double wMeanR = 0.0;
	double wVarF = 0.0;
	double wVarR = 0.0;
	int wCntF = 0;
	int wCntR = 0;
	
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
			
			wMeanF = 0.0;
			wMeanR = 0.0;
			
			for (int i = thetaF1-w+1; i < chunk-thetaR2; i++)
			{
				wMeanF = ((wMeanF * (double)(wsumF)) + ((double)(i-thetaF2) * (double)tmpF[i-thetaF2]));
				wMeanR = ((wMeanR * (double)(wsumR)) + ((double)(i+w+thetaR1-1) * (double)tmpR[i+w+thetaR1-1]));

				wsumF += (long)tmpF[i-thetaF2];      // window sum
				wsumR += (long)tmpR[i+w+thetaR1-1];  // window sum

				//cerr << wsumF << "\t" << wsumR << endl;
				//cerr << wMeanF << "\t" << wMeanR << "\t";

				if (wsumF != 0)
					wMeanF /= (double)wsumF;
				else
					wMeanF = 0.0;
				if (wsumR != 0)
					wMeanR /= (double)wsumR;
				else
					wMeanR = 0.0;
								
				if (i >= thetaF1)
				{
					*maxValF = max((int)wsumF, (*maxValF));
					*maxValR = max((int)wsumR, (*maxValR));
					
					if ((wsumF > 0) || (wsumR > 0))
					{
						//cerr << (int)wMeanF << "\t" << (int)wMeanR << "\t";
						wVarF = 0.0;
						wVarR = 0.0;
						wCntF = 0;
						wCntR = 0;
					
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

						//cerr << wVarF << "\t" << wVarR << "\t";

						if (wCntF != 0)
							wVarF /= (double)wCntF;
						else
							wVarF = 0.0;
						if (wCntR != 0)
							wVarR /= (double)wCntR;
						else
							wVarR = 0.0;

						//cerr << (int)wVarF << "\t" << (int)wVarR << endl;
						
						if ((wVarF + wVarR) > maxVarFR)
							maxVarFR = wVarF + wVarR;
					}
					
					wMeanF = (wMeanF * (double)wsumF) - ((double)tmpF[i-thetaF1] * (double)(i-thetaF1));
					wMeanR = (wMeanR * (double)wsumR) - ((double)tmpR[i+thetaR1] * (double)(i+thetaR1));

					//cerr << wMeanF << "\t" << wMeanR << "\t";
					
					// decrease the window sum with the value just before the current window.
					wsumF -= (long)tmpF[i-thetaF1];
					// decrease the window sum with the value just before the current window.
					wsumR -= (long)tmpR[i+thetaR1];

					if (wsumF != 0)
						wMeanF /= (double)wsumF;
					else
						wMeanF = 0.0;
					
					if (wsumR != 0)
						wMeanR /= (double)wsumR;
					else
						wMeanR = 0.0;

					//cerr << wMeanF << "\t" << wMeanR << endl;
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

	(*outs) << "Max Fuzzy score (var F + var R) = " << maxVarFR << " normalized to the [0,100] interval." << endl;
	(*outs) << "Max F window sum = " << *maxValF << endl;
	(*outs) << "Max R window sum = " << *maxValR << endl;

	vcerr(3) << "\t\tMax Fuzzy score (var F + var R) = " << maxVarFR << " normalized to the [0,100] interval." << endl;
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
	(*outs) << " (" << maxOddsF << " + " << maxOddsR << ") normalized to the [0,100] interval." << endl;

	vcerr(3) << "\t\tMax center position odds (odds start + odds end) = " << maxOdds;
	vcerr(3) << " (" << maxOddsF << " + " << maxOddsR << ") normalized to the [0,100] interval." << endl;
	
	int pos = thetaF1-1;
	bool newchunk = false;
	int written = thetaF1;
	double tmp;
	
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
		
		storageType discreteOdds[chunk-thetaF1-thetaR2];
		storageType discreteFuzzy[chunk-thetaF1-thetaR2];
		
		// Write min coordinate to odds files:
		outo[iter].write((char*)&minPos[iter], sizeof(int));
		// Write min coordinate to fuzzy files:
		outf[iter].write((char*)&minPos[iter], sizeof(int));

		// Fill with zeros up to thetaF1 - 1
		storageType zero = 0;
		for (int i = 0; i < thetaF1; i++)
		{
			outo[iter].write((char*)&zero, sizeof(storageType));
			outf[iter].write((char*)&zero, sizeof(storageType));
		}
		
		pos = thetaF1-1;
		newchunk = false;
		written = thetaF1;
		
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

			wMeanF = 0.0;
			wMeanR = 0.0;

			for (int i = thetaF1-w+1; i < chunk-thetaR2; i++)
			{
				wMeanF = ((wMeanF * (double)(wsumF)) + ((double)(i-thetaF2) * (double)tmpF[i-thetaF2]));
				wMeanR = ((wMeanR * (double)(wsumR)) + ((double)(i+w+thetaR1-1) * (double)tmpR[i+w+thetaR1-1]));
				
				wsumF += (long)tmpF[i-thetaF2];      // window sum
				wsumR += (long)tmpR[i+w+thetaR1-1];  // window sum

				if (wsumF != 0)
					wMeanF /= (double)wsumF;
				else
					wMeanF = 0.0;
				if (wsumR != 0)
					wMeanR /= (double)wsumR;
				else
					wMeanR = 0.0;
				
				// if we've read the window size the counts are valid.
				if (i >= thetaF1)
				{
					oddsF = oddsFs[(int)wsumF];
					oddsR = oddsRs[(int)wsumR];

					if (oddsF > 0.0 && oddsR > 0.0)
					{
						tmp = (oddsF+oddsR) / maxOdds;
						discreteOdds[i-thetaF1] = np_round(1.0 + (99.0 * tmp));
					}
					else
					{
						discreteOdds[i-thetaF1] = 0;
					}
					pos++;

					if ((wsumF > 0) || (wsumR > 0))
					{
						wVarF = 0.0;
						wVarR = 0.0;
						wCntF = 0;
						wCntR = 0;
					
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

						//cerr << wVarF << "\t" << wVarR << "\t";

						if (wCntF != 0)
							wVarF /= (double)wCntF;
						else
							wVarF = 0.0;
						if (wCntR != 0)
							wVarR /= (double)wCntR;
						else
							wVarR = 0.0;
						
						tmp = (wVarF + wVarR) / maxVarFR;
						if (tmp > 1.0)
							cerr << tmp << "\t" << wVarF << "\t" << wVarR << "\t" << maxVarFR << endl;
						discreteFuzzy[i-thetaF1] = np_round(1.0 + (99.0 * tmp));
					}
					else
					{
						discreteFuzzy[i-thetaF1] = 0;
					}
					
					wMeanF = (wMeanF * (double)wsumF) - ((double)tmpF[i-thetaF1] * (double)(i-thetaF1));
					wMeanR = (wMeanR * (double)wsumR) - ((double)tmpR[i+thetaR1] * (double)(i+thetaR1));
					
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

					if (wsumF != 0)
						wMeanF /= (double)wsumF;
					else
						wMeanF = 0.0;
					if (wsumR != 0)
						wMeanR /= (double)wsumR;
					else
						wMeanR = 0.0;
				}
			}
			
			// Write the chunk of calls to file
			if (chunk == MAXNR)
			{
				outo[iter].write((char*)discreteOdds, sizeof(discreteOdds));
				outf[iter].write((char*)discreteFuzzy, sizeof(discreteFuzzy));
			}
			else
			{
				for (int i = 0; i < chunk-thetaF1-thetaR2; i++)
				{
					outo[iter].write((char*)&discreteOdds[i], sizeof(storageType));
					outf[iter].write((char*)&discreteFuzzy[i], sizeof(storageType));
				}
			}
			written += chunk-thetaF1-thetaR2;

			if (!newchunk)
				read += chunk-thetaF1-thetaR2;
			else
				read += chunk;
			
			// Check if we can read in a full chunk of data or simply the rest of the file. 
			if (((maxPos[iter] - minPos[iter] + 1 - read) < chunk) && ((maxPos[iter] - minPos[iter] + 1 - read) >= maxD))
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
