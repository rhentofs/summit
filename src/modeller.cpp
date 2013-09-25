#include "zeb_utils.h"

int vlevel = 3;
int main(int argc, char* argv[]) 
{
	// parameters that you can set. 
	string infF    = "";
	string infR    = "";
	string infPath = "";
	string outI    = "oIter.txt";
	string outP    = "oPar.txt";
	string outCorr = "oCorr.txt";
	string outCnt  = "oCnt.txt";
	string outPath = "";
	string outStub = "";
	//bool useemp    = false;
	double truncLimit = 0.99;
	int verbose    = 3;
	int lambdaD    = 147;
	int minD       = 130;
	int maxD       = 180;
	int sampleSize = -1;
	int burnins    = 1000;
	int iterations = 10000;
	int seed       = -1;
	// Unknowns
	double lambdaF[2];
	double lambdaR[2];
	double pF[2];
	double pR[2];

	/*
	 *  Read and check input parameters
	 */

	string errorLine =  "usage " + 
		string(argv[0]) + 
		" [parameters] \n" +
		"\t-iF      <forward reads binary file> \n" +
		"\t-iR      <reverse reads binary file> \n" +
		"\t-iPath   <path to forward and reverse reads binary files \n" +
		"            (overrides iR and iF if set)> \n" +
		"\t-oPath   <outdirectory, where all output files are put. \n" + 
		"\t          NOT created - needs to exists. Defaults to where \n"+
		"\t          the program is executed from.>\n" +
		"\t-oStub   <outfile names stub, \n" +
		"\t          defaults -oIter to <-oStub_>oIter.txt, \n" + 
		"\t          defaults -oPar to <-oStub_>oPar.txt> \n" + 
		"\t          defaults -oCorr to <-oStub_>oCorr.txt> \n" +
		"\t          defaults -oCnt to <-oStub_>oCnt.txt> \n" + 
		"\t-oIter   <outfile, where to write MCMC parameter\n" +
		"\t          estimates for each iteration,default generated from -oStub>\n" +
		"\t-oPar    <parameter-file, where to write MCMC\n" +
		"\t          final parameter estimates, default generated from -oStub> \n" +
         /*
		   "\t-oCorr   <out-file, lists correlation coefficients between forward\n" +
		   "\t          and reverse reads in [mind,maxd] \n" +
		   "\t          default generated from -oStub> \n" +
		 */
		"\t-oCnt    <out-file, lists forward and reverse read count frequencies\n" +
		"\t          default generated from -oStub> \n" +
		"\t-trunc   <truncation limit in (0,1], default 0.99.> \n" +
		"\t-size    <sample size, default all observations.> \n" +
		"\t-burnin  <number of MCMC burnin iterations,\n" +
		"\t          default 1000.> \n" +
		"\t-iter    <number of MCMC iterations\n" +
		"\t          (non-burnin iterations), default 10000.> \n" +
		"\t-seed    <set seed to random number generator.> \n" +
		"\t-ld      <distance-lambda, default 147 bp.> \n" +
		"\t-mind    <min distance, default 130 bp.> \n" +
		"\t-maxd    <max distance, default 180 bp.> \n" +
		/*
		  "\t-useemp  <toggle use of data-driven distance distribution, or poisson\n" +
		  "\t          around distance-lambda. default off.>\n" +
		*/
		"\t-v       <verbose level 0-3. default 2>";
	
	bool fail = false;
	string failmessage = "";
	
	for (int i = 1; i < argc; i++)
	{
	    if(strcmp(argv[i],"-iF") == 0)
			infF.assign(argv[++i]);
	    else if(strcmp(argv[i],"-iR") == 0)
			infR.assign(argv[++i]);
		else if(strcmp(argv[i],"-iPath") == 0)
			infPath.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oPath") == 0)
			outPath.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oStub") == 0)
			outStub.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oIter") == 0)
			outI.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oPar") == 0)
			outP.assign(argv[++i]);
	    /*
		  else if(strcmp(argv[i],"-oCorr") == 0)
			outCorr.assign(argv[++i]);
		*/
	    else if(strcmp(argv[i],"-oCnt") == 0)
			outCnt.assign(argv[++i]);
	    /*
		  else if (strcmp(argv[i],"-useemp") == 0)
			useemp = true;
		*/
	    else if (strcmp(argv[i],"-v") == 0)
			verbose = atoi(argv[++i]);
	    else if (strcmp(argv[i],"-seed") == 0)
			seed = atoi(argv[++i]);
	    else if (strcmp(argv[i],"-trunc") == 0)
			truncLimit = atof(argv[++i]);
	    else if (strcmp(argv[i],"-size") == 0)
			sampleSize = atoi(argv[++i]);
	    else if (strcmp(argv[i],"-burnin") == 0)
			burnins = atoi(argv[++i]);
	    else if (strcmp(argv[i],"-iter") == 0)
			iterations = atoi(argv[++i]);
	    else if (strcmp(argv[i],"-ld") == 0)
			lambdaD = atoi(argv[++i]);
	    else if (strcmp(argv[i],"-mind") == 0)
			minD = atoi(argv[++i]);
	    else if (strcmp(argv[i],"-maxd") == 0)
			maxD = atoi(argv[++i]);
	    else
		{
			failmessage.assign("Unknown argument: ");
			failmessage.append(argv[i]);
			failmessage.append("\n");
			fail = true;
		}
	}
	
	if (truncLimit <= 0 || truncLimit > 1)
	{
	    failmessage.append("-trunc value does not make sense.\n");
	    fail = true;
	}

	bool infPathSpec = false;
	if (strcmp(infPath.c_str(), "") != 0)
	{
		infPathSpec = true;
	    DIR *d = opendir(infPath.c_str());
	    if(d)
		{
			closedir(d);
		}
	    else
		{
			failmessage.append("-iPath does not exist.\n");
			fail = true;
		}
	}
	
	if (strcmp(infF.c_str(), "") == 0)
	{
		if (!infPathSpec)
		{
			failmessage.append("-iF or -iPath must be specified.\n");
			fail = true;
		}
	}
	
	if (strcmp(infR.c_str(), "") == 0)
	{
		if (!infPathSpec)
		{
			failmessage.append("-iR or -iPath must be specified.\n");
			fail = true;
		}
	}
		
	if (strcmp(outI.c_str(), "") == 0)
	{
	    failmessage.append("invalid -oIter.\n");
	    fail = true;
	}
		
	if (strcmp(outP.c_str(), "") == 0)
	{
	    failmessage.append("invalid -oPar.\n");
	    fail = true;
	}

	if (strcmp(outCorr.c_str(), "") == 0)
	{
	    failmessage.append("invalid -oCorr.\n");
	    fail = true;
	}

	if (strcmp(outCnt.c_str(), "") == 0)
	{
	    failmessage.append("invalid -oCnt.\n");
	    fail = true;
	}
	
	if (strcmp(outPath.c_str(), "") != 0)
	{
	    DIR* d = opendir(outPath.c_str());
	    if(d)
		{
			closedir(d);
		}
	    else
		{
			failmessage.append("-oPath does not exist.\n");
			fail = true;
		}
	}

	int infFCnt = 1;
	if (infPathSpec)
	{
		infFCnt = countFiles(infPath);
		if (infFCnt < 1)
		{
			failmessage.append("ERROR: infile path \"");
			failmessage.append(infPath.c_str());
			failmessage.append("\" does not contain a positive number of F and R binary files, aborting.\n");
			fail = true;
		}
	}
	
	if (fail)
	{
		cerr << endl << failmessage.c_str() << endl << errorLine << endl;
		return(-1);
	}
	
	if(strcmp(outStub.c_str(),"") != 0)
	{
	    outI = outStub + "_" + outI;
	    outP = outStub + "_" + outP;
	    outCorr = outStub + "_" + outCorr;
		outCnt = outStub + "_" + outCnt;
	}
	
	if(strcmp(outPath.c_str(),"") != 0)
	{
	    outI = outPath + outI;
	    outP = outPath + outP;
	    outCorr = outPath + outCorr;
		outCnt = outPath + outCnt;
	}

	if (seed < -1)
		seed = -1;

	ifstream iff[infFCnt];
	ifstream ifr[infFCnt];
	string fileNames[infFCnt];
	
	if (infPathSpec)
	{
		if (openFiles(infPath, iff, ifr, fileNames) != infFCnt)
		{
			failmessage.append("ERROR: all files in \"");
			failmessage.append(infPath.c_str());
			failmessage.append("\" could not be opened, aborting.\n");
			fail = true;
		}
	}
	else
	{
		iff[0].open(infF.c_str(),ios::binary);
		ifr[0].open(infR.c_str(),ios::binary);
		if (iff[0].fail())
		{
			failmessage.append("ERROR: Forward reads binary file \"");
			failmessage.append(infF.c_str());
			failmessage.append("\" could not be opened, aborting.\n");
			fail = true;
		}
		if (ifr[0].fail())
		{
			failmessage.append("ERROR: Reverse reads binary file \"");
			failmessage.append(infR.c_str());
			failmessage.append("\" could not be opened, aborting.\n");
			fail = true;
		}
	}
	
/*
	iff.open(infF.c_str(),ios::binary);
	if (iff.fail())
	{
		failmessage.append("ERROR: Forward reads binary file \"");
		failmessage.append(infF.c_str());
		failmessage.append("\" could not be opened, aborting.\n");
		fail = true;
	}
	ifstream ifr;
	ifr.open(infR.c_str(),ios::binary);
	if (ifr.fail())
	{
		failmessage.append("ERROR: Reverse reads binary file \"");
		failmessage.append(infR.c_str());
		failmessage.append("\" could not be opened, aborting.\n");
		fail = true;
	}
*/
	
	ofstream ofi;
	ofi.open(outI.c_str(),ios::trunc);
	if (ofi.fail())
	{
		failmessage.append("ERROR: Output file \"");
		failmessage.append(outI.c_str());
		failmessage.append("\" could not be created, aborting.\n");
		fail = true;
	}

	ofstream ofc;
	ofc.open(outCorr.c_str(),ios::trunc);
	if (ofc.fail())
	{
		failmessage.append("ERROR: Output file \"");
		failmessage.append(outCorr.c_str());
		failmessage.append("\" could not be created, aborting.\n");
		fail = true;
	}

	ofstream ofcnt;
	ofcnt.open(outCnt.c_str(),ios::trunc);
	if (ofcnt.fail())
	{
		failmessage.append("ERROR: Output file \"");
		failmessage.append(outCnt.c_str());
		failmessage.append("\" could not be created, aborting.\n");
		fail = true;
	}
	ofstream ofp;
	int truncValF,truncValR;
	int estMinD = minD, estMaxD = maxD, estLambdaD = lambdaD;
	double* distDens;
	vector<string> distDensV;

	ofp.open(outP.c_str(),ios::trunc);
	if (ofp.fail())
	{
		failmessage.append("ERROR: Paramater file \"");
		failmessage.append(outP.c_str());
		failmessage.append("\" could not be created, aborting.\n");
		fail = true;
	}

	if (fail)
	{
		cerr << endl << failmessage.c_str() << endl << errorLine << endl;
		return(-1);
	}
	
	ofi << "! Command:";
	for (int i = 0; i < argc; i++)
		ofi << " " << argv[i];
	ofi << endl;
	
	vlevel = verbose;
	
	/*
	 * Check file lengths
	 */

	int minPos[infFCnt], maxPos[infFCnt], minF[infFCnt], maxF[infFCnt], minR[infFCnt], maxR[infFCnt];
	checkFileLengths(iff, ifr, minPos, maxPos, minR, minF, maxF, maxR, infFCnt);
	
	/*
	 *  Estimate parameters
	 */
	
	ofp << "! Command:";
	for (int i = 0; i < argc; i++)
		ofp << " " << argv[i];
	ofp << endl;
	
	vcerr(1) << "*** Identifying truncation limits and data distributions ***" << endl;
	calculateTruncVal(iff,ifr,&ofcnt,
					  &truncValF, &truncValR, truncLimit,
					  minPos, maxPos,
					  minR, minF,
					  maxF, maxR, infFCnt);
	ofcnt.close();
	
	distDens = new double[maxD-minD+1];

	if (infFCnt == 1)
	{
		estimateDistance(&iff[0],&ifr[0],&ofc,
						 &estMinD, &estMaxD, &estLambdaD,
						 minD, maxD,
						 minPos[0], maxPos[0],
						 minR[0], minF[0],
						 maxF[0], maxR[0],
						 truncValF, truncValR,distDens,0.3);
		ofc.close();
	}
	
	/*
	if (useemp)
	{
		vcerr(2) << "\t* Estimating forward-reverse distance" << endl;
		estimateDistance(&iff,&ifr,&ofc,
						 &estMinD, &estMaxD, &estLambdaD,
						 minD, maxD,
						 minPos, maxPos,
						 minR, minF,
						 maxF, maxR,
						 truncValF, truncValR,distDens,0.3);
		
		ofc.close();
	}
	else
	{
	*/
	for (int dist = minD; dist <= maxD; dist++)
		distDens[dist-minD] = gsl_ran_poisson_pdf(dist, lambdaD);
	/*
	  }
	*/
	
	vcerr(1) << "*** Parameter estimation ***" << endl;
	if (estimateParameters(iff,ifr,&ofi,
						   estMinD,estMaxD,estLambdaD,
						   pF, pR, lambdaF, lambdaR,
						   sampleSize,burnins,iterations,
						   seed,
						   minPos, maxPos,
						   minR, minF,
						   maxF, maxR,
						   truncValF, truncValR, infFCnt) < 0)
	{
		cerr << "ERROR: estimateParameters failed, aborting." << endl;
		
		for (int i=0; i<infFCnt; i++)
		{
			iff[i].close();
			ifr[i].close();
		}		
		ofi.close();
		ofp.close();
		delete[] distDens;
		return(-1);
	}
 	else
 	{
 		if (writeParameters(&ofp, pF, pR, lambdaF, lambdaR, lambdaD,
							truncValF, truncValR,
							minD, maxD, estLambdaD,
							estMinD, estMaxD, distDens) < 0)
 		{
 			cerr << "WARNING: could not write parameters to file, skipping." << endl;
 		}
 	}
	ofi.close();
	ofp.close();

	for (int i=0; i<infFCnt; i++)
	{
		iff[i].close();
		ifr[i].close();
	}

	delete[] distDens;
	
	vcerr(3) << setprecision(3);
	vcerr(3) << endl << "\t\tParameter estimates:" << endl;
	vcerr(3) << "\t\tpF: S = " << pF[0] << " NotS = " << pF[1] << endl;
	vcerr(3) << "\t\tpR: E = " << pR[0] << " NotE = " << pR[1] << endl;
	vcerr(3) << "\t\tlambdaF: S = " << lambdaF[0] << " NotS = " << lambdaF[1] << endl;
	vcerr(3) << "\t\tlambdaR: E = " << lambdaR[0] << " NotE = " << lambdaR[1] << endl;
	vcerr(3) << "\t\testMinD:  = " << estMinD << " estMaxd = " << estMaxD << endl;
}
