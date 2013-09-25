
#include "np_utils.h"
#include <dirent.h>

int vlevel = 3;
int main(int argc, char* argv[]) 
{
	// parameters that you can set. 
	string infF    = "";
	string infR    = "";
	string outO    = "oOdds.bin";
	string outI    = "oIter.txt";
	string outP    = "oPar.txt";
	string outB    = "oPos.txt";
	string outCorr = "oCorr.txt";
	string outCnt  = "oCnt.txt";
	string outPath = "";
	string outStub = "";
	bool readest   = false;
	bool onlyest   = false;
	bool useemp    = false;
	double truncLimit = 0.95;
	int verbose    = 3;
	int lambdaD    = 147;
	int minD       = 130;
	int maxD       = 180;
	int startCrd   = 0;
	int endCrd     = -1;
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
	  "\t-oPath   <outdirectory, where all output files are put. \n" + 
	  "\t          NOT created - needs to exists. Defaults to where the program is executed from.>\n" +
	  "\t-oStub   <outfile names stub, \n" +
	  "\t          defaults -oOdds to <-oStub_>oOdds.bin, \n" + 
	  "\t          defaults -oIter to <-oStub_>oIter.txt, \n" + 
	  "\t          defaults -oPos to <-oStub_>oPos.txt, \n" + 
	  "\t          defaults -oPar to <-oStub_>oPar.txt> \n" + 
	  "\t          defaults -oCorr to <-oStub_>oCorr.txt> \n" + 
	  "\t-oOdds   <outfile, where to write posterior odds\n" +
	  "\t          for starts,default generated from -oStub> \n" +
	  "\t-oIter   <outfile, where to write MCMC parameter\n" +
	  "\t          estimates for each iteration,default generated from -oStub> \n" +
	  "\t-oPos    <outfile, where to write predicted\n" +
	  "\t          positions and fuzzyness of boundaries,default generated from -oStub> \n" +
	  "\t-oPar    <parameter-file, where to write/read MCMC\n" +
	  "\t          final parameter estimates,default generated from -oStub> \n" +
	  "\t-oCorr    <out-file, lists correlation coefficients between forward and reverse reads in [mind,maxd]" +
	  "\t          default generated from -oStub> \n" +
	  "\t-oCnt    <out-file, lists forward and reverse read count frequencies" +
	  "\t          default generated from -oStub> \n" +
	  "\t-readest <toggles skip estimation on or off (default off),\n" +
	  "\t          if on then read parameters from\n" +
	  "\t          parameter-file instead.> \n" +
	  "\t-trunc   <truncation limit in (0,1], default 0.99.> \n" +
	  "\t-size    <sample size, default all observations.> \n" +
	  "\t-burnin  <number of MCMC burnin iterations,\n" +
	  "\t          default 1000.> \n" +
	  "\t-iter    <number of MCMC iterations\n" +
	  "\t          (non-burnin iterations), default 1000.> \n" +
	  "\t-seed    <set seed to random number generator.> \n" +
	  "\t-ld      <distance-lambda, default 147 bp.> \n" +
	  "\t-mind    <min distance, default 130 bp.> \n" +
	  "\t-maxd    <max distance, default 180 bp.> \n" +
	  "\t-sc      <limit boundary prediction from this point forward.> \n" +
	  "\t-ec      <limit boundary prediction up until this point.> \n" +
	  "\t-onlyest <toggle only estimate parameters yes/no. default no.>\n" +
	  "\t-useemp  <toggle use of data-driven distance distribution, or poisson around distance-lambda. default off.>\n" + 
	  "\t-v       <verbose level 0-3. default 2>";
	
	bool fail = false;
	string failmessage = "";
	
	for (int i = 1; i < argc; i++)
	  {
	    if(strcmp(argv[i],"-iF") == 0)
	      infF.assign(argv[++i]);
	    else if(strcmp(argv[i],"-iR") == 0)
	      infR.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oPath") == 0)
	      outPath.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oStub") == 0)
	      outStub.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oOdds") == 0)
	      outO.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oIter") == 0)
	      outI.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oPos") == 0)
	      outB.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oPar") == 0)
	      outP.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oCorr") == 0)
	      outCorr.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oCnt") == 0)
	      outCnt.assign(argv[++i]);
	    else if (strcmp(argv[i],"-readest") == 0)
	      readest = true;
	    else if (strcmp(argv[i],"-onlyest") == 0)
	      onlyest = true;
	    else if (strcmp(argv[i],"-useemp") == 0)
	      useemp = true;
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
	    else if (strcmp(argv[i],"-sc") == 0)
	      startCrd = atoi(argv[++i]);
	    else if (strcmp(argv[i],"-ec") == 0)
	      endCrd = atoi(argv[++i]);
	    else
	      {
		failmessage.assign("Unknown argument: ");
		failmessage.append(argv[i]);
		failmessage.append("\n");
		fail = true;
	      }
	  }
	
	if (startCrd > endCrd && endCrd != -1)
	  {
	    failmessage.append("-sc and -ec values do not make sense.\n");
	    fail = true;
	  }

	if (truncLimit <= 0 || truncLimit > 1)
	{
	    failmessage.append("-trunc value does not make sense.\n");
	    fail = true;
	}
	
	// Change so that instead of -oOdds, -oIter, -oPos, -oPar,
	// a -oPath could be specified with file names generated from the
	// names of -iF and -iR
	
	if (strcmp(infF.c_str(), "") == 0)
	  {
	    failmessage.append("-iF must be specified.\n");
	    fail = true;
	  }
	
	if (strcmp(infR.c_str(), "") == 0)
	  {
	    failmessage.append("-iR must be specified.\n");
	    fail = true;
	  }
	
	if (strcmp(outO.c_str(), "") == 0)
	  {
	    failmessage.append("invalid -oOdds.\n");
	    fail = true;
	  }
	
	if (strcmp(outI.c_str(), "") == 0)
	  {
	    failmessage.append("invalid -oIter.\n");
	    fail = true;
	  }
	
	if (strcmp(outB.c_str(), "") == 0)
	  {
	    failmessage.append("invalid -oPos.\n");
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
	    DIR *d;
	    d = opendir(outPath.c_str());
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
	
	if (fail)
	{
		cerr << endl << failmessage.c_str() << endl << errorLine << endl;
		return(-1);
	}
	
	if(strcmp(outStub.c_str(),"") != 0)
	  {
	    outO = outStub + "_" + outO;
	    outI = outStub + "_" + outI;
	    outP = outStub + "_" + outP;
	    outB = outStub + "_" + outB;
	    outCorr = outStub + "_" + outCorr;
		outCnt = outStub + "_" + outCnt;
	  }
	
	if(strcmp(outPath.c_str(),"") != 0)
	  {
	    outO = outPath + outO;
	    outI = outPath + outI;
	    outP = outPath + outP;
	    outB = outPath + outB;
	    outCorr = outPath + outCorr;
		outCnt = outPath + outCnt;
	  }

	if (seed < -1)
	  seed = -1;
	
	ifstream iff;
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
	ofstream ofi;
	ofi.open(outI.c_str(),ios::trunc);
	if (ofi.fail())
	{
		failmessage.append("ERROR: Output file \"");
		failmessage.append(outI.c_str());
		failmessage.append("\" could not be created, aborting.\n");
		fail = true;
	}
	ofstream ofo;
	ofo.open(outO.c_str(),ios::trunc | ios::binary);
	if (ofo.fail())
	{
		failmessage.append("ERROR: Output file \"");
		failmessage.append(outO.c_str());
		failmessage.append("\" could not be created, aborting.\n");
		fail = true;
	}
	ofstream ofb;
	ofb.open(outB.c_str(),ios::trunc);
	if (ofb.fail())
	{
		failmessage.append("ERROR: Output file \"");
		failmessage.append(outB.c_str());
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
	ifstream ifp;
	ofstream ofp;
	int truncValF,truncValR;
	int estMinD = minD, estMaxD = maxD, estLambdaD = lambdaD;
	double* distDens;
	vector<string> distDensV;

	if (readest)
	{
		ifp.open(outP.c_str(),ios::in);
		if (ifp.fail())
		{
			failmessage.append("ERROR: Parameter file \"");
			failmessage.append(outP.c_str());
			failmessage.append("\" could not be opened, aborting.\n");
			fail = true;
		}
		else
		{
		  if(readParameters(&ifp,pF,pR,lambdaF,lambdaR,&lambdaD,&minD,&maxD,
				    &estLambdaD,&estMinD,&estMaxD,&distDensV)<0)			  
		    //if (readParameters(&ifp, pF, pR, lambdaF, lambdaR) < 0)
		    {				
		      failmessage.append("ERROR: Unable to parse parameter file \"");
		      failmessage.append(outP.c_str());
		      failmessage.append("\", aborting.\n");
		      fail = true;
		    }

		  distDens = new double[distDensV.size()];
		  for (int i = 0;i<(int)distDensV.size();i++)
			  distDens[i] = atof(distDensV[i].c_str());		  
		}
	}
	else
	{
		ofp.open(outP.c_str(),ios::trunc);
		if (ofp.fail())
		{
			failmessage.append("ERROR: Paramater file \"");
			failmessage.append(outP.c_str());
			failmessage.append("\" could not be created, aborting.\n");
			fail = true;
		}
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
	int minPos, maxPos, minF, maxF, minR, maxR;
	checkFileLengths(&iff, &ifr, &minPos, &maxPos, &minR, &minF, &maxF, &maxR);
	
	/*
	 *  Estimate parameters if not read from file
	 */
	
	if (!readest)
	{
	  vcerr(1) << "*** Parameter estimation ***" << endl;
	  
	  ofp << "! Command:";
	  for (int i = 0; i < argc; i++)
	    ofp << " " << argv[i];
	  ofp << endl;
	  
	  vcerr(2) << "\t* Identifying truncation limits" << endl;
	  calculateTruncVal(&iff,&ifr,&ofcnt,
			    &truncValF, &truncValR, truncLimit,
			    minPos, maxPos,
			    minR, minF,
			    maxF, maxR);
	  ofcnt.close();

	  distDens = new double[maxD-minD+1];
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
		  for (int dist = minD; dist <= maxD; dist++)
			  distDens[dist-minD] = gsl_ran_poisson_pdf(dist, lambdaD);
	  }
	  
	  if (estimateParameters(&iff,&ifr,&ofi,
				 estMinD,estMaxD,estLambdaD,
				 pF, pR, lambdaF, lambdaR,
				 sampleSize,burnins,iterations,
				 seed,
				 minPos, maxPos,
				 minR, minF,
				 maxF, maxR,
				 truncValF, truncValR) < 0)
	    {
	      cerr << "ERROR: estimateParameters failed, aborting." << endl;
	      
	      iff.close();
	      ifr.close();
	      ofi.close();
	      if (readest)
		ifp.close();
	      else
		ofp.close();
	      
	      return(-1);
	    }
	  else
	    {
	      
	      if (writeParameters(&ofp, pF, pR, lambdaF, lambdaR,lambdaD,minD,maxD,estLambdaD,estMinD,estMaxD,distDens) < 0)
		{
		  cerr << "WARNING: could not write parameters to file, skipping." << endl;
		}
	    }
	}
	ofi.close();
	
	if (readest)
	  ifp.close();
	else
	  ofp.close();
	
	vcerr(3) << setprecision(3);
	vcerr(3) << endl << "\t\tParameter estimates:" << endl;
	vcerr(3) << "\t\tpF: S = " << pF[0] << " NotS = " << pF[1] << endl;
	vcerr(3) << "\t\tpR: E = " << pR[0] << " NotE = " << pR[1] << endl;
	vcerr(3) << "\t\tlambdaF: S = " << lambdaF[0] << " NotS = " << lambdaF[1] << endl;
	vcerr(3) << "\t\tlambdaR: E = " << lambdaR[0] << " NotE = " << lambdaR[1] << endl;
	vcerr(3) << "\t\testMinD:  = " << estMinD << " estMaxd = " << estMaxD << endl;

	// Truncate distDens to 0 outside estMinD and estMaxD
	for (int dist = minD; dist < estMinD; dist++)
		distDens[dist-minD] = 0.0;

	for (int dist = estMaxD+1; dist <= maxD; dist++)
		distDens[dist-minD] = 0.0;
	
	if(!onlyest)
	  {
	    
	    /*
	     *  Calculate odds of center positions
	     */
	    
	    int maxValF;
	    int maxValR;
	    
	    vcerr(1) << "*** Center position odds ***" << endl;
	    if (centerPositionOdds(&iff, &ifr, &ofo,
							   pF, pR, lambdaF, lambdaR,
							   estMinD, estMaxD, &maxValF, &maxValR,
							   startCrd, endCrd,
							   minPos, maxPos,
							   minR, minF,
							   maxF, maxR,
							   truncValF, truncValR) < 0)
	      {
		cerr << "ERROR: centerPositionOdds failed, aborting." << endl;
		
		iff.close();
		ifr.close();
		ofo.close();
		ofb.close();
		
		return(-1);
	      }
	    ofo.close();
		
	    /*
	     *  Predict starts and ends
	     */
	    
	    ifstream ifo;
	    ifo.open(outO.c_str(),ios::in | ios::binary);
	    if (ifo.fail())
	      {
		cerr << "ERROR: Odds file \"" << outO.c_str() << "\" could not be opened, aborting." << endl;
		iff.close();
		ifr.close();
		ofb.close();
		
		return(-1);
	      }
	    
	    vcerr(1) << "*** Boundary prediction ***" << endl;
	    if (predictBoundaries(&iff, &ifr, &ifo, &ofb,
							  pF, pR, lambdaF, lambdaR,
							  estLambdaD, estMinD, estMaxD,
							  minD, maxD,
							  maxValF, maxValR,
							  startCrd,endCrd,
							  minPos, maxPos,
							  minR, minF,
							  maxF, maxR,
							  truncValF, truncValR,
							  distDens) < 0)
	      {
		cerr << "ERROR: predictBoundaries failed, aborting." << endl;
		
		iff.close();
		ifr.close();
		ifo.close();
		ofb.close();
		
		return(-1);
	      }
	  
	    ifo.close();
	  }

	delete[] distDens;
	
	ofb.close();
	iff.close();
	ifr.close();
}
