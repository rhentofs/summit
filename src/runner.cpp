#include "zeb_utils.h"

int vlevel = 3;
int main(int argc, char* argv[]) 
{
	// parameters that you can set. 
	string infF    = "";
	string infR    = "";
	string infPath = "";
	string outO    = "oOdds.bin";
	string outF    = "oFuzzy.bin";
	string infP    = "iPar.txt";
	//string outB    = "oPos.txt";
	string outS    = "oStats.txt";
	string outPath = "";
	string outStub = "";
	int verbose    = 3;
	int lambdaD    = 147;
	int minD       = 130;
	int maxD       = 180;
	// Unknowns
	double lambdaF[2];
	double lambdaR[2];
	double pF[2];
	double pR[2];
	bool bodds = false;
	bool neg = false;
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
		"\t          NOT created - needs to exists. Defaults to where the program is executed from.>\n" +
		"\t-oStub   <outfile names stub, \n" +
		"\t          defaults -oOdds to <-oStub_>oOdds.bin, \n" +
		"\t          defaults -oFuzzy to <-oStub_>oFuzzy.bin \n" + 
		"\t-oOdds   <outfile, where to write posterior odds\n" +
		"\t          for insides, default generated from -oStub> \n" +
		"\t-oFuzzy  <outfile, where to write fuzzyness scores\n" +
		"\t          for starts, default generated from -oStub> \n" +
		"\t-oStats  <outfile, where to write stats from runner.\n" +
		/*
		  "\t-oPos    <outfile, where to write predicted\n" +
		  "\t          positions and fuzzyness of boundaries,default generated from -oStub> \n" +
		*/
		"\t-iPar    <parameter-file, where to read MCMC\n" +
		"\t          final parameter estimates, default generated from -oStub> \n" +
		"\t-bodds   <output boundary odds, default off.> \n" +
		"\t-negodds <output negative odds, default off.> \n" +
		"\t-sc      <limit prediction from this point forward.> \n" +
		"\t-ec      <limit prediction up until this point.> \n" +
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
	    else if(strcmp(argv[i],"-oOdds") == 0)
			outO.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oStats") == 0)
			outS.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oFuzzy") == 0)
			outF.assign(argv[++i]);
	    /*
		  else if(strcmp(argv[i],"-oPos") == 0)
		  outB.assign(argv[++i]);
		*/
	    else if(strcmp(argv[i],"-iPar") == 0)
			infP.assign(argv[++i]);
		else if(strcmp(argv[i],"-bodds") == 0)
			bodds = true;
		else if(strcmp(argv[i],"-negodds") == 0)
			neg = true;
	    else if (strcmp(argv[i],"-v") == 0)
			verbose = atoi(argv[++i]);
	    else
		{
			failmessage.assign("Unknown argument: ");
			failmessage.append(argv[i]);
			failmessage.append("\n");
			fail = true;
		}
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
			failmessage.append("-iF must be specified.\n");
			fail = true;
		}
	}
	
	if (strcmp(infR.c_str(), "") == 0)
	{
		if (!infPathSpec)
		{
			failmessage.append("-iR must be specified.\n");
			fail = true;
		}
	}
	
	if (strcmp(outO.c_str(), "") == 0)
	{
	    failmessage.append("invalid -oOdds.\n");
	    fail = true;
	}

	if (strcmp(outS.c_str(), "") == 0)
	{
	    failmessage.append("invalid -oStats.\n");
	    fail = true;
	}
	
	if (strcmp(outF.c_str(), "") == 0)
	{
	    failmessage.append("invalid -oFuzzy.\n");
	    fail = true;
	}

	/*
	if (strcmp(outB.c_str(), "") == 0)
	{
	    failmessage.append("invalid -oPos.\n");
	    fail = true;
	}
	*/
	
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

	string outOs[infFCnt];
	string outFs[infFCnt];

	string outOSs[infFCnt];
	string outOEs[infFCnt];
	
	if(strcmp(outStub.c_str(),"") != 0)
	{
		if (infPathSpec)
		{
			for (int i=0; i < infFCnt; i++)
			{
				outOs[i] = outStub + "_" + fileNames[i] + "_" + outO;
				if (bodds)
				{
					outOSs[i] = outStub + "_" + fileNames[i] + "_Start_" + outO;
					outOEs[i] = outStub + "_" + fileNames[i] + "_End_" + outO;
				}
				outFs[i] = outStub + "_" + fileNames[i] + "_" + outF;
			}
		}
		else
		{
			outOs[0] = outStub + "_" + outO;
			if (bodds)
			{
				outOSs[0] = outStub + "_Start_" + outO;
				outOEs[0] = outStub + "_End_" + outO;
			}
			outFs[0] = outStub + "_" + outF;
		}
		outS = outStub + "_" + outS;
	    //outB = outStub + "_" + outB;
	}
	else
	{
		if (infPathSpec)
		{
			for (int i=0; i < infFCnt; i++)
			{
				outOs[i] = fileNames[i] + "_" + outO;
				if (bodds)
				{
					outOSs[i] = fileNames[i] + "_Start_" + outO;
					outOEs[i] = fileNames[i] + "_End_" + outO;
				}
				outFs[i] = fileNames[i] + "_" + outF;
			}
		}
		else
		{
			outOs[0] = outO;
			if (bodds)
			{
				outOSs[0] = "Start_" + outO;
				outOEs[0] = "End_" + outO;
			}
			outFs[0] = outF;
		}
	}
	
	if(strcmp(outPath.c_str(),"") != 0)
	{
		for (int i=0; i < infFCnt; i++)
		{
			outOs[i] = outPath + outOs[i];
			if (bodds)
			{
				outOSs[i] = outPath + outOSs[i];
				outOEs[i] = outPath + outOEs[i];
			}
			outFs[i] = outPath + outFs[i];
		}
		outS = outPath + outS;
	    //outB = outPath + outB;
	}
		
	ofstream ofo[infFCnt];
	ofstream off[infFCnt];
	ofstream ofos[infFCnt];
	ofstream ofoe[infFCnt];
	for (int i = 0; i < infFCnt; i++)
	{
		ofo[i].open(outOs[i].c_str(),ios::trunc | ios::binary);
		off[i].open(outFs[i].c_str(),ios::trunc | ios::binary);
		if (ofo[i].fail())
		{
			failmessage.append("ERROR: Output file \"");
			failmessage.append(outOs[i].c_str());
			failmessage.append("\" could not be created, aborting.\n");
			fail = true;
		}
		if (off[i].fail())
		{
			failmessage.append("ERROR: Output file \"");
			failmessage.append(outFs[i].c_str());
			failmessage.append("\" could not be created, aborting.\n");
			fail = true;
		}
		if (bodds)
		{
			ofos[i].open(outOSs[i].c_str(),ios::trunc | ios::binary);
			ofoe[i].open(outOEs[i].c_str(),ios::trunc | ios::binary);
			if (ofos[i].fail())
			{
				failmessage.append("ERROR: Output file \"");
				failmessage.append(outOSs[i].c_str());
				failmessage.append("\" could not be created, aborting.\n");
				fail = true;
			}
			if (ofoe[i].fail())
			{
				failmessage.append("ERROR: Output file \"");
				failmessage.append(outOEs[i].c_str());
				failmessage.append("\" could not be created, aborting.\n");
				fail = true;
			}
		}
	}

	ofstream ofs;
	ofs.open(outS.c_str(),ios::trunc);
	if (ofs.fail())
	{
		failmessage.append("ERROR: Output file \"");
		failmessage.append(outS.c_str());
		failmessage.append("\" could not be created, aborting.\n");
		fail = true;
	}
	
	/*
	ofstream ofb;
	ofb.open(outB.c_str(),ios::trunc);
	if (ofb.fail())
	{
		failmessage.append("ERROR: Output file \"");
		failmessage.append(outB.c_str());
		failmessage.append("\" could not be created, aborting.\n");
		fail = true;
	}
	*/

	vlevel = verbose;
	
	ifstream ifp;
	int truncValF,truncValR;
	int estMinD, estMaxD, estLambdaD;
	double* distDens = NULL;
	vector<string> distDensV;

	ofs << "! Command:";
	for (int i = 0; i < argc; i++)
		ofs << " " << argv[i];
	ofs << endl;
	
	ifp.open(infP.c_str(),ios::in);
	if (ifp.fail())
	{
		failmessage.append("ERROR: Parameter file \"");
		failmessage.append(infP.c_str());
		failmessage.append("\" could not be opened, aborting.\n");
		fail = true;
	}
	else
	{
		vcerr(1) << "*** Reading parameters ***" << endl;
		if(readParameters(&ifp,pF,pR,lambdaF,lambdaR,&lambdaD,
						  &truncValF, &truncValR,
						  &minD,&maxD,
						  &estLambdaD,&estMinD,&estMaxD,&distDensV)<0)			  
		{				
			failmessage.append("ERROR: Unable to parse parameter file \"");
			failmessage.append(infP.c_str());
			failmessage.append("\", aborting.\n");
			fail = true;
		}

		if (!fail)
		{
			distDens = new double[distDensV.size()];
			for (int i = 0;i<(int)distDensV.size();i++)
				distDens[i] = atof(distDensV[i].c_str());
		}
	}

	ifp.close();
	
	if (fail)
	{
		cerr << endl << failmessage.c_str() << endl << errorLine << endl;
		for (int i=0; i<infFCnt; i++)
		{
			iff[i].close();
			ifr[i].close();
		}
		return(-1);
	}
		
	vcerr(3) << setprecision(3);
	vcerr(3) << "\tParameter estimates:" << endl;
	vcerr(3) << "\tpF: S = " << pF[0] << " NotS = " << pF[1] << endl;
	vcerr(3) << "\tpR: E = " << pR[0] << " NotE = " << pR[1] << endl;
	vcerr(3) << "\tlambdaF: S = " << lambdaF[0] << " NotS = " << lambdaF[1] << endl;
	vcerr(3) << "\tlambdaR: E = " << lambdaR[0] << " NotE = " << lambdaR[1] << endl;
	vcerr(3) << "\testMinD:  = " << estMinD << " estMaxd = " << estMaxD << endl;
	
	/*
	 * Check file lengths
	 */
	int minPos[infFCnt], maxPos[infFCnt], minF[infFCnt], maxF[infFCnt], minR[infFCnt], maxR[infFCnt];
	checkFileLengths(iff, ifr, minPos, maxPos, minR, minF, maxF, maxR, infFCnt);
		
	// Truncate distDens to 0 outside estMinD and estMaxD
	for (int dist = minD; dist < estMinD; dist++)
		distDens[dist-minD] = 0.0;

	for (int dist = estMaxD+1; dist <= maxD; dist++)
		distDens[dist-minD] = 0.0;
	
	/*
	 *  Calculate odds of center positions
	 */

	int maxValF;
	int maxValR;

	cerr << "bodds: " << (bodds ? "yes" : "no") << endl;
	if (bodds)
	{
		cerr << outOs[0].c_str() << endl;
		cerr << outOSs[0].c_str() << endl;
		cerr << outOEs[0].c_str() << endl;
	}
	
	vcerr(1) << "*** Center position odds ***" << endl;
	bool ret = false;
	if (neg)
	{
		ret = centerPositionOddsNeg(iff, ifr, ofo, ofos, ofoe, off, &ofs,
									bodds,
									pF, pR, lambdaF, lambdaR,
									estMinD, estMaxD, &maxValF, &maxValR,
									minPos, maxPos,
									minR, minF,
									maxF, maxR,
									truncValF, truncValR,
									infFCnt, fileNames);
	}
	else
	{
		ret = centerPositionOdds(iff, ifr, ofo, ofos, ofoe, off, &ofs,
								 bodds,
								 pF, pR, lambdaF, lambdaR,
								 estMinD, estMaxD, &maxValF, &maxValR,
								 minPos, maxPos,
								 minR, minF,
								 maxF, maxR,
								 truncValF, truncValR,
								 infFCnt, fileNames);
		
	}
	if (ret < 0)
	{
		cerr << "ERROR: centerPositionOdds failed, aborting." << endl;
		
		for (int i=0; i<infFCnt; i++)
		{
			iff[i].close();
			ifr[i].close();
			ofo[i].close();
			off[i].close();
			if (bodds)
			{
				ofos[i].close();
				ofoe[i].close();
			}
		}
		ofs.close();
		//ofb.close();
		delete[] distDens;
		return(-1);
	}
	
	for (int i=0; i<infFCnt; i++)
	{
		ofo[i].close();
		off[i].close();
		if (bodds)
		{
			ofos[i].close();
			ofoe[i].close();
		}
	}
	ofs.close();

	  /*
	  *  Predict starts and ends
	  */

	/*
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
	  ofb.close();
    */
	
	delete[] distDens;
	
	for (int i=0; i<infFCnt; i++)
	{
		iff[i].close();
		ifr[i].close();
	}
}
