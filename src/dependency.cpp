#include "zeb_utils.h"

int vlevel = 3;
int main(int argc, char* argv[]) 
{
	// parameters that you can set. 
	string infF1    = "";
	string infF2    = "";
	string infR1    = "";
	string infR2    = "";
	string infPath1 = "";
	string infPath2 = "";
	string infPathOdds1 = "";
	string infPathOdds2 = "";
	string infP1    = "iPar.txt";
	string infP2    = "iPar.txt";
	string outD     = "oDep.txt";
	string outPath  = "";
	string outStub  = "";
	int sampleSize  = 100000;
	int verbose     = 3;

	// Read from parameter files
	int lambdaD1    = 147;
	int lambdaD2    = 147;
	int minD1       = 130;
	int minD2       = 130;
	int maxD1       = 180;
	int maxD2       = 180;
	double lambdaF1[2];
	double lambdaF2[2];
	double lambdaR1[2];
	double lambdaR2[2];
	double pF1[2];
	double pF2[2];
	double pR1[2];
	double pR2[2];

	/*
	 *  Read and check input parameters
	 */
	
	string errorLine =  "usage " + 
		string(argv[0]) + 
		" [parameters] \n" +
		"\t-iF1      <forward reads binary file 1> \n" +
		"\t-iF2      <forward reads binary file 2> \n" +
		"\t-iR1      <reverse reads binary file 1> \n" +
		"\t-iR2      <reverse reads binary file 2> \n" +
		"\t-iPath1   <path 1 to forward and reverse reads binary files \n" +
		"\t           required to have names like F_<chrom>.bin or R_<chrom>.bin \n" +
		"             (overrides iR1 and iF1 if set)> \n" +
		"\t-iPath2   <path 2 to forward and reverse reads binary files \n" +
		"\t           required have names like F_<chrom>.bin or R_<chrom>.bin \n" +
		"             (overrides iR2 and iF2 if set)> \n" +
		"\t-iPathOdds1   <path 1 to start and end binary odds files> \n" +
		"\t-iPathOdds2   <path 2 to start and end binary odds files> \n" +
		"\t-oPath   <outdirectory, where all output files are put. \n" + 
		"\t          NOT created - needs to exists. Defaults to where the program is executed from.>\n" +
		"\t-oStub   <outfile names stub, \n" +
		"\t          defaults -oDep to <-oStub_>oDep.txt> \n" +
		"\t-oDep    <outfile, where to write regression variables,\n" +
		"\t          default generated from -oStub> \n" +
		"\t-iPar1   <parameter-file 1, where to read MCMC\n" +
		"\t          final parameter estimates> \n" +
		"\t-iPar2   <parameter-file 2, where to read MCMC\n" +
		"\t          final parameter estimates> \n" +
		"\t-size    <sample size. default 100000>";
		"\t-v       <verbose level 0-3. default 2>";
	
	bool fail = false;
	string failmessage = "";
	
	for (int i = 1; i < argc; i++)
	{
		if(strcmp(argv[i],"-iF1") == 0)
			infF1.assign(argv[++i]);
		else if(strcmp(argv[i],"-iF2") == 0)
			infF2.assign(argv[++i]);
	    else if(strcmp(argv[i],"-iR1") == 0)
			infR1.assign(argv[++i]);
	    else if(strcmp(argv[i],"-iR2") == 0)
			infR2.assign(argv[++i]);
		else if(strcmp(argv[i],"-iPath1") == 0)
			infPath1.assign(argv[++i]);
		else if(strcmp(argv[i],"-iPath2") == 0)
			infPath2.assign(argv[++i]);
		else if(strcmp(argv[i],"-iPathOdds1") == 0)
			infPathOdds1.assign(argv[++i]);
		else if(strcmp(argv[i],"-iPathOdds2") == 0)
			infPathOdds2.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oPath") == 0)
			outPath.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oStub") == 0)
			outStub.assign(argv[++i]);
	    else if(strcmp(argv[i],"-oDep") == 0)
			outD.assign(argv[++i]);
	    else if(strcmp(argv[i],"-iPar1") == 0)
			infP1.assign(argv[++i]);
	    else if(strcmp(argv[i],"-iPar2") == 0)
			infP2.assign(argv[++i]);
	    else if (strcmp(argv[i],"-size") == 0)
			sampleSize = atoi(argv[++i]);
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

	bool infPath1Spec = false;
	bool infPath2Spec = false;
	if (strcmp(infPath1.c_str(), "") != 0)
	{
		infPath1Spec = true;
	    DIR *d = opendir(infPath1.c_str());
	    if(d)
		{
			closedir(d);
		}
	    else
		{
			failmessage.append("-iPath1 does not exist.\n");
			fail = true;
		}
	}
	if (strcmp(infPath2.c_str(), "") != 0)
	{
		infPath2Spec = true;
	    DIR *d = opendir(infPath2.c_str());
	    if(d)
		{
			closedir(d);
		}
	    else
		{
			failmessage.append("-iPath2 does not exist.\n");
			fail = true;
		}
	}	

	if ((infPath1Spec && !infPath2Spec) ||
		(!infPath1Spec && infPath2Spec))
	{
		failmessage.append("ERROR: if infile paths are used, they must be specified for both samples, aborting.\n");
		fail = true;
	}
	
	if (strcmp(infF1.c_str(), "") == 0)
	{
		if (!infPath1Spec)
		{
			failmessage.append("-iF1 must be specified.\n");
			fail = true;
		}
	}	
	if (strcmp(infR1.c_str(), "") == 0)
	{
		if (!infPath1Spec)
		{
			failmessage.append("-iR1 must be specified.\n");
			fail = true;
		}
	}

	if (strcmp(infF2.c_str(), "") == 0)
	{
		if (!infPath2Spec)
		{
			failmessage.append("-iF2 must be specified.\n");
			fail = true;
		}
	}	
	if (strcmp(infR2.c_str(), "") == 0)
	{
		if (!infPath2Spec)
		{
			failmessage.append("-iR2 must be specified.\n");
			fail = true;
		}
	}
	
	if (strcmp(outD.c_str(), "") == 0)
	{
	    failmessage.append("invalid -oDep.\n");
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

	int infFCnt1 = 1;
	int infFCnt2 = 1;
	if (infPath1Spec)
	{
		infFCnt1 = countFiles(infPath1);
		if (infFCnt1 < 1)
		{
			failmessage.append("ERROR: infile path \"");
			failmessage.append(infPath1.c_str());
			failmessage.append("\" does not contain a positive number of F and R binary files, aborting.\n");
			fail = true;
		}
	}
	if (infPath2Spec)
	{
		infFCnt2 = countFiles(infPath2);
		if (infFCnt2 < 1)
		{
			failmessage.append("ERROR: infile path \"");
			failmessage.append(infPath2.c_str());
			failmessage.append("\" does not contain a positive number of F and R binary files, aborting.\n");
			fail = true;
		}
	}

	if (infFCnt1 != infFCnt2)
	{
		failmessage.append("ERROR: infile paths \"");
		failmessage.append(infPath1.c_str());
		failmessage.append("\" and \"");
		failmessage.append(infPath2.c_str());
		failmessage.append("\" does not contain an equal number of F and R binary files, aborting.\n");
		fail = true;
	}

	if (fail)
	{
		cerr << endl << failmessage.c_str() << endl << errorLine << endl;
		return(-1);
	}

	ifstream iff1[infFCnt1];
	ifstream iff2[infFCnt1];
	ifstream ifr1[infFCnt1];
	ifstream ifr2[infFCnt1];

	ifstream ofos1[infFCnt1];
	ifstream ofos2[infFCnt1];
	ifstream ofoe1[infFCnt1];
	ifstream ofoe2[infFCnt1];

	string fileNames[infFCnt1];
	
	if (infPath1Spec)
	{
		if (openFiles(infPath1, iff1, ifr1, fileNames) != infFCnt1)
		{
			failmessage.append("ERROR: all files in \"");
			failmessage.append(infPath1.c_str());
			failmessage.append("\" could not be opened, aborting.\n");
			fail = true;
		}
		if (openFilesNames(infPath2, iff2, ifr2, fileNames, infFCnt1) != infFCnt1)
		{
			failmessage.append("ERROR: all files in \"");
			failmessage.append(infPath2.c_str());
			failmessage.append("\" could not be opened, aborting.\n");
			fail = true;
		}
		for (int i = 0; i < infFCnt1; i++)
		{
			string fileS1 = infPathOdds1 + fileNames[i] + "_Start_oOdds.bin";
			string fileS2 = infPathOdds2 + fileNames[i] + "_Start_oOdds.bin";
			string fileE1 = infPathOdds1 + fileNames[i] + "_End_oOdds.bin";
			string fileE2 = infPathOdds2 + fileNames[i] + "_End_oOdds.bin";
			
			ofos1[i].open(fileS1.c_str(),ios::binary);
			ofos2[i].open(fileS2.c_str(),ios::binary);
			ofoe1[i].open(fileE1.c_str(),ios::binary);
			ofoe2[i].open(fileE2.c_str(),ios::binary);
			
			if (ofos1[i].fail() || ofos2[i].fail() || ofoe1[i].fail() || ofoe2[i].fail())
			{
				failmessage.append("ERROR: odds file for chrom \"");
				failmessage.append(fileNames[i].c_str());
				failmessage.append("\" could not be opened, aborting.\n");
				fail = true;
			}
		}
	}
	else
	{
		iff1[0].open(infF1.c_str(),ios::binary);
		ifr1[0].open(infR1.c_str(),ios::binary);
		if (iff1[0].fail())
		{
			failmessage.append("ERROR: Forward reads binary file \"");
			failmessage.append(infF1.c_str());
			failmessage.append("\" could not be opened, aborting.\n");
			fail = true;
		}
		if (ifr1[0].fail())
		{
			failmessage.append("ERROR: Reverse reads binary file \"");
			failmessage.append(infR1.c_str());
			failmessage.append("\" could not be opened, aborting.\n");
			fail = true;
		}
		iff2[0].open(infF2.c_str(),ios::binary);
		ifr2[0].open(infR2.c_str(),ios::binary);
		if (iff2[0].fail())
		{
			failmessage.append("ERROR: Forward reads binary file \"");
			failmessage.append(infF2.c_str());
			failmessage.append("\" could not be opened, aborting.\n");
			fail = true;
		}
		if (ifr2[0].fail())
		{
			failmessage.append("ERROR: Reverse reads binary file \"");
			failmessage.append(infR2.c_str());
			failmessage.append("\" could not be opened, aborting.\n");
			fail = true;
		}
	}

	if(strcmp(outStub.c_str(),"") != 0)
	{
		outD = outStub + "_" + outD;
	}
	
	if(strcmp(outPath.c_str(),"") != 0)
	{
		outD = outPath + outD;
	}
		
	ofstream ofd;
	ofd.open(outD.c_str(),ios::trunc);
	if (ofd.fail())
	{
		failmessage.append("ERROR: Output file \"");
		failmessage.append(outD.c_str());
		failmessage.append("\" could not be created, aborting.\n");
		fail = true;
	}
	
	vlevel = verbose;
	
	ifstream ifp1;
	ifstream ifp2;
	int truncValF1,truncValR1,truncValF2,truncValR2;
	int estMinD1, estMaxD1, estLambdaD1, estMinD2, estMaxD2, estLambdaD2;
	double* distDens1 = NULL;
	double* distDens2 = NULL;
	vector<string> distDensV1;
	vector<string> distDensV2;

	ofd << "! Command:";
	for (int i = 0; i < argc; i++)
		ofd << " " << argv[i];
	ofd << endl;

	if (!fail)
	{
		ifp1.open(infP1.c_str(),ios::in);
		if (ifp1.fail())
		{
			failmessage.append("ERROR: Parameter file \"");
			failmessage.append(infP1.c_str());
			failmessage.append("\" could not be opened, aborting.\n");
			fail = true;
		}
		else
		{
			vcerr(1) << "*** Reading parameters (file 1) ***" << endl;
			if(readParameters(&ifp1,pF1,pR1,lambdaF1,lambdaR1,&lambdaD1,
							  &truncValF1, &truncValR1,
							  &minD1,&maxD1,
							  &estLambdaD1,&estMinD1,&estMaxD1,&distDensV1)<0)			  
			{				
				failmessage.append("ERROR: Unable to parse parameter file \"");
				failmessage.append(infP1.c_str());
				failmessage.append("\", aborting.\n");
				fail = true;
			}
			
			if (!fail)
			{
				distDens1 = new double[distDensV1.size()];
				for (int i = 0;i<(int)distDensV1.size();i++)
					distDens1[i] = atof(distDensV1[i].c_str());
			}
		}
		ifp1.close();
	}

	if (!fail)
	{
		ifp2.open(infP2.c_str(),ios::in);
		if (ifp2.fail())
		{
			failmessage.append("ERROR: Parameter file \"");
			failmessage.append(infP2.c_str());
			failmessage.append("\" could not be opened, aborting.\n");
			fail = true;
		}
		else
		{
			vcerr(1) << "*** Reading parameters (file 2) ***" << endl;
			if(readParameters(&ifp2,pF2,pR2,lambdaF2,lambdaR2,&lambdaD2,
							  &truncValF2, &truncValR2,
							  &minD2,&maxD2,
							  &estLambdaD2,&estMinD2,&estMaxD2,&distDensV2)<0)			  
			{				
				failmessage.append("ERROR: Unable to parse parameter file \"");
				failmessage.append(infP2.c_str());
				failmessage.append("\", aborting.\n");
				fail = true;
			}
			
			if (!fail)
			{
				distDens2 = new double[distDensV2.size()];
				for (int i = 0;i<(int)distDensV2.size();i++)
					distDens2[i] = atof(distDensV2[i].c_str());
			}
		}
		ifp2.close();
	}

	if (!fail)
	{
		if ((minD1 != minD2) ||
			(maxD1 != maxD2) ||
			(lambdaD1 != lambdaD2) ||
			(estMinD1 != estMinD2) ||
			(estMaxD1 != estMaxD2) ||
			(estLambdaD1 != estLambdaD2))
		{
			failmessage.append("ERROR: boundary parameters specified in the two parameter files differ, aborting.");
			fail = true;
		}
	}
	
	if (fail)
	{
		cerr << endl << failmessage.c_str() << endl << errorLine << endl;
		for (int i=0; i<infFCnt1; i++)
		{
			iff1[i].close();
			ifr1[i].close();
			iff2[i].close();
			ifr2[i].close();
			ofos1[i].close();
			ofos2[i].close();
			ofoe1[i].close();
			ofoe2[i].close();
		}
		return(-1);
	}
		
	vcerr(3) << setprecision(3);
	vcerr(3) << "\tParameter estimates (file 1):" << endl;
	vcerr(3) << "\tpF: S = " << pF1[0] << " NotS = " << pF1[1] << endl;
	vcerr(3) << "\tpR: E = " << pR1[0] << " NotE = " << pR1[1] << endl;
	vcerr(3) << "\tlambdaF: S = " << lambdaF1[0] << " NotS = " << lambdaF1[1] << endl;
	vcerr(3) << "\tlambdaR: E = " << lambdaR1[0] << " NotE = " << lambdaR1[1] << endl;
	vcerr(3) << "\testMinD:  = " << estMinD1 << " estMaxd = " << estMaxD1 << endl << endl;
	vcerr(3) << "\tParameter estimates (file 2):" << endl;
	vcerr(3) << "\tpF: S = " << pF2[0] << " NotS = " << pF2[1] << endl;
	vcerr(3) << "\tpR: E = " << pR2[0] << " NotE = " << pR2[1] << endl;
	vcerr(3) << "\tlambdaF: S = " << lambdaF2[0] << " NotS = " << lambdaF2[1] << endl;
	vcerr(3) << "\tlambdaR: E = " << lambdaR2[0] << " NotE = " << lambdaR2[1] << endl;
	vcerr(3) << "\testMinD:  = " << estMinD2 << " estMaxd = " << estMaxD2 << endl;
	
	/*
	 * Check file lengths
	 */
	int minPos[infFCnt1], maxPos[infFCnt1];
	int minF1[infFCnt1], maxF1[infFCnt1], minR1[infFCnt1], maxR1[infFCnt1];
	int minF2[infFCnt2], maxF2[infFCnt2], minR2[infFCnt2], maxR2[infFCnt2];
	checkFileLengths(iff1, iff2, ifr1, ifr2, minPos, maxPos,
					 minR1, minR2, minF1, minF2, maxF1, maxF2, maxR1, maxR2,
					 infFCnt1);
		
	// Truncate distDens to 0 outside estMinD and estMaxD
	for (int dist = minD1; dist < estMinD1; dist++)
	{
		distDens1[dist-minD1] = 0.0;
		distDens2[dist-minD1] = 0.0;
	}

	for (int dist = estMaxD1+1; dist <= maxD1; dist++)
	{
		distDens1[dist-minD1] = 0.0;
		distDens2[dist-minD1] = 0.0;
	}
	
	if (dependency(iff1, iff2, ifr1, ifr2,
				   ofos1, ofos2, ofoe1, ofoe2, &ofd,
				   pF1, pF2, pR1, pR2,
				   lambdaF1, lambdaF2, lambdaR1, lambdaR2,
				   estMinD1, estMaxD1,
				   minPos, maxPos,
				   minR1, minR2, minF1, minF2,
				   maxF1, maxF2, maxR1, maxR2,
				   truncValF1, truncValF2, truncValR1, truncValR2,
				   infFCnt1, fileNames, sampleSize) < 0)
	{
		cerr << "ERROR: failed, aborting." << endl;
		
		for (int i=0; i<infFCnt1; i++)
		{
			iff1[i].close();
			ifr1[i].close();
			iff2[i].close();
			ifr2[i].close();
			ofos1[i].close();
			ofos2[i].close();
			ofoe1[i].close();
			ofoe2[i].close();
		}
		ofd.close();
		delete[] distDens1;
		delete[] distDens2;
		return(-1);
	}
	
	ofd.close();
	
	delete[] distDens1;
	delete[] distDens2;
	
	for (int i=0; i<infFCnt1; i++)
	{
		iff1[i].close();
		ifr1[i].close();
		iff2[i].close();
		ifr2[i].close();
		ofos1[i].close();
		ofos2[i].close();
		ofoe1[i].close();
		ofoe2[i].close();
	}
}
