
#include "np_utils.h"
#include <dirent.h>

int checkFileLengths(ifstream* infF, ifstream* infR,
		     int* minPos, int* maxPos,
		     int* minR, int* minF,
		     int* maxF, int* maxR)
{
  cerr << "*** Checking file lengths ***" << endl;
  
  getMinMaxF(infF, minF, maxF);
  getMinMaxF(infR, minR, maxR);
  
  if ((*minF != *minR) || (*maxF != *maxR)) // Input files have not the same coordinates
    {
      cerr << "\tFiles start at different positions. Largest common region: ";
      
      *minPos = (*minF > *minR ? *minF : *minR);
      *maxPos = (*maxF < *maxR ? *maxF : *maxR);
      
      cerr << "[ " << *minPos << ", " << *maxPos << " ]." << endl;
    }
  else // input files have the exact same coordinates
    {
      *minPos = *minF;
      *maxPos = *maxF;
    }
  
  cerr<<"*** Done ***"<<endl;
  
  return(1);
}




int detectP(ifstream* infF, ifstream* infR, ofstream* outf,
	    int lag,int minPos, int maxPos,int minF,int minR)
{
  int chunk = CHUNK_MULTI*MAXNR;
  // if necessary, reduce to file size.
  chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);
  
  storageType* tmpF = new storageType[chunk];
  storageType* tmpR = new storageType[chunk];
	
  int read = 0;


  while (((maxPos - minPos + 1) - read >= chunk) && (chunk - lag >= 0))
    {
      cerr << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)read / (double)(maxPos - minPos + 1) << "  \t% complete.\r";
      
      // Read data from F file
      infF->seekg(sizeof(storageType) * (read + (minPos-minF)) + sizeof(int), ios::beg);
      infF->read((char*)tmpF, sizeof(storageType) * chunk);
      // Read data from R file
      infR->seekg(sizeof(storageType) * (read + (minPos-minR)) + sizeof(int), ios::beg);
      infR->read((char*)tmpR, sizeof(storageType) * chunk);
      
      //Loop all F positions
      bool fabove = false;
      bool rabove = false;
      for (int i = 0; i <= chunk - lag; i++)
	{
	  fabove = (int)tmpF[i] > 0;
	  rabove = (int)tmpR[i+lag-1] > 0;
	  if (fabove && rabove)
	    {
	      *(outf)<<i+read+minPos<< "\t" << tmpF[i] << "\t" << tmpR[i+lag-1]<<endl;
	    }
	}
      read += chunk-lag+1;
		
      // Check if we can read in a full chunk of data or simply the rest of the file. 
      if (((maxPos - minPos + 1 - read) < chunk) && ((maxPos - minPos + 1 - read) >= lag))
	chunk = (maxPos - minPos + 1) - read;
      
    }
  cerr << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;
  delete[] tmpF;
  delete[] tmpR;
  return(1);
}



int main(int argc, char* argv[]) 
{
  // parameters that you can set. 
  string ifF     = "";
  string ifR     = "";
  string ouf     = "";
  int lag = 50;

  /*
   *  Read and check input parameters
   */
  
  string errorLine =  "usage " + 
    string(argv[0]) + 
    " [parameters] \n" +
    "\t-iF      <input binary file> \n" + 
    "\t-iR      <input binary file> \n" + 
    "\t-o      <output file> \n" + 
    "\t-l      <lag to check, def 50> \n";
  
  bool fail = false;
  string failmessage = "";
  
  for (int i = 1; i < argc; i++)
    {
      if(strcmp(argv[i],"-iF") == 0)
	ifF.assign(argv[++i]);
      else if(strcmp(argv[i],"-iR") == 0)
	ifR.assign(argv[++i]);
      else if (strcmp(argv[i],"-o") == 0)
	ouf.assign(argv[++i]);
      else if (strcmp(argv[i],"-l") == 0)
	lag = atoi(argv[++i]);
      else
	{
	  failmessage.assign("Unknown argument: ");
	  failmessage.append(argv[i]);
	  failmessage.append("\n");
	  fail = true;
	}
    }
	
  if (strcmp(ifF.c_str(), "") == 0)
    {
      failmessage.append("-iF must be specified.\n");
      fail = true;
    }
 
 if (strcmp(ifR.c_str(), "") == 0)
   {
     failmessage.append("-iR must be specified.\n");
     fail = true;
   }

  if (strcmp(ouf.c_str(), "") == 0)
    {
      failmessage.append("-o must be specified.\n");
      fail = true;
    }

  if (lag < 0)
    {
      failmessage.append("negative lags are pointless.\n");
      fail = true;
    }
  
  ifstream iff;
  iff.open(ifF.c_str(),ios::binary);
  if (iff.fail())
    {
      failmessage.append("ERROR: input binary file \"");
      failmessage.append(ifF.c_str());
      failmessage.append("\" could not be opened, aborting.\n");
      fail = true;
    }

  ifstream ifr;
  ifr.open(ifR.c_str(),ios::binary);
  if (ifr.fail())
    {
      failmessage.append("ERROR: input binary file \"");
      failmessage.append(ifR.c_str());
      failmessage.append("\" could not be opened, aborting.\n");
      fail = true;
    }
  
  ofstream ofp;
  ofp.open(ouf.c_str(),ios::trunc);
  if (ofp.fail())
    {
      failmessage.append("ERROR: Output file \"");
      failmessage.append(ouf.c_str());
      failmessage.append("\" could not be created, aborting.\n");
      fail = true;
    }
  
  if (fail)
    {
      cerr << endl << failmessage.c_str() << endl << errorLine << endl;
      return(-1);
    }
  
  
  /*
   * Check file lengths
   */
  int minPos, maxPos, minF, maxF, minR, maxR;
  checkFileLengths(&iff, &ifr, &minPos, &maxPos, &minR, &minF, &maxF, &maxR);
  
  detectP(&iff,&ifr,&ofp,lag,minPos,maxPos,minF,minR);

  iff.close();
  ifr.close();
  ofp.close();
  return(0);
}
