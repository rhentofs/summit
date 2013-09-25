
#include "np_utils.h"
#include <dirent.h>

int checkFileLengths(ifstream* infO,
		     int* minPos, int* maxPos)
		     
{
  cerr << "*** Checking file lengths ***" << endl;
  
  getMinMaxF(infO, minPos, maxPos);
  
  cerr<<"*** Done ***"<<endl;
  
  return(1);
}

int checkFile(ifstream* infO, ofstream* outf,ofstream* outfd,int minPos, int maxPos,int gap)
{
  int chunk = CHUNK_MULTI*MAXNR;
  // if necessary, reduce to file size.
  chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);
  
  storageType* tmpO = new storageType[chunk];
  
  int read = 0;
  int stretches[1000];
  int maxStretch = 1000;
  for(int i = 0;i<maxStretch;i++)
    stretches[i] = 0;
  bool lastZero = true;
  int strechStart=0,strechEnd=0;
  int gapCnt = 0;
  
  while (((maxPos - minPos + 1) - read > chunk))
    {
      cerr << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)read / (double)(maxPos - minPos + 1) << "  \t% complete.\r";
      // Read data from file
      infO->seekg(sizeof(storageType) * (read) + sizeof(int), ios::beg);
      infO->read((char*)tmpO, sizeof(storageType) * chunk);
      
      for (int i = 0; i <= chunk; i++)
	{
	  // putative center? 
	  if(tmpO[i] > 0)
	    {
	      gapCnt = 0;
	      // first in strech?
	      if(lastZero)
		{
		  strechStart = minPos + read;
		  strechEnd = minPos + read;
		  lastZero = false;
		} 
	      else
		{
		  strechEnd = minPos + read;
		}
	    }
	  else
	    {
	      if(gapCnt>gap)
		{
		  gapCnt = 0;
		  if(!lastZero) //strechEnd 
		    {
		      stretches[(strechEnd-strechStart > maxStretch ? maxStretch-1 : strechEnd-strechStart)]++;
		      (*outf)<<strechStart<<"\t"<<strechEnd<<endl;
		    }
		  lastZero = true;
		}
	      else
		{
		  gapCnt++;
		}
	    }
	  read++;
	}
    
       // Check if we can read in a full chunk of data or simply the rest of the file. 
      if (((maxPos - minPos + 1 - read) < chunk))
	chunk = (maxPos - minPos + 1) - read;
    }
  
  if(!lastZero)
    {
      stretches[(strechEnd-strechStart > maxStretch ? maxStretch-1 : strechEnd-strechStart)]++;
      (*outf)<<strechStart<<"\t"<<strechEnd<<endl;
    }
  for(int i = 0;i<maxStretch;i++)
    (*outfd)<<i+1<<"\t"<<stretches[i]<<endl;
  
  cerr << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;
  delete[] tmpO;
  return(1);
}


int main(int argc, char* argv[]) 
{
  // parameters that you can set. 
  string ifO     = "";
  string ouf     = "";
  string oufd     = "";
  int gap = 0;
  
  /*
   *  Read and check input parameters
   */
  
  string errorLine =  "usage " + 
    string(argv[0]) + 
    " [parameters] \n" +
    "\t-i      <input log odds binary file> \n" + 
    "\t-g      <strech gap allowed, default = 0> \n" + 
    "\t-o      <strech output file> \n" +
    "\t-os      <strech stats output file> \n";
  
  bool fail = false;
  string failmessage = "";
  
  for (int i = 1; i < argc; i++)
    {
      if(strcmp(argv[i],"-i") == 0)
	ifO.assign(argv[++i]);
      else if (strcmp(argv[i],"-o") == 0)
	ouf.assign(argv[++i]);
      else if (strcmp(argv[i],"-os") == 0)
	oufd.assign(argv[++i]);
      else if (strcmp(argv[i],"-g") == 0)
	gap = atoi(argv[++i]);
      else
	{
	  failmessage.assign("Unknown argument: ");
	  failmessage.append(argv[i]);
	  failmessage.append("\n");
	  fail = true;
	}
    }
	
  if (strcmp(ifO.c_str(), "") == 0)
    {
      failmessage.append("-i must be specified.\n");
      fail = true;
    }
  if (strcmp(ouf.c_str(), "") == 0)
    {
      failmessage.append("-o must be specified.\n");
      fail = true;
    }

  if (strcmp(oufd.c_str(), "") == 0)
    {
      failmessage.append("-os must be specified.\n");
      fail = true;
    }

  ifstream iff;
  iff.open(ifO.c_str(),ios::binary);
  if (iff.fail())
    {
      failmessage.append("ERROR: input binary file \"");
      failmessage.append(ifO.c_str());
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

  ofstream ofd;
  ofd.open(oufd.c_str(),ios::trunc);
  if (ofd.fail())
    {
      failmessage.append("ERROR: Output file \"");
      failmessage.append(oufd.c_str());
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
  int minPos, maxPos;
  checkFileLengths(&iff,&minPos, &maxPos);
  
  checkFile(&iff,&ofp,&ofd,minPos,maxPos,gap);

  iff.close();
  ofp.close();
  ofd.close();
  return(0);
}
