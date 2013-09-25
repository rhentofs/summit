
#include "np_utils.h"
#include <dirent.h>



int main(int argc, char* argv[]) 
{
  cout<<sizeof(int)<<endl;
	// parameters that you can set. 
	string inf     = "";

	/*
	 *  Read and check input parameters
	 */

	string errorLine =  "usage " + 
	  string(argv[0]) + 
	  " [parameters] \n" +
	  "\t-i      <input binary file> \n";

	bool fail = false;
	string failmessage = "";
	
	for (int i = 1; i < argc; i++)
	  {
	    if(strcmp(argv[i],"-i") == 0)
	      inf.assign(argv[++i]);
	    else
	      {
		failmessage.assign("Unknown argument: ");
		failmessage.append(argv[i]);
		failmessage.append("\n");
		fail = true;
	      }
	  }
	
	if (strcmp(inf.c_str(), "") == 0)
	  {
	    failmessage.append("-i must be specified.\n");
	    fail = true;
	  }
	
	ifstream iff;
	iff.open(inf.c_str(),ios::binary);
	if (iff.fail())
	  {
	    failmessage.append("ERROR: input binary file \"");
	    failmessage.append(inf.c_str());
	    failmessage.append("\" could not be opened, aborting.\n");
	    fail = true;
	  }
	
	if (fail)
	{
		cerr << endl << failmessage.c_str() << endl << errorLine << endl;
		return(-1);
	}
	
	
	// get min/max postions of  the input file. 
	int minPos=0,maxPos = 0;
	getMinMaxF(&iff, &minPos, &maxPos);
	cerr<<"Input range:\t"<<(long)minPos<<"\t"<<(long)maxPos<<endl;
	
	// start reading the input
	
	int chunk = CHUNK_MULTI*MAXNR;
	int read = 0;
	double total = 0;

	// if necessary, reduce to file size.
	chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);
	
	storageType* inpData = new storageType[chunk];
	storageType* tmpData = new storageType[1];
	cerr << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete.";
	
	while ((maxPos - minPos + 1) - read >= chunk && (chunk > 0))
	  {
	    // Read data from input file
	    iff.seekg(sizeof(storageType)*read + sizeof(int), ios::beg);
	    iff.read((char*)inpData, sizeof(storageType) * chunk);
	    
	    //Loop all positions
	    for (int i = 0; i < chunk; i++)
	      total += (double)inpData[i];
	    
	    read += chunk;
	    
	    // Check if we can read in a full chunk of data or simply the rest of the file. 
	    if (((maxPos - minPos + 1 - read) < chunk))
	      chunk = (maxPos - minPos + 1) - read;
	    
	    cerr << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)read / (double)(maxPos - minPos + 1) << "  \t% complete.\r";
	  }
	cerr << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;
	
	delete[] inpData;
	delete[] tmpData;
	
	iff.close();
	cerr<<(long)total<<endl;
	cout<<(long)total<<endl;
}
