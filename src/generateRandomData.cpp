
#include "np_utils.h"
#include <dirent.h>


/*
 * create a empty (full of zeroes) binary file with the given number of positions.
 */

void createEmptyBinary (fstream * fs,int minc,int maxc,bool progress)
{

  int nrPos = maxc-minc;
  int wrPos = 0;
  storageType tmpShort = 0;

  storageType * tmpSignal;
  /*
   * Write the start coordinate offset and fill with zeroes. 
   * Large chunks when possible and then one at a time. 
   */ 
  fs->write((char*)&minc,sizeof(int));
  tmpSignal = new storageType[CHUNK_SIZE];
  for (int i=0;i<CHUNK_SIZE;i++)
    tmpSignal[i] = 0;
  
  while ((nrPos - wrPos >= CHUNK_SIZE))
    {
      fs->write((char*)tmpSignal,CHUNK_SIZE*sizeof(storageType));
      wrPos += CHUNK_SIZE;
      if(progress) 
	cerr<<setw(4)<<setfill(' ')<<setprecision(2)<<(int)(100*((double)wrPos/(double)nrPos))<<" % complete.\r";
    }
  delete[] tmpSignal;
  for (;wrPos<=nrPos;wrPos++)
    {
      if(progress && (wrPos % 1000) == 0) 
	cerr<<setw(4)<<setfill(' ')<<setprecision(2)<<(int)(100*((double)wrPos/(double)nrPos))<<" % complete.\r";
      fs->write((char*)&tmpShort,sizeof(storageType));
    }
  // just to make sure. 
  fs->flush();
  if(progress)
    cerr<<" 100 % complete."<<endl;
}



int main(int argc, char* argv[]) 
{
	// parameters that you can set. 
	string inf     = "";
	string outf    = "";

	/*
	 *  Read and check input parameters
	 */

	string errorLine =  "usage " + 
	  string(argv[0]) + 
	  " [parameters] \n" +
	  "\t-i      <input binary file> \n" +
	  "\t-o      <output binary file> \n";
	
	bool fail = false;
	string failmessage = "";
	
	for (int i = 1; i < argc; i++)
	  {
	    if(strcmp(argv[i],"-i") == 0)
	      inf.assign(argv[++i]);
	    else if(strcmp(argv[i],"-o") == 0)
	      outf.assign(argv[++i]);
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
	
	if (strcmp(outf.c_str(), "") == 0)
	  {
	    failmessage.append("-o must be specified.\n");
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
	
	fstream of;
	of.open(outf.c_str(),ios::out | ios::in | ios::trunc | ios::binary);
	if (of.fail())
	{
		failmessage.append("ERROR: output binary file \"");
		failmessage.append(outf.c_str());
		failmessage.append("\" could not be created, aborting.\n");
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
	cerr<<"Input range:\t"<<minPos<<"\t"<<maxPos<<endl;
	
	// create an empty (zero-filles) binary file of the same size.
	cerr<<"Creating empty binary file."<<endl;
	createEmptyBinary(&of,minPos,maxPos,true);
	getMinMaxF( (ifstream*) (&of), &minPos, &maxPos);
	cerr<<"Ouput range:\t"<<minPos<<"\t"<<maxPos<<endl;
	
	// start reading the input and randomly write the non-zero positions in the output file.
	

	int chunk = CHUNK_MULTI*MAXNR;
	int read = 0;
	int range = maxPos-minPos +1;
	int pos;
	bool written = false;

	// if necessary, reduce to file size.
	chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);
	
	storageType* inpData = new storageType[chunk];
	storageType* tmpData = new storageType[1];
	cerr << "\t\t" << setprecision(3) << setw(4) << 0.0 << "  \t% complete." << endl;
	
	while ((maxPos - minPos + 1) - read >= chunk && (chunk > 0))
	  {
	    // Read data from input file
	    iff.seekg(sizeof(storageType) * (read + (minPos-minPos)) + sizeof(int), ios::beg);
	    iff.read((char*)inpData, sizeof(storageType) * chunk);
	    
	    //Loop all positions
	    for (int i = 0; i < chunk; i++)
	      {
		if(inpData[i] != 0) // only pertubate non-zeroes.
		  {
		    written = false;
		    while (!written)
		      {
			// where to write this piece of data. 
			pos = minPos + (int)((double)range * rand()/(RAND_MAX + 1.0));
			// make sure that this postion is nut used before. 
			of.seekg(sizeof(storageType) * pos + sizeof(int), ios::beg);
			of.read((char*)tmpData, sizeof(storageType));
			
			if(of.eof()) // moved to first/last postion. still ok. 
			  of.clear();
			if(tmpData[0] == 0) // ok to write. 
			  {
			    of.seekg(sizeof(storageType) * pos + sizeof(int), ios::beg);
			    if(of.eof()) // moved to first/last postion. still ok. 
			      of.clear();
			    of.write((char*)(&inpData[i]), sizeof(storageType));
			    written = true;
			  }
		      }
		  }
	      }
	    
	    read += chunk;
	    
	    // Check if we can read in a full chunk of data or simply the rest of the file. 
	    if (((maxPos - minPos + 1 - read) < chunk))
	      chunk = (maxPos - minPos + 1) - read;
	    
	    cerr << "\t\t" << setprecision(3) << setw(4) << 100.0 * (double)read / (double)(maxPos - minPos + 1) << "  \t% complete.\r";
	  }
	cerr << "\t\t" << setprecision(3) << setw(4) << 100.0 << "  \t% complete." << endl;
	
	delete[] inpData;
	delete[] tmpData;

	of.close();
	iff.close();
}
