#include "zeb_utils.h"

int bin2bed(ifstream *ifi, ofstream *ofi, string chrom)
{
	int minPos, maxPos;
	getMinMaxF(ifi, &minPos, &maxPos);

	int read = 0;
	int chunk = CHUNK_MULTI*MAXNR;

	// if necessary, reduce to file size.
	chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);

	storageType *tmp = new storageType[chunk];

	int maxVal = 0;
	
	while ((maxPos - minPos + 1 - read) >= chunk && chunk > 0)
	{
		ifi->seekg(sizeof(storageType) * (read) + sizeof(int), ios::beg);
		ifi->read((char*)tmp, sizeof(storageType) * chunk);

		for (int i = 0; i < chunk; i++)
		{
			if ((int)tmp[i] > 0)
			{
				if ((int)tmp[i] > maxVal)
					maxVal = (int)tmp[i];
			}
		}
		
		read += chunk;
		
		// Check if we can read in a full chunk of data or simply the rest of the file. 
		if (((maxPos - minPos + 1 - read) < chunk))
			chunk = (maxPos - minPos + 1) - read;
	}
	
	if (ifi->eof())
		ifi->clear();

	read = 0;
	chunk = CHUNK_MULTI*MAXNR;

	// if necessary, reduce to file size.
	chunk = ((maxPos - minPos + 1) < chunk ? (maxPos - minPos + 1) : chunk);

	while ((maxPos - minPos + 1 - read) >= chunk && chunk > 0)
	{
		ifi->seekg(sizeof(storageType) * (read) + sizeof(int), ios::beg);
		ifi->read((char*)tmp, sizeof(storageType) * chunk);

		for (int i = 0; i < chunk; i++)
		{
			if ((int)tmp[i] > 0)
			{
				(*ofi) << chrom << "\t" << minPos + read + i - 1 << "\t" << minPos + read + i << "\tpos\t";
				(*ofi) << np_round(1.0 + (999.0 * ((double)tmp[i] / (double)maxVal))) << endl;
			}
		}
		
		read += chunk;
		
		// Check if we can read in a full chunk of data or simply the rest of the file. 
		if (((maxPos - minPos + 1 - read) < chunk))
			chunk = (maxPos - minPos + 1) - read;
	}
	
	if (ifi->eof())
		ifi->clear();

	delete[] tmp;
	
	return(1);
}

int vlevel = 3;
int main(int argc, char* argv[]) 
{
	string inf   = "";
	string outf  = "";
	string name  = "";
	string chrom = "";

	string errorLine =  "usage " + 
		string(argv[0]) + 
		" [parameters] \n" +
		"\t-i\t<binary infile>\n" +
		"\t-o\t<bed outfile>\n" +
		"\t-c\t<chromosome name>\n" +
		"\t-n\t<track name>";

	vlevel = 3;

	bool fail = false;
	string failmessage = "";
	
	for (int i = 1; i < argc; i++)
	{
	    if(strcmp(argv[i],"-i") == 0)
			inf.assign(argv[++i]);
	    else if(strcmp(argv[i],"-o") == 0)
			outf.assign(argv[++i]);
		else if(strcmp(argv[i],"-n") == 0)
			name.assign(argv[++i]);
		else if(strcmp(argv[i],"-c") == 0)
			chrom.assign(argv[++i]);
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

	if (strcmp(name.c_str(), "") == 0)
	{
		failmessage.append("-n must be specified.\n");
		fail = true;
	}

	if (strcmp(chrom.c_str(), "") == 0)
	{
		failmessage.append("-c must be specified.\n");
		fail = true;
	}
	
	if (fail)
	{
		cerr << endl << failmessage.c_str() << endl << errorLine << endl;
		return(-1);
	}

	ifstream ifi;
	ofstream ofi;

	ifi.open(inf.c_str(),ios::binary);
	ofi.open(outf.c_str(),ios::trunc);

	if (ifi.fail())
	{
		failmessage.append("ERROR: Input file \"");
		failmessage.append(inf.c_str());
		failmessage.append("\" could not be opened, aborting.\n");
		fail = true;
	}
	
	if (ofi.fail())
	{
		failmessage.append("ERROR: Output file \"");
		failmessage.append(outf.c_str());
		failmessage.append("\" could not be created, aborting.\n");
		fail = true;
	}

	if (fail)
	{
		cerr << endl << failmessage.c_str() << endl << errorLine << endl;
		return(-1);
	}

	ofi << "track name=\"" << name << "\" description=\"\" visibility=2" << endl;
	
	bin2bed(&ifi,&ofi,chrom);

	ifi.close();
	ofi.close();
	
	return(1);
}


	
		

	
