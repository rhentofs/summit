#ifndef INFO_H
#define INFO_H

#include <iostream>

using namespace std;

class Info
{
public:
	Info(double first, double second, double third, double fourth=0.0)
		{
			info = new double[4];
			info[0] = first;
			info[1] = second;
			info[2] = third;
			info[3] = fourth;
		}
	~Info()
		{
			delete[] info;
		}
	double get(int pos)
		{
			return(info[pos]);
		}
	void print()
		{
			cerr << "***Info***" << endl;
			cerr << "\t" << info[0] << endl;
			cerr << "\t" << info[1] << endl;
			cerr << "\t" << info[2] << endl;
			cerr << "\t" << info[3] << endl;
			cerr << "**********" << endl;
		}

private:
	double *info;
};

#endif
