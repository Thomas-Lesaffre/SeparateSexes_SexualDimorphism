#include "header.h"
#include "mt.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <future>
#include <thread>
#include <algorithm>
#include <iterator>

using namespace std;

MTRand eng;

// Pointers on input and output files:

FILE * fileP;

int main()
{ 	
    cout << "Program initialization" << "\n";

	parameters par;
	
	openfileP();
	
	bool end { false };
	int k;
	
	cout << "coucou" << endl;
	
	end = readpar(par);
	
	
	for(k=0; k<par.n_it; k++)
	{
    	recursion(par, k+1); // Run the recursion as many times as indicated by the parameter file.
	}	
		
	return 0;
}


