// Functions to open input and output files,
// read parameter values and write them in output file

#include "header.h"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

extern FILE * fileP;

//Opens input files:

// Parameters file
void openfileP()
{
	fileP = fopen(filePar,"r");
	if (!fileP)
		cout << "The file " << filePar << " doesn't exist!"  << endl;
}

// Reads parameter values from input file,
// returns 1 if it reaches the end of input file, else returns 0:

bool readpar(parameters &parr)
{
	int z;
	bool term;
	do {z = fgetc(fileP);} while (!((z == '*') || (z == EOF)));
		// Lines with parameter sets must begin with *
	if (z == EOF)
	{
		cout << "\nEnd of input file\n";
		term = true;
	}
	else
	{	
		// Parameters
	
		// General parameters
		fscanf(fileP," %d",&parr.n);				
		fscanf(fileP," %d",&parr.tfinal);
		fscanf(fileP," %d",&parr.tmes);
		fscanf(fileP," %d",&parr.nmes);
				
		// Initial conditions
		fscanf(fileP," %lf",&parr.a0);					
		fscanf(fileP," %lf",&parr.x0);			
		fscanf(fileP," %lf",&parr.z0);

		// Sex-specific functionsB	
		fscanf(fileP," %lf",&parr.Opf);			
		fscanf(fileP," %lf",&parr.wof);			
		
		fscanf(fileP," %lf",&parr.Opm);			
		fscanf(fileP," %lf",&parr.wom);									

		// Genetic parameters
		fscanf(fileP," %lf",&parr.ux);	
		fscanf(fileP," %lf",&parr.uz);			
		fscanf(fileP," %lf",&parr.sigma);			
					
		// Repeats
		fscanf(fileP," %d",&parr.n_it);	
		fscanf(fileP," %d",&parr.crit);			
								
																					
        term = false;
	}
	
	return term;
}

