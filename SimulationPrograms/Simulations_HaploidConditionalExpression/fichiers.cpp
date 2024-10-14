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
		fscanf(fileP," %lf",&parr.x0);			
		fscanf(fileP," %lf",&parr.a0);
		fscanf(fileP," %lf",&parr.b0);
		
		// Cost of plasticity
		fscanf(fileP," %lf",&parr.c);
				
		// Sex-specific functions	
		fscanf(fileP," %lf",&parr.Qf);			
		fscanf(fileP," %lf",&parr.wf);			
		
		fscanf(fileP," %lf",&parr.Qm);			
		fscanf(fileP," %lf",&parr.wm);									

		// Genetic parameters
		fscanf(fileP," %lf",&parr.ux);	
		fscanf(fileP," %lf",&parr.uab);									
		fscanf(fileP," %lf",&parr.sigma);			

		// Repeats
		fscanf(fileP," %d",&parr.n_it);	
								
																					
        term = false;
	}
	
	return term;
}

