#ifndef DEPRESSION_H
#define DEPRESSION_H

#include <vector>
#include <iostream>
#include <future>
#include <thread>
#include <iterator>
#include <random>

using namespace std;

// global variables

#define filePar "par" 

struct parameters
{
	// General parameters
	int n; // Population size		
	int tfinal; // # of generations 
	int tmes; // # of generations after which dominance evolution is triggered
	int nmes; // # of generations after which dominance evolution is triggered
		
	// Initial conditions
	double x0; // Initial sex allocation
	double z0; // Initial conflict trait
	double r0; // Initial recombination rate
	
	// Sex-specific functions
	double Opf; // Female optimum
	double wof; // Strength of selection through female function
	
	double Opm; // Male optimum
	double wom; // Strength of selection through male function
	
	// Genetic parameters
	double ux; // Sex allocation mutation rate
	double uz; // Conflict trait mutation rate
	double ur; // Recombination mutation rate		
	double sigma;	// Size of mutation steps for x and z (StDev of the Gaussian)
	double sigmarec; // Size of mutation steps for the recombination rate (StDev of the Gaussian)
	
	// Repeats
	int n_it; // # of repeats.
};

// Structure describing an haplotype, which encodes a sex allocation strategy 'x', a pollen export strategy 'e' and has a promoter affinity of 'a'.
struct ind
{
	double x; // Sex allocation allele carried
	double z; // Conflict trait allele
	double r; // Allele at the recombination rate modifier
};

// Prototypes of functions


void openfileP(); // Open the parameters file

bool readpar(parameters &parr); // Read parameter values

double Vof(double z, parameters parv); // Function capturing the effect of z on female fecundity
double Vom(double z, parameters parv); // Function capturing the effect of z on male fecundity

void recursion(parameters parv, int it); // Main recursion function

void cntl_c_handler(int bidon);

// Distributions

double gammln(const double xx);
double poisdev(const double xm);
double binldev(const double pp, const int n);
double gaussdev();

#endif
