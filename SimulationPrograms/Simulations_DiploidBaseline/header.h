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
	int n; // Deme size		
	int tfinal; // # of generations 
	int tmes; // # of generations after which dominance evolution is triggered
	int nmes; // # of generations after which dominance evolution is triggered
		
	// Initial conditions
	double a0;
	double x0;
	double z0;
	
	// Sex-specific functions
	double Opf;
	double wof;
	
	double Opm;
	double wom;
	
	// Genetic parameters
	double ux; // Sex allocation mutation rate
	double uz; // SA trait mutation rate
	double sigma;	// Size of mutation steps (StDev of the Gaussian)
	
	// Repeats
	int n_it; // # of repeats.
	int crit;
};

// Structure describing an haplotype, which encodes a sex allocation strategy 'x', a pollen export strategy 'e' and has a promoter affinity of 'a'.
struct hap
{
	double ax;
	double x; // Sex allocation allele carried
	
	double az;
	double z;
};

// Prototypes of functions


void openfileP(); // Open the parameters file

bool readpar(parameters &parr); // Read parameter values

double Vof(hap h1, hap h2, parameters parv);
double Vom(hap h1, hap h2, parameters parv);
double alloc(hap h1, hap h2);

void recursion(parameters parv, int it);

void cntl_c_handler(int bidon);

// Distributions

double gammln(const double xx);
double poisdev(const double xm);
double binldev(const double pp, const int n);
double gaussdev();

#endif
