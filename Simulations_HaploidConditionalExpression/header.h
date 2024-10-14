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
	double a0; // Initial sex allocation-dependent trait component
	double b0; // Initial sex allocation independent trait component
	
	// Cost of plasticity
	double c;
	
	// Sex-specific functions
	double Qf; // Female optimum
	double wf; // Strength of selection through female function
	
	double Qm; // Male optimum
	double wm; // Strength of selection through male function
	
	// Genetic parameters
	double ux; // Sex allocation mutation rate
	double uab; // Mutation rate for components of z		
	double sigma;	// Size of mutation steps (StDev of the Gaussian)
	
	// Repeats
	int n_it; // # of repeats.
};

// Structure describing an haplotype, which encodes a sex allocation strategy 'x', a pollen export strategy 'e' and has a promoter affinity of 'a'.
struct ind
{
	double x; // Sex allocation allele carried
	
	// Sex allocation-dependent and independent conflict trait values
	double a;
	double b;
};

// Prototypes of functions


void openfileP(); // Open the parameters file

bool readpar(parameters &parr); // Read parameter values

double RN(ind Y); // Function that gives the conflict trait value expressed by individual Y given its sex allocation and values of a and b.
double Vof(double z, parameters parv); // Fecundity effect of trait z on female function
double Vom(double z, parameters parv); // Fecundity effect of trait z on male function

void recursion(parameters parv, int it); // Main recursion function.

void cntl_c_handler(int bidon);

// Distributions

double gammln(const double xx);
double poisdev(const double xm);
double binldev(const double pp, const int n);
double gaussdev();

#endif
