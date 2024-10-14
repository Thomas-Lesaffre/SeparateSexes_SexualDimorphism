#include "header.h"
#include "mt.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <future>
#include <chrono>
#include <thread>

using namespace std;

extern MTRand eng;

// Effect of z on female fecundity
double Vof(double z, parameters parv){	
	double w = exp( - parv.wof*pow(z - parv.Opf, 2.0) ); // Gaussian function
	return w;
}

// Effect of z on male fecundity
double Vom(double z, parameters parv){	
	double w = exp( - parv.wom*pow(z - parv.Opm, 2.0) ); // Gaussian function
	return w;
}
