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

// Function that gives the conflict trait value expressed by individual Y given its sex allocation and values of a and b.
double RN(ind Y){	
	double z = Y.a*Y.x + Y.b;
	return z;
}

// Fecundity effect of trait z on female function
double Vof(double z, parameters parv){	
	double w = exp( - parv.wf*pow(z - parv.Qf, 2.0) );
	return w;
}

// Fecundity effect of trait z on male function
double Vom(double z, parameters parv){	
	double w = exp( - parv.wm*pow(z - parv.Qm, 2.0) );
	return w;
}
