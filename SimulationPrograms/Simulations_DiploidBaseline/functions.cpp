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

double Vof(hap h1, hap h2, parameters parv){	
	double z = (h1.z*h1.az + h2.z*h2.az)/(h1.az + h2.az);
	double w = exp( - parv.wof*pow(z - parv.Opf, 2.0) );
	
	return w;
}

double Vom(hap h1, hap h2, parameters parv){	
	double z = (h1.z*h1.az + h2.z*h2.az)/(h1.az + h2.az);
	double w = exp( - parv.wom*pow(z - parv.Opm, 2.0) );
	
	return w;
}

double alloc(hap h1, hap h2){	
	double x = (h1.x*h1.ax + h2.x*h2.ax)/(h1.ax + h2.ax);
	return x;
}
