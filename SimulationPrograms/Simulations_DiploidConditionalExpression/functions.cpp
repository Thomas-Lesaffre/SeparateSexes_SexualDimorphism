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
	
	double A = (h1.a*h1.pa + h2.a*h2.pa)/(h1.pa + h2.pa);
	double z = A*((h1.x*h1.px + h2.x*h2.px)/(h1.px + h2.px)) + (h1.b*h1.pb + h2.b*h2.pb)/(h1.pb + h2.pb);
	double w = exp( - parv.wof*pow(z - parv.Opf, 2.0) )*exp( -parv.c*pow(A, 2.0) );
	
	return w;
}

double Vom(hap h1, hap h2, parameters parv){	
	double A = (h1.a*h1.pa + h2.a*h2.pa)/(h1.pa + h2.pa);
	double z = A*((h1.x*h1.px + h2.x*h2.px)/(h1.px + h2.px)) + (h1.b*h1.pb + h2.b*h2.pb)/(h1.pb + h2.pb);
	double w = exp( - parv.wom*pow(z - parv.Opm, 2.0) )*exp( -parv.c*pow(A, 2.0) );
	
	return w;
}

double alloc(hap h1, hap h2){	
	double x = (h1.x*h1.px + h2.x*h2.px)/(h1.px + h2.px);
	return x;
}
