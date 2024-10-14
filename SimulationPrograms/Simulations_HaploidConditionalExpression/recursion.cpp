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
#include <random>


using namespace std;

extern MTRand eng;

// Recursion function

void recursion(parameters parv, int it)
{
	
	// Creating output files for phenotypes, sex allocation genotypes and promoter affinities
	
	char filename_x[256];
	stringstream sstream_x;
	sstream_x << "res_rn_x_Qf" << parv.Qf << "_Qm" << parv.Qm << "_wf" << parv.wf << "_wm" << parv.wm << "_ux" << parv.ux << "_ur" << parv.ur << "-" << it;
	sstream_x >> filename_x;
	ofstream foutx(filename_x);
	
	char filename_a[256];
	stringstream sstream_a;
	sstream_a << "res_rn_a_Qf" << parv.Qf << "_Qm" << parv.Qm << "_wf" << parv.wf << "_wm" << parv.wm << "_ux" << parv.ux << "_ur" << parv.ur << "-" << it;
	sstream_a >> filename_a;
	ofstream fouta(filename_a);
	
	char filename_b[256];
	stringstream sstream_b;
	sstream_b << "res_rn_b_Qf" << parv.Qf << "_Qm" << parv.Qm << "_wf" << parv.wf << "_wm" << parv.wm << "_ux" << parv.ux << "_ur" << parv.ur << "-" << it;
	sstream_b >> filename_b;
	ofstream foutb(filename_b);
   
   	// Declare variables
   	
	ind tmp; // Dummy haplotype

	int i, t, k, l, d, p, stop;
	double r, s, tm, tf;
	
	vector<ind> pop(parv.n); // Vector containing the population
	vector<ind> par(parv.n); // Vector containing the parental population
		
	vector<double> vf(parv.n), vm(parv.n); // Vectors of male and female fecundities

	// Initialising the population
	
	tmp.x = parv.x0;
	tmp.a = parv.a0;
	tmp.b = parv.b0;
	
	for(i=0; i< parv.n; i++) // For each chromosome (2 per individual)
	{
		pop[i] = tmp; // Push the chromosome into the population and parental population vectors
		par[i] = tmp;
	}
	
	for(t=0; t<parv.tfinal; t++) // For each time step,
	{		
		// cout << t << endl;
		
		// Compute fecundities and store them in cumulated vectors
			
		for(i=0; i<parv.n; i++)
		{
			par[i] = pop[i];
				
			if(i == 0)
			{
				vf[i] = pop[i].x*Vof(RN(pop[i]), parv)*exp( -parv.c*pow(pop[i].a, 2.0) );	
				vm[i] = (1.0 - pop[i].x)*Vom(RN(pop[i]), parv)*exp( -parv.c*pow(pop[i].a, 2.0) );									
			}
			else
			{
				vf[i] = vf[i-1] + pop[i].x*Vof(RN(pop[i]), parv)*exp( -parv.c*pow(pop[i].a, 2.0) );
				vm[i] = vm[i-1] + (1.0 - pop[i].x)*Vom(RN(pop[i]), parv)*exp( -parv.c*pow(pop[i].a, 2.0) );				
			}	
		}
		
		tm = vm.back(); // Max male fecundity
		tf = vf.back(); // Max female fecundity
		
		// Create the next generation
		
		for(i=0; i<parv.n; i++)
		{
			// Find a mother
			k=-1;
			s=eng.rand(tf);
			do{
				k++;
			}while(vf[k] < s);
			
			// Find a father
			l=-1;
			s=eng.rand(tm);
			do{
				l++;
			}while(vm[l] < s);
			
			// Loci are all freely recombining.		
			
			if( eng.rand() < 0.5 ){
				pop[i].x = par[k].x;
				
			}else{
				pop[i].x = par[l].x;
			}

			if( eng.rand() < 0.5 ){
				pop[i].a = par[k].a;
				
			}else{
				pop[i].a = par[l].a;
			}
			
			if( eng.rand() < 0.5 ){
				pop[i].b = par[k].b;
				
			}else{
				pop[i].b = par[l].b;
			}
									
			// Mutation
			
			if( eng.rand() < parv.ux ){
				pop[i].x += parv.sigma*gaussdev();
		
				if(pop[i].x < 0){ pop[i].x = 0; }
				if(pop[i].x > 1){ pop[i].x = 1; }		
			}

			if( eng.rand() < parv.uab ){
				pop[i].a += parv.sigma*gaussdev();
				pop[i].b += parv.sigma*gaussdev();
			}	
		}

		// cout << "r = " << rec << endl;		
		// Measurements
		
		if(t % parv.tmes == 0) // Every 'tmes' generations
		{										
			for(i=0; i<parv.nmes; i++) // Record the values encoded by the haplotypes of 'nmes' individuals in output files.
			{				
				if(i == 0)
				{
					foutx << pop[i].x;
					fouta << pop[i].a;
					foutb << pop[i].b;					
				}
				else
				{
					foutx << " " << pop[i].x;
					fouta << " " << pop[i].a;
					foutb << " " << pop[i].b;	
				}
			}

			foutx << endl;
			fouta << endl;	
			foutb << endl;	
			
			stop=0; // This can be set to 1 for the simulation to automatically stop when no individual expresses a sex allocation value between 0.1 and 0.9. T
			// This option is handy for trials involving many replicates
			for(i=0; i<parv.n; i++){
				if( pop[i].x > 0.10 && pop[i].x < 0.90  ){
					stop=0;
					break;
				}
			}
		} // End of 'if' for writing
		
		if(stop == 1){
			break;
		}			
	} // End of loop over time	
} // End of function
