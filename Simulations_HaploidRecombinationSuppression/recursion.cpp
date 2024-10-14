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
	sstream_x << "res_rec_x_Of" << parv.Opf << "_Om" << parv.Opm << "_wf" << parv.wof << "_wm" << parv.wom << "_ux" << parv.ux << "_uz" << parv.uz << "_ur" << parv.ur << "-" << it;
	sstream_x >> filename_x; 
	ofstream foutx(filename_x); // Sex allocation values
	
	char filename_z[256];
	stringstream sstream_z;
	sstream_z << "res_rec_z_Of" << parv.Opf << "_Om" << parv.Opm << "_wf" << parv.wof << "_wm" << parv.wom << "_ux" << parv.ux << "_uz" << parv.uz << "_ur" << parv.ur<< "-" << it;
	sstream_z >> filename_z;
	ofstream foutz(filename_z); // Conflict trait values

	char filename_r[256];
	stringstream sstream_r;
	sstream_r << "res_rec_r_Of" << parv.Opf << "_Om" << parv.Opm << "_wf" << parv.wof << "_wm" << parv.wom << "_ux" << parv.ux << "_uz" << parv.uz << "_ur" << parv.ur << "-" << it;
	sstream_r >> filename_r;
	ofstream foutr(filename_r); // Recombination rates
   
	ind tmp; // Dummy haplotype
	
	// Declaring variables
	
	int i, t, k, l, d, p;
	double r, s, tm, tf, rec;
	
	vector<ind> pop(parv.n); // Vector containing the population
	vector<ind> par(parv.n); // Vector containing the parental population
		
	vector<double> vf(parv.n), vm(parv.n); // Vectors of male and female fecundities

	// Initialising the population
	
	tmp.x = parv.x0;
	tmp.z = parv.z0;
	tmp.r = parv.r0;
	
	for(i=0; i< parv.n; i++) // For each chromosome (one per individual)
	{
		pop[i] = tmp; // Push the chromosome into the population and parental population vectors
		par[i] = tmp;
	}
	
	for(t=0; t<parv.tfinal; t++) // For each time step,
	{		
		// cout << t << endl;
		
		// Compute female and male fecundities and store them in vectors of cumulated fecundities. 
			
		for(i=0; i<parv.n; i++)
		{
			par[i] = pop[i];
				
			if(i == 0)
			{
				vf[i] = pop[i].x*Vof(pop[i].z, parv);	
				vm[i] = (1.0 - pop[i].x)*Vom(pop[i].z, parv);									
			}
			else
			{
				vf[i] = vf[i-1] + pop[i].x*Vof(pop[i].z, parv);
				vm[i] = vm[i-1] + (1.0 - pop[i].x)*Vom(pop[i].z, parv);				
			}	
		}
		
		tm = vm.back(); // Max male fecundity
		tf = vf.back(); // Male female fecundity
		
		// Create the next generation
		rec=0; // Set number of recombinant offspring to zero.
		
		for(i=0; i<parv.n; i++) // For each recruited offspring
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
			
			// The recombination modifier is freely recombining with the other loci.			
			
			if( eng.rand() < 0.5 ){
				pop[i].r = par[k].r;
			}else{
				pop[i].r = par[l].r;
			}

			// Recombination occurs between the loci encoding x and z with a probability given by the diploid genotype of the zygote.
			
			if( eng.rand() < (par[k].r + par[l].r)/2.0 ){					
				rec += 1; // If recombination occurs, count one additional recombinant.		
				if(eng.rand() < 0.5){
					pop[i].x = par[k].x;
					pop[i].z = par[l].z;
				}else{
					pop[i].x = par[l].x;
					pop[i].z = par[k].z;
				}
			}else{				
				if(eng.rand() < 0.5){
					pop[i].x = par[k].x;
					pop[i].z = par[k].z;
				}else{
					pop[i].x = par[l].x;
					pop[i].z = par[l].z;	
				}
			}
			
			// Mutation occurs for each trait depending on their respective mutation rates
			
			if( eng.rand() < parv.ux ){
				pop[i].x += parv.sigma*gaussdev();
		
				if(pop[i].x < 0){ pop[i].x = 0; }
				if(pop[i].x > 1){ pop[i].x = 1; }		
			}
			
			if( eng.rand() < parv.uz ){
				pop[i].z += parv.sigma*gaussdev();			
			}
			
			if( eng.rand() < parv.ur ){
				
				pop[i].r += parv.sigmarec*gaussdev();
														
				if(pop[i].r < 0){ pop[i].r = 0; }
				if(pop[i].r > 0.5){ pop[i].r = 0.5; }
			}	
		}

		// Measurements
		
		if(t % parv.tmes == 0) // Every ten generations
		{	
			rec /= parv.n; // Get the proportion of recombinants in the new generation (proxy for recombination rate).
									
			for(i=0; i<parv.nmes; i++) // Record the sex allocation and conflict trait values of 'nmes' haploid individuals.
			{				
				if(i == 0)
				{
					foutx << pop[i].x;
					foutz << pop[i].z;
				}
				else
				{
					foutx << " " << pop[i].x;
					foutz << " " << pop[i].z;
				}
			}

			foutx << endl;
			foutz << endl;	
			foutr << rec << endl;	
			
			// If the recombination rate is lower than the chosen value, stop the simulation. This is useful only when running trials for recombination suppression. 
			// The threshold is set to -1 such that this can never be satisfied.
			if(rec < -1){ 
				break;
			}
		} // End of 'if' for writing
					
	} // End of loop over time	
} // End of function
