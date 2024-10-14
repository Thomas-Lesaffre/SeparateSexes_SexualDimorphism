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
	sstream_x << "res_diplo_rec_dom_x_Of" << parv.Opf << "_Om" << parv.Opm << "_wf" << parv.wof << "_wm" << parv.wom << "_ux" << parv.ux << "_uz" << parv.uz << "-" << it;
	sstream_x >> filename_x;
	ofstream foutx(filename_x);
	
	char filename_ax[256];
	stringstream sstream_ax;
	sstream_ax << "res_diplo_rec_dom_ax_Of" << parv.Opf << "_Om" << parv.Opm << "_wf" << parv.wof << "_wm" << parv.wom << "_ux" << parv.ux << "_uz" << parv.uz << "-" << it;
	sstream_ax >> filename_ax;
	ofstream foutax(filename_ax);
	
	char filename_z[256];
	stringstream sstream_z;
	sstream_z << "res_diplo_rec_dom_z_Of" << parv.Opf << "_Om" << parv.Opm << "_wf" << parv.wof << "_wm" << parv.wom << "_ux" << parv.ux << "_uz" << parv.uz << "-" << it;
	sstream_z >> filename_z;
	ofstream foutz(filename_z);

	char filename_az[256];
	stringstream sstream_az;
	sstream_az << "res_diplo_rec_dom_az_Of" << parv.Opf << "_Om" << parv.Opm << "_wf" << parv.wof << "_wm" << parv.wom << "_ux" << parv.ux << "_uz" << parv.uz << "-" << it;
	sstream_az >> filename_az;
	ofstream foutaz(filename_az);
   
	hap tmp; // Dummy haplotype

	int i, t, k, l, d, p, nb, val;
	double r, s, tm, tf, rec;
	
	int tn = 2*parv.n;
	
	vector<hap> pop(tn); // Vector containing the population
	vector<hap> par(tn); // Vector containing the parental population
		
	vector<double> vf(parv.n), vm(parv.n); // Vectors of male and female fecundities

	// Initialising the population
	
	tmp.x = parv.x0;
	tmp.ax = parv.a0;
	
	tmp.z = parv.z0;
	tmp.az = parv.a0;
	
	for(i=0; i< tn; i++) // For each haplotype
	{
		pop[i] = tmp; // Push itinto the population and parental population vectors
		par[i] = tmp;
	}
	
	for(t=0; t<parv.tfinal; t++) // For each time step,
	{		
		// cout << t << endl;
		
	
		// Compute fecundities
			
		for(i=0; i<parv.n; i++)
		{
			nb=2*i;		
			
			par[nb] = pop[nb];
			par[nb+1] = pop[nb+1];
							
			if(nb == 0)
			{
				vf[i] = Vof(pop[nb], pop[nb+1], parv)*alloc(pop[nb], pop[nb+1]);	
				vm[i] = Vom(pop[nb], pop[nb+1], parv)*(1.0 - alloc(pop[nb], pop[nb+1]));									
			}
			else
			{
				vf[i] = vf[i-1] + Vof(pop[nb], pop[nb+1], parv)*alloc(pop[nb], pop[nb+1]);	
				vm[i] = vm[i-1] + Vom(pop[nb], pop[nb+1], parv)*(1.0 - alloc(pop[nb], pop[nb+1]));					
			}	
		}
		
		tm = vm.back();
		tf = vf.back();

		// cout << "a" << endl;
				
		// Create the next generation
		rec=0;
		
		for(i=0; i<parv.n; i++)
		{
			// cout << "i=" << i << endl;		
			nb=2*i;
			
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
			
			
			//// Maternally inherited haplotype
				
			if(eng.rand() < 0.5){
				pop[nb].ax = par[2*k].ax;				
				pop[nb].x = par[2*k].x;
				
				pop[nb].az = par[2*k].az;										
				pop[nb].z = par[2*k].z;
			}else{
				pop[nb].ax = par[2*k + 1].ax;										
				pop[nb].x = par[2*k + 1].x;
				
				pop[nb].az = par[2*k + 1].az;					
				pop[nb].z = par[2*k + 1].z;	
			}
			
			if( eng.rand() < parv.ux ){
				pop[nb].x += parv.sigma*gaussdev();				
				if(pop[nb].x < 0){ pop[nb].x = 0; }
				if(pop[nb].x > 1){ pop[nb].x = 1; }		
				
				pop[nb].ax += parv.sigma*gaussdev();
				if(pop[nb].ax < 0.000001){ pop[nb].ax = 0.000001; }
			}
	
			if( eng.rand() < parv.uz ){
				pop[nb].z += parv.sigma*gaussdev();		

				pop[nb].az += parv.sigma*gaussdev();
				if(pop[nb].az < 0.000001){ pop[nb].az = 0.000001; }					
			}
			
			//// Paternally inherited haplotype	
			
			if(eng.rand() < 0.5){
				pop[nb+1].ax = par[2*l].ax;				
				pop[nb+1].x = par[2*l].x;
				
				pop[nb+1].az = par[2*l].az;										
				pop[nb+1].z = par[2*l].z;
			}else{
				pop[nb+1].ax = par[2*l + 1].ax;										
				pop[nb+1].x = par[2*l + 1].x;
				
				pop[nb+1].az = par[2*l + 1].az;					
				pop[nb+1].z = par[2*l + 1].z;	
			}

			
			if( eng.rand() < parv.ux ){
				pop[nb+1].x += parv.sigma*gaussdev();				
				if(pop[nb+1].x < 0){ pop[nb+1].x = 0; }
				if(pop[nb+1].x > 1){ pop[nb+1].x = 1; }		
				
				pop[nb+1].ax += parv.sigma*gaussdev();
				if(pop[nb+1].ax < 0.000001){ pop[nb+1].ax = 0.000001; }
			}
	
			if( eng.rand() < parv.uz ){
				pop[nb+1].z += parv.sigma*gaussdev();		

				pop[nb+1].az += parv.sigma*gaussdev();
				if(pop[nb+1].az < 0.000001){ pop[nb+1].az = 0.000001; }					
			}			
		}
		
		// Measurements
		
		if(t % parv.tmes == 0) // Every tmes generations
		{										
			for(i=0; i<parv.nmes; i++) // Record the phenotype of nmes individual in output files.
			{		
				nb = 2*i;		
				if(i == 0)
				{
					foutx << pop[nb].x << " " << pop[nb+1].x;
					foutax << pop[nb].ax << " " << pop[nb+1].ax;
					
					foutz << pop[nb].z << " " << pop[nb+1].z;
					foutaz << pop[nb].az << " " << pop[nb+1].az;					
				}
				else
				{
					foutx << " " << pop[nb].x << " " << pop[nb+1].x;
					foutax << " " << pop[nb].ax << " " << pop[nb+1].ax;
					
					foutz << " " << pop[nb].z << " " << pop[nb+1].z;
					foutaz << " " << pop[nb].az << " " << pop[nb+1].az;					
				}
			}

			foutx << endl;
			foutax << endl;
						
			foutz << endl;	
			foutaz << endl;

			
			if(parv.crit > 0.50){ // If crit > 0.50, then we have a stopping criterion acting: if no individual expresses a sex allocation between 0.05 and 0.95, the program automatically stops.
				val=0;
				for(i=0; i<parv.n; i++){
					nb=2*i;					
					if( alloc(pop[nb], pop[nb+1]) < 0.95 && alloc(pop[nb], pop[nb+1]) > 0.05 ){
						val++;
						break;
					}
				}
				
				if(val < 0.5){
					break;
				}
			}
		} // End of 'if' for writing
					
	} // End of loop over time	
} // End of function
