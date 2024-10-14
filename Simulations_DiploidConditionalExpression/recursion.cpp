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
	sstream_x << "res_diplo_rec_dom_x_Of" << parv.Opf << "_Om" << parv.Opm << "_wf" << parv.wof << "_wm" << parv.wom << "_ux" << parv.ux << "_uz" << parv.uz<< "-" << it;
	sstream_x >> filename_x;
	ofstream foutx(filename_x);
	
	char filename_px[256];
	stringstream sstream_px;
	sstream_px << "res_diplo_rec_dom_px_Of" << parv.Opf << "_Om" << parv.Opm << "_wf" << parv.wof << "_wm" << parv.wom << "_ux" << parv.ux << "_uz" << parv.uz << "-" << it;
	sstream_px >> filename_px;
	ofstream foutpx(filename_px);
	
	char filename_a[256];
	stringstream sstream_a;
	sstream_a << "res_diplo_rec_dom_a_Of" << parv.Opf << "_Om" << parv.Opm << "_wf" << parv.wof << "_wm" << parv.wom << "_ux" << parv.ux << "_uz" << parv.uz << "-" << it;
	sstream_a >> filename_a;
	ofstream fouta(filename_a);

	char filename_pa[256];
	stringstream sstream_pa;
	sstream_pa << "res_diplo_rec_dom_pa_Of" << parv.Opf << "_Om" << parv.Opm << "_wf" << parv.wof << "_wm" << parv.wom << "_ux" << parv.ux << "_uz" << parv.uz << "-" << it;
	sstream_pa >> filename_pa;
	ofstream foutpa(filename_pa);
	
	char filename_b[256];
	stringstream sstream_b;
	sstream_b << "res_diplo_rec_dom_b_Of" << parv.Opf << "_Om" << parv.Opm << "_wf" << parv.wof << "_wm" << parv.wom << "_ux" << parv.ux << "_uz" << parv.uz << "-" << it;
	sstream_b >> filename_b;
	ofstream foutb(filename_b);

	char filename_pb[256];
	stringstream sstream_pb;
	sstream_pb << "res_diplo_rec_dom_pb_Of" << parv.Opf << "_Om" << parv.Opm << "_wf" << parv.wof << "_wm" << parv.wom << "_ux" << parv.ux << "_uz" << parv.uz << "-" << it;
	sstream_pb >> filename_pb;
	ofstream foutpb(filename_pb);
   
	hap tmp; // Dummy haplotype

	int i, t, k, l, d, p, nb, val;
	double r, s, tm, tf, rec;
	
	int tn = 2*parv.n;
	
	vector<hap> pop(tn); // Vector containing the population
	vector<hap> par(tn); // Vector containing the parental population
		
	vector<double> vf(parv.n), vm(parv.n); // Vectors of male and female fecundities
	
	// Initialising the population
	
	tmp.x = parv.x0;
	tmp.px = parv.p0;
	
	tmp.a = parv.a0;
	tmp.pa = parv.p0;

	tmp.b = parv.b0;
	tmp.pb = parv.p0;
		
	for(i=0; i< tn; i++) // For each haplotype
	{
		pop[i] = tmp; // Push itinto the population and parental population vectors
		par[i] = tmp;
	}
	
	for(t=0; t<parv.tfinal; t++) // For each time step,
	{		
		cout << t << endl;
		
	
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
				
		// Create the next generation
				
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
			
			// All loci are freely recombining.			
			
			if(eng.rand() < 0.5){
				pop[nb].x = par[2*k].x;
				pop[nb].px = par[2*k].px;
			}else{
				pop[nb].x = par[2*k+1].x;
				pop[nb].px = par[2*k+1].px;			
			}			
				
			if(eng.rand() < 0.5){
				pop[nb].a = par[2*k].a;
				pop[nb].pa = par[2*k].pa;
			}else{
				pop[nb].a = par[2*k+1].a;
				pop[nb].pa = par[2*k+1].pa;			
			}		
			
			if(eng.rand() < 0.5){
				pop[nb].b = par[2*k].b;
				pop[nb].pb = par[2*k].pb;
			}else{
				pop[nb].b = par[2*k+1].b;
				pop[nb].pb = par[2*k+1].pb;			
			}					
								
			if( eng.rand() < parv.ux ){
				pop[nb].x += parv.sigma*gaussdev();				
				if(pop[nb].x < 0){ pop[nb].x = 0; }
				if(pop[nb].x > 1){ pop[nb].x = 1; }		
				
				pop[nb].px += parv.sigma*gaussdev();
				if(pop[nb].px < 0.000001){ pop[nb].px = 0.000001; }
			}
	
			if( eng.rand() < parv.uz ){
				pop[nb].a += parv.sigma*gaussdev();		
				pop[nb].pa += parv.sigma*gaussdev();
				if(pop[nb].pa < 0.000001){ pop[nb].pa = 0.000001; }	
				
				pop[nb].b += parv.sigma*gaussdev();		
				pop[nb].pb += parv.sigma*gaussdev();
				if(pop[nb].pb < 0.000001){ pop[nb].pb = 0.000001; }				
			}

			
			//// Paternally inherited haplotype	
			
			if(eng.rand() < 0.5){
				pop[nb+1].x = par[2*l].x;
				pop[nb+1].px = par[2*l].px;
			}else{
				pop[nb+1].x = par[2*l+1].x;
				pop[nb+1].px = par[2*l+1].px;			
			}			
				
			if(eng.rand() < 0.5){
				pop[nb+1].a = par[2*l].a;
				pop[nb+1].pa = par[2*l].pa;
			}else{
				pop[nb+1].a = par[2*l+1].a;
				pop[nb+1].pa = par[2*l+1].pa;			
			}		
			
			if(eng.rand() < 0.5){
				pop[nb+1].b = par[2*l].b;
				pop[nb+1].pb = par[2*l].pb;
			}else{
				pop[nb+1].b = par[2*l+1].b;
				pop[nb+1].pb = par[2*l+1].pb;			
			}	
								
			if( eng.rand() < parv.ux ){
				pop[nb+1].x += parv.sigma*gaussdev();				
				if(pop[nb+1].x < 0){ pop[nb+1].x = 0; }
				if(pop[nb+1].x > 1){ pop[nb+1].x = 1; }		
				
				pop[nb+1].px += parv.sigma*gaussdev();
				if(pop[nb+1].px < 0.000001){ pop[nb+1].px = 0.000001; }
			}
	
			if( eng.rand() < parv.uz ){
				pop[nb+1].a += parv.sigma*gaussdev();		
				pop[nb+1].pa += parv.sigma*gaussdev();
				if(pop[nb+1].pa < 0.000001){ pop[nb+1].pa = 0.000001; }	
				
				pop[nb+1].b += parv.sigma*gaussdev();		
				pop[nb+1].pb += parv.sigma*gaussdev();
				if(pop[nb+1].pb < 0.000001){ pop[nb+1].pb = 0.000001; }				
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
					foutpx << pop[nb].px << " " << pop[nb+1].px;
					
					fouta << pop[nb].a << " " << pop[nb+1].a;
					foutpa << pop[nb].pa << " " << pop[nb+1].pa;
					
					foutb << pop[nb].b << " " << pop[nb+1].b;
					foutpb << pop[nb].pb << " " << pop[nb+1].pb;									
				}
				else
				{
					foutx << " " << pop[nb].x << " " << pop[nb+1].x;
					foutpx << " " << pop[nb].px << " " << pop[nb+1].px;
					
					fouta << " " << pop[nb].a << " " << pop[nb+1].a;
					foutpa << " " << pop[nb].pa << " " << pop[nb+1].pa;
					
					foutb << " " << pop[nb].b << " " << pop[nb+1].b;
					foutpb << " " << pop[nb].pb << " " << pop[nb+1].pb;				
				}
			}

			foutx << endl;
			foutpx << endl;

			fouta << endl;
			foutpa << endl;
			
			foutb << endl;
			foutpb << endl;
		} // End of 'if' for writing
					
	} // End of loop over time	
} // End of function
