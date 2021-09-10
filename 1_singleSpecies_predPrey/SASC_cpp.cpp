//============================================================================
// Name        : SASC
// Author      : Emanuel A. Fronhofer
// Date	       : 2021
//============================================================================

/*
	Copyright (C) 2021  Emanuel A. Fronhofer

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
*/

#include <iostream>
#include <cstdlib>								//standard C library
#include <ctime>								//access system time library
#include <fstream>								//file streaming library
#include <string>								//string library included
#include <sstream>								//string streaming for reading numbers from

#include <vector>
#include <cmath>								//standard math library

#include <gsl/gsl_rng.h>						//gsl random number generator
#include <gsl/gsl_randist.h>					//gsl random distributions
#include <gsl/gsl_statistics.h>					//gsl some statistical methods
#include <gsl/gsl_statistics_double.h> 			//gsl some double statistical methods
#include <gsl/gsl_sort_double.h> 				//gsl sort double arrays

#include <algorithm>

using namespace std;

#include "include/procedures.h"					//procedure simplifications
#include "include/classes.h"					//class definitions

//_____________________________________________________________________________
//------------------------------------------------------------ global variables
unsigned int time_max;															//actual time in simulation
int max_runs;																	//maximal number of runs

double e0;																		//consumer efficiency
double h0;																		//consumer handling time
double searchEff0;																//consumer search efficiency
double d0_cons;																	//consumer death rate

double alpha;																	//resource competition coefficient
double B0_res;																	//resource birth rate
double d0_res;																	//resource death rate

unsigned int cons0;																//starting conditions consumer
unsigned int res0;																//starting conditions resources

unsigned int act_N_resources;													//resource pop sizes per patch
unsigned int act_N_consumers;													//consumer pop sizes per patch

unsigned int sum_all_indivs;													//community size

double cbirth = 0;																//maximum rate constant for birth
double cdeath = 0;																//maximum rate constant for death

int drf;																		//resource density regulation function;
int fr;																			//consumer functional response
double eta_Allee;																//Allee effect parameter 1 (N=0 effect) for logistic and BH model
double gamma_Allee;																//Allee effect parameter 2 (effect on higher densities) for logistic and BH model

double apecc;																	//additive process error constant for consumers

//_____________________________________________________________________________
//------------------------------------------------------------------ procedures

//------------------------------------------------------------- read parameters
void readParameters(){
	ifstream parinfile("input/parameters.in");							//parameter input file
	string buffer;
	istringstream is;

	getline(parinfile,buffer); getline(parinfile,buffer);
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> time_max;																		//simulation time
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> max_runs;																		//number of replicates
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> e0;																			//consumer efficiency
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> h0;																			//consumer handling time
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> searchEff0;																	//consumer search efficiency
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> d0_cons;																		//consumer death rate
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> alpha;																		//resource competition coefficient
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> B0_res;																		//resource birth rate
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> d0_res;																		//resource death rate
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> eta_Allee;																	//strength of Allee effect 1
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> gamma_Allee;																	//strength of Allee effect 2
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> cons0;																		//starting conditions consumer
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> res0;																			//starting conditions resources
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> drf;																			//resource density regulation function
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> fr;																			//consumer functional response
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> apecc;																		//additive process error constant for consumers
	parinfile.close();
}


//------------------------------------------------------------------------ initialize
void initialize(){

	// initialize the consumer population size
	act_N_consumers = cons0;

	// initialize the resource population size
	act_N_resources = res0;
	
	// get population size sum
	sum_all_indivs = act_N_consumers + act_N_resources;
}

//------------------------------------------------------------------------ calculation of maximum rate constants
void calc_max_rate_consts(){

	// declare variables
	vector <double> cbirth_all;
	cbirth_all.clear();
	vector <double> cdeath_all;
	cdeath_all.clear();

	double help_c = 0;
	
	// birth rate consumers
	if(cons0 != 0){
		switch (fr) {
			case 1:
				// type 1 FR
				help_c = e0*searchEff0 * double(act_N_resources) + apecc;
				cbirth_all.push_back(help_c);
				break;
			case 2:
				// type 2 FR
				help_c = e0*searchEff0 * double(act_N_resources) /(1+searchEff0*h0 * double(act_N_resources)) + apecc;
				cbirth_all.push_back(help_c);
				break;
			default:
				cout << "Functional response (input parameter: fr) not correctly defined!" << endl;
				exit(0);
		}

		// death rate consumers
		cdeath_all.push_back(d0_cons + apecc);
	}

	// birth rate resources
	switch (drf) {
		case 0:
			// r-alpha logistic model
			help_c = B0_res - alpha * double(act_N_resources);
			cbirth_all.push_back(help_c);
			break;
		case 1:
			// continuous-time Beverton-Holt model
			help_c = B0_res/(1+ alpha * double(act_N_resources));
			cbirth_all.push_back(help_c);
			break;
		case 2:
			// r-alpha logistic model with Allee effect (allee effect is a mortality term; see there!)
			help_c = (B0_res - alpha * double(act_N_resources));// - eta_Allee/(1 + gamma_Allee * double(act_N_resources))) ;
			cbirth_all.push_back(help_c);
			break;
		case 3:
			// continuous-time Beverton-Holt model with Allee effect (allee effect is a mortality term; see there!)
			help_c = B0_res/(1+ alpha * double(act_N_resources));// - eta_Allee/(1 + gamma_Allee * double(act_N_resources)));
			cbirth_all.push_back(help_c);
			break;
		default:
			cout << "Density-regulation function (input parameter: drf) not correctly defined!" << endl;
			exit(0);
	}

	// calculate once death of resources by consumers
	double consumer_driven_resource_death_rate = 0;

	if (cons0 != 0){
		switch (fr) {
			case 1:
				// type 1 FR
				consumer_driven_resource_death_rate = searchEff0 * double(act_N_consumers);
				break;
			case 2:
			// type 2 FR
				consumer_driven_resource_death_rate = searchEff0 * double(act_N_consumers) /(1+searchEff0*h0 * double(act_N_resources));
				break;
			default:
				cout << "Functional response (input parameter: fr) not correctly defined!" << endl;
				exit(0);
		}
	}
	// death rate resources
		switch (drf) {
		case 0:
			// r-alpha logistic model
			help_c = d0_res + consumer_driven_resource_death_rate;
			cdeath_all.push_back(help_c);
			break;
		case 1:
			// continuous-time Beverton-Holt model (same as case 0 as density regulation is assumed to work on birth rates
			help_c = d0_res + consumer_driven_resource_death_rate;
			cdeath_all.push_back(help_c);
			break;
		case 2:
			// r-alpha logistic model with Allee effect
			//help_c = d0_res + consumer_driven_resource_death_rate;
			help_c = d0_res + consumer_driven_resource_death_rate + eta_Allee/(1 + gamma_Allee * double(act_N_resources));
			cdeath_all.push_back(help_c);
			break;
		case 3:
			// continuous-time Beverton-Holt model with Allee effect
			//help_c = d0_res + consumer_driven_resource_death_rate;
			help_c = d0_res + consumer_driven_resource_death_rate + eta_Allee/(1 + gamma_Allee * double(act_N_resources));
			cdeath_all.push_back(help_c);
			break;
		default:
			cout << "Density-regulation function (input parameter: drf) not correctly defined!" << endl;
			exit(0);
	}
	
	

	// get the maxima
	cbirth = max(cbirth_all);
	cdeath = max(cdeath_all);
}


//------------------------------------------------------------------------- potential consumer birth
bool consumer_birth(){
	// set event to zero
	bool event = 0;

	// calculate birth rate
	double focal_birth = 0;

	switch (fr) {
		case 1:
			// type 1 FR
			focal_birth = e0*searchEff0*double(act_N_resources) + apecc;
			break;
		case 2:
			// type 2 FR
			focal_birth = e0*searchEff0*double(act_N_resources)/(1+searchEff0*h0*double(act_N_resources)) + apecc;
			break;
		default:
			cout << "Functional response (input parameter: fr) not correctly defined!" << endl;
			exit(0);
	}

	if(ran() < (focal_birth/cbirth)){
		//execute birth
		act_N_consumers = act_N_consumers + 1;
		// increase community size
		sum_all_indivs = sum_all_indivs + 1;

		// set event to true
		event = 1;
	}
	return(event);
}

//------------------------------------------------------------------------- potential resource birth
bool resource_birth(){
	// set event to zero
	bool event = 0;
	
	// calculate birth rate
	double focal_birth = 0;
	
		// birth rate resources
	switch (drf) {
		case 0:
			// r-alpha logistic model
			focal_birth = B0_res - alpha * double(act_N_resources);
			break;
		case 1:
			// continuous-time Beverton-Holt model
			focal_birth = B0_res/(1+ alpha * double(act_N_resources));
			break;
		case 2:
			// r-alpha logistic model with Allee effect
			focal_birth = (B0_res - alpha * double(act_N_resources)); //- eta_Allee/(1 + gamma_Allee * double(act_N_resources))) ;
			break;
		case 3:
			// continuous-time Beverton-Holt model with Allee effect
			focal_birth = B0_res/(1+ alpha * double(act_N_resources)); //- eta_Allee/(1 + gamma_Allee * double(act_N_resources))) ;
			break;
		default:
			cout << "Density-regulation function (input parameter: drf) not correctly defined!" << endl;
			exit(0);
	}
	
	
	if(ran() < (focal_birth/cbirth)){
		// execute birth
		act_N_resources = act_N_resources+1;
		// increase community size
		sum_all_indivs = sum_all_indivs + 1;

		// set event to true
		event = 1;
	}
	return(event);
}


//------------------------------------------------------------------------potential consumer death
bool consumer_death(){
	//set event to zero
	bool event=0;
	
	// calculate death rate
	double focal_death = d0_cons + apecc;
	if(ran() < (focal_death/cdeath)){
		//execute death
		act_N_consumers = act_N_consumers -1;
		// decrease community size
		sum_all_indivs = sum_all_indivs - 1;
		
		// set event to true
		event = 1;
	}
	return(event);
}

//------------------------------------------------------------------------ potential resource death
bool resource_death(){
	// set event to zero
	bool event = 0;
	
		// calculate once death of resources by consumers
	double consumer_driven_resource_death_rate = 0;

	if (cons0 != 0){
		switch (fr) {
			case 1:
				// type 1 FR
				consumer_driven_resource_death_rate = searchEff0 * double(act_N_consumers);
				break;
			case 2:
				// type 2 FR
				consumer_driven_resource_death_rate = searchEff0 * double(act_N_consumers) /(1+searchEff0*h0 * double(act_N_resources));
				break;
			default:
				cout << "Functional response (input parameter: fr) not correctly defined!" << endl;
				exit(0);
		}
	}
	
	double focal_death = 0;
	if(act_N_consumers == 0){
		// calculate death rate
		focal_death = d0_res;
		if (drf == 2 || drf == 3) focal_death = d0_res + eta_Allee/(1 + gamma_Allee * double(act_N_resources));
	}else{
		// calculate death rate
		switch (drf) {
		case 0:
			// r-alpha logistic model
			focal_death = d0_res + consumer_driven_resource_death_rate;
			break;
		case 1:
			// continuous-time Beverton-Holt model (same as case 0 as density regulation is assumed to work on birth rates
			focal_death = d0_res + consumer_driven_resource_death_rate;
			break;
		case 2:
			// r-alpha logistic model with Allee effect
			focal_death = d0_res  + consumer_driven_resource_death_rate + eta_Allee/(1 + gamma_Allee * double(act_N_resources));
			//focal_death = d0_res  + consumer_driven_resource_death_rate;
			//cout << d0_res  + consumer_driven_resource_death_rate << "     " << d0_res  + consumer_driven_resource_death_rate + eta_Allee/(1 + gamma_Allee * double(act_N_resources)) << endl;
			break;
		case 3:
			// continuous-time Beverton-Holt model with Allee effect
			//focal_death = d0_res + consumer_driven_resource_death_rate ;
			focal_death = d0_res + consumer_driven_resource_death_rate + eta_Allee/(1 + gamma_Allee * double(act_N_resources));
			break;
		default:
			cout << "Density-regulation function (input parameter: drf) not correctly defined!" << endl;
			exit(0);
		}
	}

	if(ran() < (focal_death/cdeath)){
		// execute death
		act_N_resources = act_N_resources-1;
		// decrease community size
		sum_all_indivs = sum_all_indivs - 1;
		
		// set event to true
		event = 1;
	}
	return(event);
}

//_____________________________________________________________________________
//------------------------------------------------------------------------ main

int main() {
	// random number generator
	//specify_rng(time(NULL));
	specify_rng(RS);

	//read parameters for all simulation runs
	readParameters();

	// repeat loop
	for (int actrun = 0; actrun < max_runs; ++actrun) {

		// initialize the community
		initialize();

		// define some time variables locally
		double act_time = 0;																//present time
		unsigned long int time_cnt = 0;														//counter for time
		double cum_time_interval = 0;														//cumulative time interval

		// calculate the maximum rate constants
		calc_max_rate_consts();

		//-------------------------------------------------------------------------------------------------------------------------------
		// output file: community dynamics
		stringstream outputindiv_path_stream;
		outputindiv_path_stream << "output/output_run" << actrun << ".out";
		string outputindiv_path = outputindiv_path_stream.str();
		ofstream outputindiv(outputindiv_path.c_str());
		// headers
		outputindiv << "run" << "    " << "time" << "    " << "cons" << "    " << "res" << endl;
		//-------------------------------------------------------------------------------------------------------------------------------

		// time loop; the appropriateness of it_max has to be checked!
		for (unsigned int actit = 0; actit < it_max; ++actit) {

			// calculate time increment
			double lambda = (cbirth+cdeath)*double(sum_all_indivs);

			// end simulation if lambda = 0, i.e. nothing is happening
			if (lambda == 0){
				act_time = time_max;
				break;
			}

			// calculate time interval
			double time_interval = expo(lambda);
			cum_time_interval = cum_time_interval + time_interval;

			// reset event to false
			bool event = 0;

			// randomly choose an event
			double event_ran = ran();

			// now choose the event according to the maximum rate constants (birth or death)
			if (event_ran < (cbirth/(cbirth+cdeath))) {
				// potential birth
				// decide on consumer or resources as a weighted lottery of population size
				if(ran() < double(act_N_consumers)/double(sum_all_indivs)){
					// consumers
					event = consumer_birth();
				}else{
					//resources
					event = resource_birth();
				}
			} else {
				// death
				// decide on consumer or resources as a weighted lottery of population size
				if(ran() <  (double(act_N_consumers)/double(sum_all_indivs))){
					//consumers
					event = consumer_death();
				}else{
					//resources
					event = resource_death();
				}
			}

			// update time and community if an event was executed
			if (event == 1){
				// time
				time_cnt = time_cnt + 1;
				act_time = act_time + cum_time_interval;
				cum_time_interval = 0;

				// recalculate maximal rate constants
				calc_max_rate_consts();

				//--------------------------------------------------------------------------------------------
				// some output to file
				outputindiv << actrun << "    " << act_time << "    " << act_N_consumers << "    " << act_N_resources << endl;
				// some output to console
				cout << actrun << "    " << act_time << "    " << act_N_consumers << "    " << act_N_resources << endl;
				//----------------------------------------------------------------------------------------------
			}

			// end simulation if the maximal time has been reached
			if(act_time >= time_max){
				break;
			}
		}
		// close output file
		outputindiv.close();
	}
	return 0;
}
