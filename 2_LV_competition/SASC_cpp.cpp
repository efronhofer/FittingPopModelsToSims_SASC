//============================================================================
// Name        : SASC
// Author      : Emanuel A. Fronhofer
// Version     : v0
// Date	       : April 2019
//============================================================================

/*
	Copyright (C) 2019  Emanuel A. Fronhofer

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

double alpha_11;																//resource1 intaspecific competition coefficient
double alpha_22;																//resource2 intaspecific competition coefficient

double alpha_21;																//resource1 interspecific competition coefficient
double alpha_12;																//resource2 interspecific competition coefficient

double B0_res_1;																//resource1 birth rate
double d0_res_1;																//resource1 death rate

double B0_res_2;																//resource2 birth rate
double d0_res_2;																//resource2 death rate

unsigned int res0_1;															//starting conditions resources1
unsigned int res0_2;															//starting conditions resources2

unsigned int act_N_resources_1;													//resource1 pop size
unsigned int act_N_resources_2;													//resource2 pop size

unsigned int sum_all_indivs;													//community size

double cbirth = 0;																//maximum rate constant for birth
double cdeath = 0;																//maximum rate constant for death

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
	is >> alpha_11;																		//resource1 intaspecific competition coefficient
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> alpha_22;																		//resource2 intaspecific competition coefficient
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> alpha_21;																		//resource1 interpecific competition coefficient
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> alpha_12;																		//resource2 interpecific competition coefficient
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> B0_res_1;																		//resource1 birth rate
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> B0_res_2;																		//resource2 birth rate
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> d0_res_1;																		//resource1 death rate
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> d0_res_2;																		//resource2 death rate
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> res0_1;																			//starting conditions resources1
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> res0_2;																			//starting conditions resources1

	parinfile.close();
}


//------------------------------------------------------------------------ initialize
void initialize(){

	// initialize the resource population size
	act_N_resources_1 = res0_1;
	act_N_resources_2 = res0_2;

	// get population size sum
	sum_all_indivs =  act_N_resources_1 + act_N_resources_2;
}

//------------------------------------------------------------------------ calculation of maximum rate constants
void calc_max_rate_consts(){

	// declare variables
	vector <double> cbirth_all;
	cbirth_all.clear();
	vector <double> cdeath_all;
	cdeath_all.clear();

	double help_c = 0;

	// birth rate resources1
	help_c = B0_res_1 - alpha_11 * double(act_N_resources_1) - alpha_12 * double(act_N_resources_2);
	cbirth_all.push_back(help_c);

	// birth rate resources2
	help_c = B0_res_2 - alpha_22 * double(act_N_resources_2) - alpha_21 * double(act_N_resources_1);
	cbirth_all.push_back(help_c);

	// death rate resources1
	help_c = d0_res_1;
	cdeath_all.push_back(help_c);

	// death rate resources2
	help_c = d0_res_2;
	cdeath_all.push_back(help_c);

	// get the maxima
	cbirth = max(cbirth_all);
	cdeath = max(cdeath_all);
}


//------------------------------------------------------------------------- potential resource1 birth
bool resource1_birth(){
	// set event to zero
	bool event = 0;
	
	// calculate birth rate
	double focal_birth = 0;

	// birth rate resources
	focal_birth = B0_res_1 - alpha_11 * double(act_N_resources_1) - alpha_12 * double(act_N_resources_2);

	if(ran() < (focal_birth/cbirth)){
		// execute birth
		act_N_resources_1 = act_N_resources_1+1;
		// increase community size
		sum_all_indivs = sum_all_indivs + 1;

		// set event to true
		event = 1;
	}
	return(event);
}

//------------------------------------------------------------------------- potential resource2 birth
bool resource2_birth(){
	// set event to zero
	bool event = 0;
	
	// calculate birth rate
	double focal_birth = 0;

	// birth rate resources
	focal_birth = B0_res_2 - alpha_22 * double(act_N_resources_2) - alpha_21 * double(act_N_resources_1);

	if(ran() < (focal_birth/cbirth)){
		// execute birth
		act_N_resources_2 = act_N_resources_2+1;
		// increase community size
		sum_all_indivs = sum_all_indivs + 1;

		// set event to true
		event = 1;
	}
	return(event);
}

//------------------------------------------------------------------------ potential resource1 death
bool resource1_death(){
	// set event to zero
	bool event = 0;

	double focal_death = d0_res_1;

	if(ran() < (focal_death/cdeath)){
		// execute death
		act_N_resources_1 = act_N_resources_1-1;
		// decrease community size
		sum_all_indivs = sum_all_indivs - 1;
		
		// set event to true
		event = 1;
	}
	return(event);
}

//------------------------------------------------------------------------ potential resource2 death
bool resource2_death(){
	// set event to zero
	bool event = 0;

	double focal_death = d0_res_2;

	if(ran() < (focal_death/cdeath)){
		// execute death
		act_N_resources_2 = act_N_resources_2-1;
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
		outputindiv << "run" << "    " << "time" << "    " << "res1" << "    " << "res2" << endl;
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
				// decide on res1 or res2 as a weighted lottery of population size
				if(ran() < double(act_N_resources_1)/double(sum_all_indivs)){
					// consumers
					event = resource1_birth();
				}else{
					//resources
					event = resource2_birth();
				}
			} else {
				// death
				// decide on res1 or res2 as a weighted lottery of population size
				if(ran() <  (double(act_N_resources_1)/double(sum_all_indivs))){
					//consumers
					event = resource1_death();
				}else{
					//resources
					event = resource2_death();
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
				outputindiv << actrun << "    " << act_time << "    " << act_N_resources_1 << "    " << act_N_resources_2 << endl;
				// some output to console
				cout << actrun << "    " << act_time << "    " << act_N_resources_1 << "    " << act_N_resources_2 << endl;
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
