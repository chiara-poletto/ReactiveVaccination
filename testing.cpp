/**
    Simulation of an ABM for COVID-19
    @file screening.h
    @author Pierre-Yves Boelle, Chiara Poletto
    @version 1.0 2020-04-16
    @license GPL-3.0-or-later
*/

#include <stdio.h>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <array>
#include <string>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <set>
#include <unordered_set>
#include <tuple>
#include <map>
#include <cerrno>
#include <cstring>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <utility>

#include "defs.h"
#include "routine_miscellanea.hpp"

using namespace std;

#define DEBUG

void initTesting(Testing & testing, Population & city) {
    // Nodes (once detected) for which manual contact tracing was done (1= mct done)
    testing.ct_done.resize(city.N_use);
    fill(testing.ct_done.begin(),testing.ct_done.end(), 0);
    // Nodes (once detected) for which manual contact tracing was done (1= mct done)
    testing.found_in.resize(TESTING_MODE, vector<int>(city.N_use));

    for (int loc=0; loc<TESTING_MODE;loc++) {
    	std::fill(testing.found_in[loc].begin(),testing.found_in[loc].end(), 0);
    	testing.test_in[loc]=0;
    }
}

/**
 * case detection
 * case are detected if symptomatic (12/912) or asymptomatic(13/913) with appropriate rate
 * returns :	1 if detected symptomatic
 * 				2 if detected asymptomatic
 * 				-1 if undetected
 */

void caseDetection(Population & city, Params &params, default_random_engine & generator) {
	uniform_real_distribution<double> distribution(0,1);

	// run through nodes to determine if they will be detected
	// if !detected yet and will_be_detected=0 we check if it is going to be detected
	// then set will_be_detected!=0 so we won't check again
	for (int i{0} ; i<city.N_use ; ++i)
	{
		// not detected or scheduled for detection
		if ((city.detected[i] == 0) && (city.will_be_detect[i]==0)) {
			//symptomatic & will be detected
			if ((city.status[i]==12) || (city.status[i]==912)) {
				if (distribution(generator)<params.will_be_d_symp) {
					// scheduled for symptomatic
					city.will_be_detect[i]=1;
				} else {
					city.will_be_detect[i]=-1;
				}
			}
			// asymptomatic and will be detected
			else if ((city.status[i]==13) || (city.status[i]==913)) {
				if (distribution(generator)<params.will_be_d_asymp) {
					city.will_be_detect[i]=2;
				}
				// will not be detected
				else {
					city.will_be_detect[i]=-1;
				}
			}
		}
	}
}

void findContactsInLayersAtI(
		set<int> & encounters, // store results
		int i, //index of candidate
		int it, // current time (iteration)
		int location, // look for contacts in the household only (layer 0)
		Networks & contacts, // contacts
		Population & city, // to look at isolation status
		Params & params //params
) {
	//identify contacts and returns filled encounters
	int minlayer;
	int maxlayer;
	if (location == HOUSEHOLD ) {
		minlayer=0;
		maxlayer=1;
	} else {
		minlayer=1;
		maxlayer=5;
	}

	// Go back in time
	for (int k{0}; k<params.tracing_lookback; ++k) {
        // if this back is ok
		if (it > k )
        {
			// this is the index of the network at the time k days before
            int timeAtContact= (it-k);
            int net_timeAtContact= timeAtContact % int(params.number_networks);
            for (int layer = minlayer; layer < maxlayer; layer ++) {
            	for (int j: contacts.layers[layer][net_timeAtContact][i]) {
            		if (location==HOUSEHOLD) {
            				// contacts are always inserted in the household
            			encounters.insert(j);
            		} else  // not in Household
            			// contacts are only added if
            			if (city.isolation[timeAtContact][j]==0) {// j is not isolated
            				if (city.isolation[timeAtContact][i]==0) {// i is is no isolated
								if (city.onset_time[i]>=timeAtContact || city.onset_time[i]==0)  {// i is is not symptomatic
									encounters.insert(j);
								}
            				}
            			}
            	}
            }
        }
	}
}

/**
 * detection is done in the place if not already planned
 * return 1 if positive, 0 otherwise
 *
 */
int testThisPerson(int j, int loc, Population & city, Testing & testing, Compartments & upd) {
	// in Household = 1,  : params.p_c_HH
	// in Household = 0, acquaintance  : params.p_rA
	// in Household = 0, not acquaintance  : params.p_rS
	uniform_real_distribution<double> distribution(0,1);
	// if contact was not already detected or scheduled for detection

	if ((city.detected[j]==0) && (city.will_be_detect[j]==0))
	{
		//symptomatic
		if (city.status[j]==12 || city.status[j]==912 || city.status[j]==101 || city.status[j]==9101)
		{
			upd.detected.push_back(j); // this will put detected to 1
			//it is detected and traced
			city.will_be_detect[j]=1; // this is done in place
			//if (loc == HOUSEHOLD || loc==TRACING) testing.ct_done[j]=1; // if not found in household, then contact tracing should be done
			testing.found_in[loc][j]=1;
			testing.test_in[loc]++;
			return(1);
		}
		//asymptomatic
		else if (city.status[j]==13 || city.status[j]==913 || city.status[j]==102 || city.status[j]==9102)
		{
			upd.detected.push_back(j); // updated later
			//it is detected and traced
			city.will_be_detect[j]=2;
            //if (loc==HOUSEHOLD || loc==TRACING) testing.ct_done[j]=1; // if not found in household, then contact tracing should be done
			testing.found_in[loc][j]=1;
			testing.test_in[loc]++;
			return(1);
		}
	}
	return(0);
}

/**
 *
 * we look for contacts of the person in its household over the last days (params.tracing_lookback)
 * then persons are isolated if they agree to
 *
 *
 */

void isolateContactsInHousehold (int i, Testing & testing, Population & city, Networks & contacts, Compartments & upd, Params & params,
		default_random_engine & generator, SimulationFiles & simulationFiles,
		int real, int it) {
	// HH contacts are isolated, HH infected contacts are detected (detected=1)
	// Household isolation
	// Initialisation of immunity at random
	uniform_real_distribution<double> distribution(0,1);

	set<int> encounters;
	findContactsInLayersAtI(encounters, i, it, HOUSEHOLD, contacts, city, params);

	// Ask to each contact j
	for (int j: encounters) {
		// isolate this contact in the household
		double probIsolate=params.p_c_HH;
		if (distribution(generator)<probIsolate)
		{
			//isolate contact - if not isolated
			//if (city.isolation[it-1][j] ==0 ) { // this leads to iso_time set at time of isolation
            upd.isolated.push_back(j); // isolate or reset isolation
			//}
			testThisPerson(j,HOUSEHOLD, city, testing, upd);
			// Print out isolation details
			if (params.print_isolations==1) { // records details of the node to write in output file
				simulationFiles.isolation_file << i << "," << city.Age[i] << "," << city.status[i] << "," << j << "," << city.Age[j] << "," << city.status[j] << "," << it << "," << real << endl;
			}
		}
	}
}


/**
 * CASE ISOLATION && HH ISOLATION
 * we go through each person and if it is detected now,
 	// if a node that will be detected is really detected at this time step
	// - it is isolated
	// - its HH contacts are isolated - those willing to isolate
	// - among the HH contacts that are isolated the infected are detected
	// except for the ones with will_be_detected=1 or 2 that will be detected as index nodes and enable contact tracing shortly after
	// - contact tracing is limited to contact in the last days
	// on the infected detected among the isolated
 *
 */
void caseAndHouseholdIsolation(int real, int it, Testing & testing, Population & city, Networks & contacts ,
		Compartments & upd,
		Params & params, default_random_engine & generator,SimulationFiles & simulationFiles ) {
	//
	uniform_real_distribution<double> distribution(0,1);
	for (int i{0} ; i<city.N_use ; ++i)
	{
		if (city.detected[i]==0) {  //not detected already
			if (// scheduled for detection symptomatic
					((city.will_be_detect[i]==1) && (distribution(generator)<params.r_d_symp)) ||
					// scheduled for detection asymptomatic
					((city.will_be_detect[i]==2) && (distribution(generator)<params.r_d_asymp))) {
				// node i will be detected now
				upd.detected.push_back(i); // city.detected[i]=1; marked for update
				//if (city.isolation[it-1][i] == 0) { //marked for isolation if not isolated
                upd.isolated.push_back(i); //isolate or reseat isolation
				//}
				// HH isolation is done:
				if (params.HHisolation_inplace==1) {
					isolateContactsInHousehold( i, testing, city, contacts, upd,
							params, generator, simulationFiles, real,  it);
				}
			}
		}
	}
}

/**
 * MANUAL CONTACT TRACING
 		// Run through nodes. If they have detected = 1:
		// - see if manual contact tracing will be done (ct_done==0)
		// - if manual contact tracing is done : among contacts that are isolated the infected are detected
		// except for the ones with will_be_detected=1 or 2, these will be detected as index nodes
		// with respective rate, to enable contact tracing
		// Contact tracing may be done on the people who are detected at this time
		 * later
 *
 */
void manualContactTracing(int real, int it, Testing & testing, Population & city, Networks & contacts, Compartments & upd,
		Params & params, default_random_engine & generator,SimulationFiles & simulationFiles )
{
	//
	if (params.manualCT_inplace==1) {

		//cout << "Manual contact tracing" << endl;
		uniform_real_distribution<double> distribution(0,1);

		for (int i{0} ; i<city.N_use ; ++i)
		{
			if (
					// detected
					(city.detected[i]==1) &&
					//not already ct
					(testing.ct_done[i]==0) &&
					// contact traced now
					(distribution(generator)<params.r_mCT ))
			{
				//update ct_done index
				testing.ct_done[i]= 1;
				// Find contacts of past days outside the household
				int timeIsolationAtI= city.iso_time[i];
				set<int> encounters;
				findContactsInLayersAtI(encounters, i, timeIsolationAtI, TRACING, contacts, city, params);
				//loop over contacts

				for (int j: encounters) {
					// define close contact
					int isNotAcquaintance = (contacts.acquaintances[i].find(j) == contacts.acquaintances[i].end());
					if (isNotAcquaintance==1) { // it is not an acquaintance
						double probIsolate=params.p_rS; //accept if not acquaintance
						if (distribution(generator)<probIsolate) {
							//isolate contact
							//if not already isolated
							//if (city.isolation[it-1][j]==0) {//city.isolation[it][j]=1;
                            upd.isolated.push_back(j); // isolate or reset isolation
							//}
							testThisPerson(j, TRACING, city,testing,upd);
							if (params.print_isolations==1) { // records details of the node to write in output file
								simulationFiles.isolation_file <<  i << "," << city.Age[i] << "," << city.status[i] << "," << j << "," << city.Age[j] << "," << city.status[j] << "," << it << "," << real << endl;
							}
						}
					} else if (isNotAcquaintance==0) { //it is an acquaintance
						double probIsolate=params.p_rA; //accept if acquaintance
						if (distribution(generator)<probIsolate) {
							//isolate contact - note that this updates isolation time ()
							upd.isolated.push_back(j); //city.isolation[it][j]=1;
							testThisPerson(j, TRACING, city,testing,upd);
							if (params.print_isolations==1) { // records details of the node to write in output file
								simulationFiles.isolation_file <<  i << "," << city.Age[i] << "," << city.status[i] << "," << j << "," << city.Age[j] << "," << city.status[j] << "," << it << "," << real << endl;
							}
						}
					}
				}
			}
		}
	}
}


void testingStep(int real, int it, Testing & testing, Population & city, Networks &contacts,
		Compartments &upd, Params & params,
		default_random_engine & generator,SimulationFiles & simulationFiles) {

	// ISOLATION && HH ISOLATION && MANUAL CONTACT TRACING

	// sets who will_be_detected with the new status
	// some infections are never detected, other with rates depending on symptoms
	caseDetection(city, params, generator);
    
	// case and Household isolation
	// look whether index case is detected now, isolate in the household
    caseAndHouseholdIsolation(real, it, testing, city, contacts, upd, params,generator, simulationFiles);

    // manual contact tracing
    manualContactTracing(real, it, testing, city,  contacts, upd, params, generator, simulationFiles);
    
}

