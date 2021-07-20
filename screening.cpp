/**
    Simulation of an ABM for COVID-19
    @file screening.cpp
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
#include "testing.h"


using namespace std;



void initScreening(Screening & screening) {
    // List of Workplaces/Schools where testing campaign was done
	screening.WPs_regular_screening.clear();
    // List of Workplaces/Schools where testing campaign was done
	screening.Ss_regular_screening.clear();
	    // List of Workplaces/Schools where testing campaign was done
	screening.HHs_regular_screening.clear();
    //date regular screening
	screening.date_regular_screening.clear();
    //date regular screening
	screening.date_reactive_screening.clear();
}


void initRegularScreeningInPlace(set<string> &chosenPlaces, set<string> &candidatePlaces,
		Places & places,	double probScreening, int intervalScreening ,
		Screening & screening, default_random_engine & generator, SimulationFiles & sf) {
	// choose the places that will be submitted to screening
	uniform_real_distribution<double> distribution(0,1);
	uniform_int_distribution<int> distribution_int(0,intervalScreening-1);
	// go through all candidate places
	int too_young=0;
	int too_small=0;
	int inserted=0;
	int discarded =0;
	for (string idPlace: candidatePlaces) {
		if (places.too_young_to_test.find(idPlace) == places.too_young_to_test.end()) { // not found in too young
			if (places.too_small_to_test.find(idPlace) == places.too_small_to_test.end()) { // not found in not too small
				if (distribution(generator)<probScreening) {
					chosenPlaces.insert(idPlace);
							// choose starting date
					inserted++;
					int start=distribution_int(generator);
					screening.date_regular_screening.insert({idPlace,start});
					sf.miscellaneous_file << idPlace<<" "<< start <<" ";
				} else discarded++;
			} else too_small++;
		} else too_young++;
	}
	sf.miscellaneous_file <<endl;
	sf.miscellaneous_file << " total " << candidatePlaces.size() << " inserted " << inserted << " excluded " << discarded << " (not willing) " << too_young << " (too_young) " << too_small << " (too_small)" <<endl;
}
void initRegularScreening(Screening & screening, Places & places, Params & params, default_random_engine generator, SimulationFiles & sf) {

	// in WPs
	if (params.regular_screening_WP_inplace==1) {
		sf.miscellaneous_file << "WPs : "<<endl;
		initRegularScreeningInPlace(screening.WPs_regular_screening, places.WPs, places,
				params.f_WP_screening, params.interval_screening,  screening, generator,sf);
	}
	// in Ss
	if (params.regular_screening_S_inplace==1) {
		sf.miscellaneous_file << "Ss : "<<endl;
		initRegularScreeningInPlace(screening.Ss_regular_screening, places.Ss, places,
				params.f_S_screening,params.interval_screening,screening,generator,sf);
	}
	// in HHs
	if (params.regular_screening_HH_inplace==1) {
		sf.miscellaneous_file << "HHs : "<<endl;
		initRegularScreeningInPlace(screening.HHs_regular_screening, places.HHs, places,
				params.f_HH_screening,params.interval_screening,screening, generator,sf);
	}
	sf.miscellaneous_file << "inserted "<< screening.date_regular_screening.size() <<" discarded " << places.too_small_to_test.size() << " (too small) and " << places.too_young_to_test.size() << " (too young) for screening"<<endl;
}





int screenThisPlace(int typeScreening, string & idPlace,  map<string,vector<int>> & mapPlaceNodes, double probTesting,
		Testing & testing, Population & city, Networks & contacts, Compartments & upd,
		Params & params, default_random_engine & generator, SimulationFiles & simulationFiles, int real, int it) {

	uniform_real_distribution<double> distribution(0,1);

	int tests_done=0;
	// find people in the place; test them
	// if they are positive, isolate Household; then invoke contact tracing as well
	vector<int> vec; //
	vec= mapPlaceNodes[idPlace]; // get index of people to test

	// Loop on nodes belonging to the place
	for (int j: vec) {
		//aged >min
		if (city.Age[j]>params.age_min_testing) {
			//	if in household, can be tested if isolated; othrewise tested only if not isoalted
			if ((typeScreening==HOUSEHOLD) || ((typeScreening != HOUSEHOLD ) && (city.isolation[it-1][j]==0))) {
				//agrees to be tested
				if (distribution(generator)<probTesting) {
					//resultat_test is 1 if infected , 0 otherwise
					int resultat_test = testThisPerson(j, typeScreening, city,testing,upd);
					tests_done++;
					// if it is tested positive, detected will be updated in next iteration and CT will be done
					// if the test result is positive, this person should be put in isolation now if it was not.
					if ((resultat_test == 1) && (city.isolation[it-1][j]==0)) {
						upd.isolated.push_back(j);
					}
					// we put household on isolation in the same round if the case is positve
					if ((params.HHisolation_inplace==1) && (resultat_test == 1)) {
							isolateContactsInHousehold( j, testing, city, contacts, upd,
								params, generator, simulationFiles, real,  it);
					}
				}
			}
		}
	}
	return(tests_done);
}

void regularScreening(int real, int it,
		Screening & screening,
		Testing & testing,
		Population & city,
		Places & places,
		Networks & contacts,
		Compartments & upd,
		Params & params,
		default_random_engine & generator,
		SimulationFiles & simulationFiles) {
	//in WPs that were chosen, conduct regular screening every x days
	if (params.regular_screening_HH_inplace == 1) {
		double probTesting=params.testing_proportion_HH; // individual prob agrees to test
		// we screen all people in the WPs, and then isolate those who are found
		for (string idPlace: screening.HHs_regular_screening) {
			// idPlace agrees to screening
			if ((it - screening.date_regular_screening[idPlace]) % params.interval_screening ==0) {
				// it is the date for screening
				int tests_here = screenThisPlace(HOUSEHOLD, idPlace, places.HH_nodes, probTesting, testing, city, contacts,
						upd, params, generator,simulationFiles,real,it);
				simulationFiles.out_regular_screening_file <<real<<","<<it<<"," << idPlace << "," <<tests_here <<","<<places.HH_nodes[idPlace].size()<<endl;
			}
		}
	}
	if (params.regular_screening_S_inplace == 1) {
		double probTesting=params.testing_proportion_S; // individual prob agrees to test
		// we screen all people in the WPs, and then isolate those who are found
		for (string idPlace: screening.Ss_regular_screening) {
			// idPlace agrees to screening
			if ((it - screening.date_regular_screening[idPlace]) % params.interval_screening ==0) {
				// it is the date for screening
				int tests_here=screenThisPlace(REGULAR_SCREENING, idPlace, places.S_nodes, probTesting, testing, city, contacts,
										upd, params, generator,simulationFiles,real,it);
				simulationFiles.out_regular_screening_file <<real<<","<<it<<"," << idPlace << "," <<tests_here <<","<<places.S_nodes[idPlace].size() <<endl;

			}
		}
	}
	if (params.regular_screening_WP_inplace == 1) {
		double probTesting=params.testing_proportion_WP; // individual prob agrees to test
		// we screen all people in the WPs, and then isolate those who are found
		for (string idPlace: screening.WPs_regular_screening) {
			// idPlace agrees to screening
			if ((it - screening.date_regular_screening[idPlace]) % params.interval_screening == 0) {
				// it is the date for screening
				int tests_here = screenThisPlace(REGULAR_SCREENING, idPlace, places.WP_nodes, probTesting, testing, city, contacts,
														upd, params, generator,simulationFiles,real,it);
				simulationFiles.out_regular_screening_file <<real<<","<<it<<"," << idPlace << "," <<tests_here <<","<< places.WP_nodes[idPlace].size() <<endl;
			}
		}
	}
}



int checkPlaceForReactiveScreening(int it, string & wps,  set<string> & listOfRegularPlaces, map<string, int> & date_reactive_screening,
		 Params &params, default_random_engine & generator) {
	// returns 1 :
	// - if place is not in regular screening
	// - if last screening more than params.interval_screening away
	// - with rate r_reactive
	// returns -1 if:
	// - regular screening in place
	// returns -2 if :
	// - last screening less than params.interval_screening away
	// returns 0 otherwise : it will be tested later
	uniform_real_distribution<double> distribution(0,1);

	if (listOfRegularPlaces.find(wps)==listOfRegularPlaces.end()) {// not in regular screening
		int max_date_for_last_screening = it - params.interval_screening;
		if (date_reactive_screening.find(wps)==date_reactive_screening.end()) {//not screened yet
			if (distribution(generator)<params.r_react_screen) { //enter with rate
				return(1);
			} else {
				return(0);
			}
		} else if (date_reactive_screening[wps] < max_date_for_last_screening) { //screened more than interval_screening ago
			if (distribution(generator)<params.r_react_screen) { //enter with rate
				return(1);
			} else {
				return(0);
			}
		} else { // screened not long ago
			return(-2);
		}
	} else {
		return(-1); // in regular screening
	}
	return(0);
}

/*
 *
 *
 */
void reactiveScreening(int real,int it,
		Screening & screening,
		Testing & testing,
		Population & city,
		Places & places,
		Networks & contacts,
		Compartments & upd,
		Params & params,
		default_random_engine& generator,
		SimulationFiles & simulationFiles) {

	int tests_today;
	tests_today=0;
    // WP and S to test now
    set<string> places_to_screen;
    set<string> places_to_delete;
    // If the number of cases in WP/S is higher than a threshold (th_cluster_screen)
    // I check that the testing is not done on a regular basis and not already done there too closely (interval_screening)
    // If not, with a certain delay (rate screen_rate) I do the testing (I update the date of last screening)
    if (params.reactive_screening_WP_inplace==1) {
    	for (string wps: screening.place_reactive_screening) { // this list is checked for size
    		// is is a WP or a school ?
    		if (places.WPs.find(wps) != places.WPs.end()) { // it is a WP
    			int insert_or_delete = checkPlaceForReactiveScreening(it, wps, screening.WPs_regular_screening,
    					screening.date_reactive_screening,
    					params,generator); // 1 for insertion, 0 to keep in the list, -1 or -2 for deletion
    			if ( insert_or_delete == 1) {
    				places_to_screen.insert(wps); // place is kept for screening NOW
    			} else if (insert_or_delete <0 ) {
    				places_to_delete.insert(wps); //place will not be screened (regular -1 | too close -2)
    	    		simulationFiles.out_reactive_screening_file << real<<","<<it<<"," << wps << "," <<insert_or_delete<<endl;
    			}
    		}
    	}
    	for (string wps : places_to_screen) { // places can be screened now
    		if (testing.test_in[REACTIVE_SCREENING]>params.reactive_tests_limit || tests_today> params.reactive_tests_daily)
    			break; // others will be done later
    		int tests_here = screenThisPlace(REACTIVE_SCREENING, wps, places.WP_nodes, params.testing_proportion_WP,
    				testing, city, contacts, upd, params, generator, simulationFiles,real, it);
    		tests_today = tests_today +tests_here;
    		screening.place_reactive_screening.erase(wps); // this place is done
    		screening.date_reactive_screening[wps]=it; // store date
    		simulationFiles.out_reactive_screening_file << real<<","<<it<<"," << wps << "," <<tests_here <<endl;
    	}
    	for (string wps : places_to_delete) { // places should be deleted from list
    		screening.place_reactive_screening.erase(wps);
    	}


    }
    places_to_screen.clear();
    places_to_delete.clear();
    if (params.reactive_screening_S_inplace==1) {
    	for (string wps:  screening.place_reactive_screening) {
    		if (places.Ss.find(wps) != places.Ss.end()) { // it is a School
    			int insert_or_delete = checkPlaceForReactiveScreening(it, wps,  screening.Ss_regular_screening,
    					screening.date_reactive_screening,
    				 params,generator);
    			if (insert_or_delete== 1) {
    				places_to_screen.insert(wps);
    			} else if (insert_or_delete < 0 ) {
    				places_to_delete.insert(wps);
    	    		simulationFiles.out_reactive_screening_file << real<<","<<it<<"," << wps << "," <<insert_or_delete<<endl;
    			}
    		}
    	}
    	for (string wps : places_to_screen) {
    		if (testing.test_in[REACTIVE_SCREENING]>params.reactive_tests_limit || tests_today> params.reactive_tests_daily)
    			break;
    		int tests_here = screenThisPlace(REACTIVE_SCREENING, wps, places.S_nodes, params.testing_proportion_S,
    				testing, city, contacts, upd, params, generator, simulationFiles,real, it);
    		tests_today= tests_today + tests_here;
    		screening.place_reactive_screening.erase(wps); // place was screened
    		screening.date_reactive_screening[wps]=it; // store date
    		simulationFiles.out_reactive_screening_file << real<<","<<it<<"," << wps << "," <<tests_here <<endl;
    	}
    	for (string wps : places_to_delete) { // places should be deleted from list
    		screening.place_reactive_screening.erase(wps);
    	}
    }
}


void screeningStep(int real, int it, Screening & screening, Testing & testing, Population & city, Networks &contacts, Places & places,
		Compartments &upd, Params & params,
		default_random_engine & generator,SimulationFiles & simulationFiles) {

    // regaularScreening
	regularScreening(real, it, screening, testing, city, places, contacts, upd,params, generator,simulationFiles);

	// reactiveScreening : in places where some cases have been identified
	reactiveScreening(real, it, screening, testing, city, places, contacts, upd,params, generator,simulationFiles);

}

