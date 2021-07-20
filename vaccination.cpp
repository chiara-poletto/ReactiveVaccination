/**
    Simulation of an ABM for COVID-19
    @file vaccination.cpp
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

void setVaccinationOpinion(Vaccination & vaccination, Population & city, Params & params,default_random_engine& generator) {
	uniform_real_distribution<double> distribution(0,1);

    for (int i{0}; i<city.N_use ;++i) {
    if(city.Age[i] > 65)
    {
        if (distribution(generator)> params.coverage_vacc_higher_65/100.)
        	vaccination.vaccination_opinion[i] = 1;
        else vaccination.vaccination_opinion[i] = 0;
    }
    else
        if (distribution(generator)> params.coverage_vacc_lower_65/100.)
        	vaccination.vaccination_opinion[i] = 1;
        else vaccination.vaccination_opinion[i] = 0;
    }
}

void vaccinate_elderly(Vaccination & vaccination, Population & city, Params & params,  default_random_engine& generator)
{
    uniform_real_distribution<double> distribution(0,1);

    for (int i{0}; i<city.N_use ;++i)
    if (city.Age[i]> params.age_vax_elderly_th && vaccination.vaccination_opinion[i]==0)
        if (vaccination.vaccinated[i]==VAX_NOT && distribution(generator)<params.frac_elderly_vacc/100.) {
            if (city.status[i]==0) {
            	city.status[i]= 90 ;
            	vaccination.vaccinated[i]= VAX_FUL;
            }
            else
            	vaccination.vaccinated[i]= VAX_FUL;
            vaccination.vacc_init_elderly++;
        }
}

void vaccinate_adults(Vaccination & vaccination, Population & city, Params & params, default_random_engine& generator)
{
    uniform_real_distribution<double> distribution(0,1);

    for (int i{0}; i<city.N_use ;++i) {
    	if (city.Age[i]<= params.age_vax_elderly_th && city.Age[i]> params.age_th && vaccination.vaccination_opinion[i]==0)
    		if (vaccination.vaccinated[i]==VAX_NOT && distribution(generator)<params.frac_init_vacc/100.) {
    			if (city.status[i]==0) {
    				city.status[i]= 90;
    				vaccination.vaccinated[i]= VAX_FUL;
    			}
    			else { // si pas susceptible , mis comme vacciné
    				vaccination.vaccinated[i]= VAX_FUL;
    			}
            vaccination.vacc_init_adults++;
        }
    }
}

void initVaccination(Vaccination & vaccination, Population & city, Places & places, Params &params, default_random_engine& generator,SimulationFiles & sf)
{
    uniform_real_distribution<double> distribution(0,1);

    // People vaccination opinion (1 yes, 0 no)
    vaccination.vaccination_opinion.resize(city.N_use);
    fill(vaccination.vaccination_opinion.begin(),vaccination.vaccination_opinion.end(),0);

    // Nodes (among the detected) for which step vaccination in HH was triggered
    vaccination.vaccination_triggered.resize(city.N_use);
    fill(vaccination.vaccination_triggered.begin(),vaccination.vaccination_triggered.end(),0);

    // People who have been vaccinated (0= no vacc; 1= vacc; 9= vacc waiting)
    vaccination.vaccinated.resize(city.N_use);
    fill(vaccination.vaccinated.begin(),vaccination.vaccinated.end(),VAX_NOT);


    // Reset vaccines' counting
    vaccination.vacc_total_rand=0;
    vaccination.vacc_total_react=0;
    vaccination.vacc_init_elderly=0;
    vaccination.vacc_init_adults=0;

    // sets vaccination Opinion
    setVaccinationOpinion(vaccination, city, params, generator);

	// Vaccine elderlies
	if (params.vac_elderly_inplace==1)
		vaccinate_elderly(vaccination, city, params, generator);

	if (params.vac_init_inplace==1)
			vaccinate_adults(vaccination, city, params, generator);

	vaccination.peopleToVaccinateRandom.clear(); // clear people for random vax
	vaccination.peopleToVaccinateReactive.clear(); // clear people for reactive vax


	vaccination.S_is_not_vaccinated = places.Ss;
	// fill in the places not vaccinated, filtering on size
	for (string wp : places.WPs) {
		if ((int) places.WP_nodes.size() > params.vacc_min_WP_size) {
			vaccination.WP_is_not_vaccinated.insert(wp);
		}
	}
	sf.miscellaneous_file << "initialized " << vaccination.WP_is_not_vaccinated.size() << " WPs size>"<<params.vacc_min_WP_size<<endl;

}


set<string>::iterator select_random(const set<string> &s, int n) {
  set<string>::iterator it = std::begin(s);
  // 'advance' the iterator n times
  advance(it,n);
  return it;
}


void choosePeopleToVaccinate(Vaccination & vaccination, Params & params, default_random_engine & generator) {
	// this is for random vaccination
	// shuffle id list of candidates for vaccination and insert in peopleToVaccinate
	shuffle(vaccination.allNotVaccinated.begin(), vaccination.allNotVaccinated.end(), generator);
	// now get from the end
	while (vaccination.vacc_total_rand< params.vacc_limit_rand_pop &&
			(int) vaccination.peopleToVaccinateRandom.size() < params.vac_daily_randl_pop &&
			vaccination.allNotVaccinated.size()> 0 ) {
		// to generate a random number between 0 and N_use-1 (used to select individuals randomly)
				vaccination.peopleToVaccinateRandom.push_back(vaccination.allNotVaccinated.back());
				vaccination.allNotVaccinated.pop_back();
	}
}

void doVaccination(vector<int> & peopleToVaccinate, int & total_counter, int & limit_daily, int & limit_pop, Vaccination & vaccination, Population & city) {
	// this vaccinate people waiting to be vaccinated
	// people here are not clinical / nor vaccinated as it is checked before.
	int nv=0;
	// now we vaccinate
	while (total_counter< limit_pop &&
			nv < limit_daily && peopleToVaccinate.size()>0) {
		int idxToVaccinate = peopleToVaccinate.back();
		if (city.status[idxToVaccinate]==0)
		{ // si susceptible, entre dans la vaccination
			//PYB - primoVaccination, put 8
			city.status[idxToVaccinate] = 908 ;
			vaccination.vaccinated[idxToVaccinate] = VAX_LOW;
		} else
		{ // infecté non vacciné met vacciné, les autres sont déjà vaccinés
			vaccination.vaccinated[idxToVaccinate] = VAX_FUL;
		}
		peopleToVaccinate.pop_back();

		nv++; // daily number
		total_counter++; //total number
	}

}


#define VAX_RAND_POP 1
#define VAX_RAND_WORKERS 2
#define VAX_RAND_WORKERS_SCHOOLS 3
#define VAX_WHOLE_WORKPLACE 4
#define VAX_WHOLE_SCHOOL 5

void random_vaccination(int t, Vaccination & vaccination, Population & city, Places & places, Params & params, default_random_engine& generator, SimulationFiles & sf)
{
	// fill a list with ids of peaople to vaccinate
	// in the 3 firsts only introduce so many as we can vaccinate daily
	// in the others take care of backlog of people already scheduled to keep order
	if ((int) vaccination.peopleToVaccinateRandom.size() < params.vac_daily_randl_pop) {
		// if the list can get more people today
		if (params.random_type == VAX_RAND_POP) { // random pick in the whole population
			// empty list of unvaccinated
			vaccination.allNotVaccinated.clear();
			// get a list of all unvaccinated
			for (int i=0; i < city.N_use; i++) {
				if (vaccination.vaccinated[i] == VAX_NOT && city.status[i]!=12 && city.Age[i]> params.age_th && vaccination.vaccination_opinion[i]==0) {
					if (find(vaccination.peopleToVaccinateRandom.begin(), vaccination.peopleToVaccinateRandom.end(),i) == vaccination.peopleToVaccinateRandom.end())  {
						vaccination.allNotVaccinated.push_back(i);
					}
				}
			}
			choosePeopleToVaccinate(vaccination, params, generator);
			sf.miscellaneous_file << "time "<<t<<" added for vax " << vaccination.peopleToVaccinateRandom.size() << "out of " << vaccination.allNotVaccinated.size() << " unvaccinated"<<endl;
		}  else if (params.random_type == VAX_RAND_WORKERS) { // random pick in WPs
			// empty list of unvaccinated
			vaccination.allNotVaccinated.clear();
			// Filling vector of individuals to vaccinate for random vaccination of the workers
			for(map<int,string>::iterator iter = places.node_WP.begin(); iter != places.node_WP.end(); ++iter) {
				int i= iter->first;
				if (vaccination.vaccinated[i] == VAX_NOT && city.status[i]!=12 && city.Age[i]> params.age_th && vaccination.vaccination_opinion[i]==0)
					vaccination.allNotVaccinated.push_back(i);
			}
			choosePeopleToVaccinate(vaccination, params, generator);
			sf.miscellaneous_file << "time "<<t<<" added for vax " << vaccination.peopleToVaccinateRandom.size() << "out of " << vaccination.allNotVaccinated.size() << " unvaccinated" <<endl;
		} else if (params.random_type == VAX_RAND_WORKERS_SCHOOLS) {
			// empty list of unvaccinated
			vaccination.allNotVaccinated.clear();
			// Filling vector of individuals to vaccinate for random vaccination of the workers and students
			for(map<int,string>::iterator iter = places.node_WP.begin(); iter != places.node_WP.end(); ++iter) {
				int i= iter->first;
				if (vaccination.vaccinated[i] == VAX_NOT && city.status[i]!=12 && city.Age[i]> params.age_th && vaccination.vaccination_opinion[i]==0)
					if (find(vaccination.peopleToVaccinateRandom.begin(), vaccination.peopleToVaccinateRandom.end(),i) == vaccination.peopleToVaccinateRandom.end())  {
						// not to be vaccinated already
						vaccination.allNotVaccinated.push_back(i);
					}
			}
			for(map<int,string>::iterator iter = places.node_S.begin(); iter != places.node_S.end(); ++iter) {
				int i= iter->first;
				if (vaccination.vaccinated[i] == VAX_NOT && city.status[i]!=12 && city.Age[i]> params.age_th && vaccination.vaccination_opinion[i]==0)
					if (find(vaccination.peopleToVaccinateRandom.begin(), vaccination.peopleToVaccinateRandom.end(),i) == vaccination.peopleToVaccinateRandom.end())  {
						vaccination.allNotVaccinated.push_back(i);
					}
			}
			choosePeopleToVaccinate(vaccination, params, generator);
			sf.miscellaneous_file << "time "<<t<<" added for vax " << vaccination.peopleToVaccinateRandom.size()<< "out of " << vaccination.allNotVaccinated.size() << " unvaccinated" <<endl;
		} else if (params.random_type == VAX_WHOLE_WORKPLACE){ // choose people at random and vaccinate their workplaces
			int nb_vaxx=vaccination.peopleToVaccinateRandom.size(); // get the backlog
			vector<int> backlog;
			backlog = vaccination.peopleToVaccinateRandom; // copy those that are already listed
			vaccination.peopleToVaccinateRandom.clear(); // clear
//			if (nb_vaxx >0) {sf.miscellaneous_file << "time t " <<t << "backlog " << nb_vaxx<< " for vax"<<endl;}
//			sf.miscellaneous_file << " workplaces " <<vaccination.WP_is_not_vaccinated.size() <<  " (not vaccinated) " << vaccination.WPS_is_vaccinated.size() <<  " (vaccinated) " <<endl;
			int cont=0;
			int n_WPs=0;
			while (cont==0 && nb_vaxx < params.vac_daily_randl_pop && vaccination.WP_is_not_vaccinated.size()>0) {
				uniform_int_distribution<int> distribution_int(0,vaccination.WP_is_not_vaccinated.size()-1);  // choose place at random
				set<string>::iterator it = select_random(vaccination.WP_is_not_vaccinated, distribution_int(generator)); //it  pointeur vers l'élément du set
				vaccination.WPS_is_vaccinated.insert(*it);
				vaccination.WP_is_not_vaccinated.erase(*it);
				n_WPs++;
				if (vaccination.WP_is_not_vaccinated.size()==0) cont=1; // no more WPs to vaccinate
				// put people in the vaccination list
				for (int i: places.WP_nodes[*it]) {
 					if (vaccination.vaccinated[i] == VAX_NOT && city.status[i]!=12 && city.Age[i]> params.age_th && vaccination.vaccination_opinion[i]==0 && city.isolation[t-1][i]==0) {
						if (find(vaccination.peopleToVaccinateRandom.begin(), vaccination.peopleToVaccinateRandom.end(),i) == vaccination.peopleToVaccinateRandom.end())  {
							vaccination.peopleToVaccinateRandom.push_back(i);
							nb_vaxx++;
						}
					}
				}
			}
			// now put back the original so that they are vaccinated the first
			for (int j: backlog)
				vaccination.peopleToVaccinateRandom.push_back(j);
			sf.miscellaneous_file << "time " << t <<" added for vax " << nb_vaxx << " in " << n_WPs << " WPs (+backlog ) "<<backlog.size() <<endl;
		} else if (params.random_type == VAX_WHOLE_SCHOOL) { // choose schools at random, vaccinate people in households in schools !!!
			int nb_vaxx=vaccination.peopleToVaccinateRandom.size(); // get the backlog
			vector<int> backlog;
			backlog = vaccination.peopleToVaccinateRandom; // copy those that are already listed
			vaccination.peopleToVaccinateRandom.clear(); // clear everythong
//			if (nb_vaxx >0) {sf.miscellaneous_file << "time t " <<t << "backlog " << nb_vaxx<< " for vax"<<endl;}
//			sf.miscellaneous_file << " schools " <<vaccination.S_is_not_vaccinated.size() <<  " (not vaccinated) " << vaccination.WPS_is_vaccinated.size() <<  " (vaccinated) " <<endl;
			int cont=0;
			int n_Ss=0;
			while (cont==0 && nb_vaxx < params.vac_daily_randl_pop && vaccination.S_is_not_vaccinated.size()>0) {
				uniform_int_distribution<int> distribution_int(0,vaccination.S_is_not_vaccinated.size()-1);  // choose place at random
				set<string>::iterator it = select_random(vaccination.S_is_not_vaccinated, distribution_int(generator)); //it  pointeur vers l'élément du set
				vaccination.WPS_is_vaccinated.insert(*it);
				vaccination.S_is_not_vaccinated.erase(*it);
				n_Ss++;
	//			sf.miscellaneous_file << "time " <<t << " added "<< *it << " for vaccination"<<endl;
				if (vaccination.S_is_not_vaccinated.size()==0) cont=1; //no more schools
				// put people in the vaccination list
				for (int j: places.S_nodes[*it]) {
					// we need the houshold of this person
					string HH = places.node_HH[j];
					for (int i : places.HH_nodes[HH]) { // loop on people in the HH
						if (vaccination.vaccinated[i] == VAX_NOT && city.status[i]!=12 && city.Age[i]> params.age_th && vaccination.vaccination_opinion[i]==0) {
							if (find(vaccination.peopleToVaccinateRandom.begin(), vaccination.peopleToVaccinateRandom.end(),i) == vaccination.peopleToVaccinateRandom.end()) {
								// not to be vaccinated already
								vaccination.peopleToVaccinateRandom.push_back(i);
								nb_vaxx++;
							}
						}
					}
				}
			}
			for (int j: backlog)
				vaccination.peopleToVaccinateRandom.push_back(j);
			sf.miscellaneous_file << "time " << t <<" added for vax " << nb_vaxx << " in " << n_Ss << " Ss(+backlog) " << backlog.size() <<endl;
		}
	}
}


void vaccinateInWorkingPLacesAndSchools(int real, int it, set<string> WPS_to_vaccinate, Vaccination & vaccination, Places & places, Population & city,
		Params & params, SimulationFiles & simulationFiles) {
	// this is for reactive vaccination
	// this here fills up peopleToVaccinate, but does not actually change status
	int vacc_given=0;//number daily
	// Loop on all WP/S to vaccinate
    for (string WPS: WPS_to_vaccinate) {
        int detectisol=0, previously_v=0, clinical=0, young=0, vacc=0, v_wasted=0, v_success=0, not_willing_to_v=0;
        //vacc: vaccines administered in the workplace/school
        //detectisol: people not vaccinated because detected as positive and in isolation at that time
        //previously_v: people who were already vaccinated at that time
        //clinical: people who were not vaccinated because clinical infectious
        //not_willing_to_v: people who decline vaccine
        //v_wasted: people vaccinated who were not susceptible
        //v_success: susceptible who got the vaccine
        //young: too young to be vaccinated
        // Vaccination campaign done in the S/W
        vaccination.WPS_is_vaccinated.insert(WPS);

        if (places.too_young_to_vaccinate.find(WPS) == places.too_young_to_vaccinate.end()) { // all are too young for vaccinationexit
        	// List of nodes belonging to the place
        	vector<int> vec;
        	if (places.Ss.find(WPS)!=places.Ss.end())
        		vec= places.S_nodes[WPS];
        	else if (places.WPs.find(WPS)!= places.WPs.end())
        		vec= places.WP_nodes[WPS];

        	// Loop on nodes belonging to the place
        	for (int i: vec) {
        		// I count it as a young if it is a young
        		if (city.Age[i]<= params.age_th)
        			++ young;
        		// I don't vaccinate if it is detected & isolated (I vaccinate an isolated that is not detected as a case)
        		if (city.detected[i]==1 && city.isolation[it-1][i]==1) {
        			++ detectisol;
        			if (vaccination.vaccinated[i]!=VAX_NOT)
        				++ previously_v;
        		} else {
        			// I don't vaccinate if it has been already vaccinated
        			if (vaccination.vaccinated[i]!=VAX_NOT) {
        				++ previously_v;
        			} else {
        				// I don't vaccinate if it is a clinical case
        				if (city.status[i]==12) {
        					++ clinical;
        				} else {
        					// If it is not too young
        					if (city.Age[i] > params.age_th) {
        						// if we run with increasedAcceptance , we dont care about opinion
        						if((params.doIncreasedVaccination ==1) || (vaccination.vaccination_opinion[i] == 0)) {
        							vaccination.peopleToVaccinateReactive.push_back(i);
        							// I VACCINATE
        							// If it is susceptible I change status to vaccinated waiting
        							if (city.status[i]==0) {
//        								city.status[i] = 909 ;
//        								vaccination.vaccinated[i] = 9;
        								++ v_success;
        							}
        							// If it is not susceptible I assume the vaccine does not have effect
        							else {
//        								vaccination.vaccinated[i] = 1;
        								++ v_wasted;
        							}
        							++vacc;
        							++ vacc_given; // today
//        							++ vaccination.vacc_total_react; //count in all vaccinations
        						}
        						else
                                ++ not_willing_to_v;
        					}
        				}
        			}
        		}
        	}
        	simulationFiles.outv_detailed_WPS_file << WPS << "," << vec.size() << "," << detectisol << "," << previously_v << "," << clinical << "," << young << "," << vacc << "," << v_wasted << "," << v_success << "," <<  not_willing_to_v << "," << it << "," <<  real << endl;
        } else {
        	simulationFiles.outv_detailed_WPS_file << WPS << "," << 0 << "," << 0 << "," << 0 << "," << 0 << "," << "all" << "," << 0 << "," << 0 << "," << 0 << "," <<  0 << "," << it << "," <<  real << endl;
        }
    }
 }

void vaccinateInHouseholds(int real, int it, set <string> HH_to_vaccinate, Vaccination & vaccination, Places & places, Population & city,
		Params & params, SimulationFiles & simulationFiles) {
	// this is for reactive vaccination
    int vacc_given=0;

    // Loop on all HHs to vaccinate
    int detectisol=0, previously_v=0, clinical=0, vacc=0, young=0, v_wasted=0, v_success=0, not_willing_to_v=0;
    // Loop on all HH to vaccinate
    for (string HH: HH_to_vaccinate)
    {
        // Loop on all nodes in the HH
        for (int i: places.HH_nodes[HH])
        {
            // I count it as a young if it is a young
            if (city.Age[i]<= params.age_th)
                ++ young;

            // I don't vaccinate if it is detected & isolated (I vaccinate an isolated that is not detected as a case)
            if (city.detected[i]==1 && city.isolation[it-1][i]==1) {
                ++ detectisol;
                if (vaccination.vaccinated[i]!=VAX_NOT)
                    ++ previously_v;
            }
            else
            {
                // I don't vaccinate if it has been already vaccinated
                if (vaccination.vaccinated[i]!=VAX_NOT)
                    ++ previously_v;
                else
                {
                    // I don't vaccinate if it is a clinical case
                    if (city.status[i]==12)
                        ++ clinical;
                    else
                    {
                        // If it is not too young
                        if (city.Age[i] > params.age_th)
                        {
                            // if it wants to be vaccinated
                        	if((params.doIncreasedVaccination ==1) || (vaccination.vaccination_opinion[i] == 0))
                            {
                                // I VACCINATE
    							vaccination.peopleToVaccinateReactive.push_back(i);
                                // If it is susceptible I change status to vaccinated waiting
                                if (city.status[i]==0)
                                {
//                                    city.status[i] = 909 ;
//                                    vaccination.vaccinated[i] = 9;
                                    ++ v_success;
                                }
                                // If it is not susceptible I assume the vaccine does not have effect
                                else
                                {
//                                	vaccination.vaccinated[i] = 1;
                                    ++ v_wasted;
                                }

                                ++ vacc;
                                ++ vacc_given;
//                                ++ vaccination.vacc_total_react;
                            }
                            else
                                ++ not_willing_to_v;
                        }
                    }
                }
            }
        }
    }
    simulationFiles.outv_detailed_HH_file << HH_to_vaccinate.size() << "," << detectisol << "," << previously_v << "," << clinical << "," << young << "," << vacc << "," << v_wasted << "," << v_success << "," <<  not_willing_to_v << "," << it << "," <<  real << endl;

}

// different vaccinations protocols
void step1_vaccination(int it, int real, Vaccination & vaccination, vector<int> &vaccination_triggers, Population & city, Places & places, Params & params, default_random_engine& generator, SimulationFiles &simulationFiles)
{
	//we don't take new candidates if not capacity enough

	if ((int) vaccination.peopleToVaccinateReactive.size() < params.vacc_daily_react) {

		// backlog to be put in the back again to keep order of vaccintaiton
		// it will be put back in the end
		vector<int> backlog;
		backlog=vaccination.peopleToVaccinateReactive;
		vaccination.peopleToVaccinateReactive.clear();


    	uniform_real_distribution<double> distribution(0,1);

    	// WP and S to vaccinate (list of WP/S to be vaccinated now)
    	set<string> WPS_to_vaccinate;
    	// HH to vaccinate (list of HH to be vaccinated now)
    	set<string> HH_to_vaccinate;
    	// change we fill up vaccination.peopleToVaccinate
    // skip adding more people if already over capacity
    	// there is room for reactive vaccination
    	// If the number of cases in WP/S is higher than a threshold (th_cluster_vacc)
    	// I check that the vaccination campaign was not already done there
    	// If not, with a certain delay (rate vacc_rate) I do the vaccination (I put that WP/S in the list of WP/S to vaccinate)
    	for (string wps: vaccination.place_reactive_vaccination) {// this has been checked for size
    		if (vaccination.WPS_is_vaccinated.find(wps)==vaccination.WPS_is_vaccinated.end() &&
    				distribution(generator)<params.r_V) {
    			WPS_to_vaccinate.insert(wps); // schedule for vaccination
    		}
    	}
    	// split in 2 steps
    	for (string wps: WPS_to_vaccinate) {
			vaccination.place_reactive_vaccination.erase(wps); // erase from list
    	}

    	// For each vaccination trigger (defined in the main) I put its HH in the list of HH to vaccinate
    	// (here delay is accounted for in the definition of vaccination trigger)
    	for (int i: vaccination_triggers)
    	{
    		string HH= places.node_HH[i];
    		HH_to_vaccinate.insert(HH);
    	}

    	// put people in the Queue
    	// Loop on all WP/S to vaccinate
    	vaccinateInWorkingPLacesAndSchools(real, it, WPS_to_vaccinate, vaccination,  places,  city, params, simulationFiles);
    	// vaccinate inHouseholds
    	vaccinateInHouseholds(real, it, HH_to_vaccinate, vaccination,  places,  city,  params, simulationFiles);

    	//put back people in the queue to keep order
    	for (int j: backlog) {
    		vaccination.peopleToVaccinateReactive.push_back(j);
    	}

    }
}

void step2_vaccination(int it, int real, Vaccination & vaccination, vector<int> &vaccination_triggers, Population & city, Places & places, Params & params, default_random_engine& generator, SimulationFiles &simulationFiles)
{
	if ((int) vaccination.peopleToVaccinateReactive.size() < params.vacc_daily_react) {

		// backlog to be put in the back again
		vector<int> backlog;
		backlog=vaccination.peopleToVaccinateReactive;
		vaccination.peopleToVaccinateReactive.clear();


		uniform_real_distribution<double> distribution(0,1);

		set<string> WPS_to_vaccinate;
		set<string> HH_to_vaccinate;

		//vaccinate WPs
		for (string wps: vaccination.place_reactive_vaccination)
			if (vaccination.WPS_is_vaccinated.find(wps)==vaccination.WPS_is_vaccinated.end() && distribution(generator)<params.r_V)
			{
				WPS_to_vaccinate.insert(wps);
				//step 2 vaccination : HH of employes of wps
				// include HH for vaccination,
				if (places.WPs.find(wps) != places.WPs.end()) {
					for (int i: places.WP_nodes.at(wps))
					{
						string HH= places.node_HH[i];
						HH_to_vaccinate.insert(HH);
					}
				} else if (places.Ss.find(wps) != places.Ss.end()) { // it is a school
					for (int i: places.S_nodes.at(wps))
					{
						string HH= places.node_HH[i];
						HH_to_vaccinate.insert(HH);
					}
				}
			}

    	// split in 2 steps
    	for (string wps: WPS_to_vaccinate) {
			vaccination.place_reactive_vaccination.erase(wps); // erase from list
    	}


		for (int i: vaccination_triggers)
		{
			string HH= places.node_HH[i];
			HH_to_vaccinate.insert(HH);
		}

		// Loop on all WP/S to vaccinate
		vaccinateInWorkingPLacesAndSchools(real, it, WPS_to_vaccinate, vaccination,  places,  city, params, simulationFiles);
		// vaccinate inHouseholds
		vaccinateInHouseholds(real, it, HH_to_vaccinate, vaccination,  places,  city,  params, simulationFiles);
		//put back already	listed for vaccination
		// put back in queue
		for (int j: backlog) {
			vaccination.peopleToVaccinateReactive.push_back(j);
		}
	}
}


void step3_vaccination(int it, int real, Vaccination & vaccination, vector<int> &vaccination_triggers, Population & city, Places & places, Params & params, default_random_engine& generator, SimulationFiles &simulationFiles)
{
	if ((int) vaccination.peopleToVaccinateReactive.size() < params.vacc_daily_react) {

		// backlog to be put in the back again
		vector<int> backlog;
		backlog=vaccination.peopleToVaccinateReactive;
		vaccination.peopleToVaccinateReactive.clear();

		uniform_real_distribution<double> distribution(0,1);

		set<string> WPS_to_vaccinate;
		set<string> HH_to_vaccinate;


		for (string wps: vaccination.place_reactive_vaccination) {
			if (vaccination.WPS_is_vaccinated.find(wps)==vaccination.WPS_is_vaccinated.end() && distribution(generator)<params.r_V)
			{
				WPS_to_vaccinate.insert(wps);
				//step 2 vaccination : HH of employes of wps
				// include HH for vaccination,
				if (places.WPs.find(wps) != places.WPs.end()) { // if it is a WP
					for (int i: places.WP_nodes.at(wps))
					{
						string HH= places.node_HH[i];
						HH_to_vaccinate.insert(HH);
					}
				} else if (places.Ss.find(wps) != places.Ss.end()) { // it is a school
					for (int i: places.S_nodes.at(wps))
					{
						string HH= places.node_HH[i];
						HH_to_vaccinate.insert(HH);
					}
				}
			}
		}

    	// split in 2 steps
    	for (string wps: WPS_to_vaccinate) {
			vaccination.place_reactive_vaccination.erase(wps); // erase from list
    	}

		for (int i: vaccination_triggers)
		{
			string HH= places.node_HH[i];
			HH_to_vaccinate.insert(HH);

			for (int j: places.HH_nodes[HH])
				if (places.node_WP.find(j)!=places.node_WP.end())
				{
					string wps= places.node_WP[j];
					if (vaccination.WPS_is_vaccinated.find(wps)==vaccination.WPS_is_vaccinated.end())
						WPS_to_vaccinate.insert(wps);
				}
		}

		// Loop on all WP/S to vaccinate
		vaccinateInWorkingPLacesAndSchools(real, it, WPS_to_vaccinate, vaccination,  places,  city, params, simulationFiles);
		// vaccinate inHouseholds
		vaccinateInHouseholds(real, it, HH_to_vaccinate, vaccination,  places,  city,  params, simulationFiles);
		// put back in queue
		for (int j: backlog) {
			vaccination.peopleToVaccinateReactive.push_back(j);
		}
	}

}


void step4_vaccination(int it, int real, Vaccination & vaccination, vector<int> &vaccination_triggers, Population & city, Places & places, Params & params, default_random_engine& generator, SimulationFiles &simulationFiles)

{
	if ((int) vaccination.peopleToVaccinateReactive.size() < params.vacc_daily_react) {

		// backlog to be put in the back again
		vector<int> backlog;
		backlog=vaccination.peopleToVaccinateReactive;
		vaccination.peopleToVaccinateReactive.clear();

		uniform_real_distribution<double> distribution(0,1);

		set<string> WPS_to_vaccinate;
		set<string> HH_to_vaccinate;

		for (string wps: vaccination.place_reactive_vaccination) {
			if (vaccination.WPS_is_vaccinated.find(wps)==vaccination.WPS_is_vaccinated.end() && distribution(generator)<params.r_V)
			{
				WPS_to_vaccinate.insert(wps);
				//step 2 vaccination : HH of employes of wps
				// include HH for vaccination,
				if (places.WPs.find(wps) != places.WPs.end()) { // if it is a WP
					for (int i: places.WP_nodes.at(wps))
					{
						string HH= places.node_HH[i];
						HH_to_vaccinate.insert(HH);
					}
				} else if (places.Ss.find(wps) != places.Ss.end()) { // it is a school
					for (int i: places.S_nodes.at(wps))
					{
						string HH= places.node_HH[i];
						HH_to_vaccinate.insert(HH);
					}
				}
			}
		}

    	// split in 2 steps
    	for (string wps: WPS_to_vaccinate) {
			vaccination.place_reactive_vaccination.erase(wps); // erase from list
    	}

		for (int i: vaccination_triggers)
		{
			string HH= places.node_HH[i];
			HH_to_vaccinate.insert(HH);

			for (int j: places.HH_nodes[HH])
			{
				if (places.node_WP.find(j)!=places.node_WP.end())
				{
					string wps= places.node_WP[j];
					if (vaccination.WPS_is_vaccinated.find(wps)==vaccination.WPS_is_vaccinated.end())
						WPS_to_vaccinate.insert(wps);
				}
				if (places.node_S.find(j)!=places.node_S.end())
				{
					string wps= places.node_S[j];
					if (vaccination.WPS_is_vaccinated.find(wps)==vaccination.WPS_is_vaccinated.end())
						WPS_to_vaccinate.insert(wps);
				}
			}
		}

		// Loop on all WP/S to vaccinate
		vaccinateInWorkingPLacesAndSchools(real, it, WPS_to_vaccinate, vaccination,  places,  city, params, simulationFiles);
		// vaccinate inHouseholds
		vaccinateInHouseholds(real, it, HH_to_vaccinate, vaccination,  places,  city,  params, simulationFiles);
		// put back in line
		for (int j: backlog) {
			vaccination.peopleToVaccinateReactive.push_back(j);
		}
	}
}


void vaccinationStep(int real, int it, Vaccination & vaccination, Population & city, Places &places, Compartments & epid, Params & params,
		default_random_engine & generator, SimulationFiles & simulationFiles)
{
	 
	 // VACCINATION
	 
	uniform_real_distribution<double> distribution(0,1);

    // Random vaccination of the population
	// we fill in a list of people (peopleToVaccinate) that will be vaccinated these are chosen at random
    if (params.vac_rand_inplace==1 && vaccination.vacc_total_rand < params.vacc_limit_rand_pop)
    {
        random_vaccination(it, vaccination, city, places, params, generator, simulationFiles);
        // performs vaccination
        int candidates_vax_rand = vaccination.peopleToVaccinateRandom.size();
    	doVaccination(vaccination.peopleToVaccinateRandom, vaccination.vacc_total_rand, params.vac_daily_randl_pop,params.vacc_limit_rand_pop,  vaccination, city);
        epid.vaccines_used_random[it] = candidates_vax_rand - vaccination.peopleToVaccinateRandom.size();
        simulationFiles.miscellaneous_file << "r " << real << " t " << it << " random vax " << epid.vaccines_used_random[it] << " backlog " << vaccination.peopleToVaccinateRandom.size() <<endl;
    }
    
    // Reactive vaccination
    if (params.reactive_vacc_inplace==1 && vaccination.vacc_total_react<params.vacc_limit_react)
    {
        // If a node is detected after a certain delay (measured by rate r_V)
        // it triggers the vaccination in HH (becames a vaccination trigger)
        // once it does that it is marked as vaccination_triggered=1 so it cannot trigger the vaccination anymore
        vector<int> vaccination_triggers;
        for (int i{0} ; i<city.N_use ; ++i)
        if (city.detected[i]==1)
            if ((vaccination.vaccination_triggered[i]==0) && distribution(generator)<params.r_V)
            {
                vaccination.vaccination_triggered[i]=1;
                vaccination_triggers.push_back(i);
            }
        //Step 1 vaccination
        if (params.react_type == 1) step1_vaccination( it, real, vaccination, vaccination_triggers, city, places, params, generator, simulationFiles);
        //Step 2 vaccination
        else if (params.react_type == 2) step2_vaccination( it, real, vaccination,  vaccination_triggers, city, places, params, generator, simulationFiles);
        //Step 3 vaccination
        else if (params.react_type == 3) step3_vaccination( it, real, vaccination, vaccination_triggers, city, places, params, generator, simulationFiles);
        //Step 4 vaccination
        else if (params.react_type == 4) step4_vaccination( it, real, vaccination, vaccination_triggers, city, places, params, generator, simulationFiles);

        // do vaccination and prints to miscellaneous
        int candidates_vax_react = vaccination.peopleToVaccinateReactive.size();
        doVaccination(vaccination.peopleToVaccinateReactive, vaccination.vacc_total_react, params.vacc_daily_react,params.vacc_limit_react, vaccination, city);
        epid.vaccines_used_reactive[it] = candidates_vax_react - vaccination.peopleToVaccinateReactive.size();
        simulationFiles.miscellaneous_file << "r " << real << " t " << it << " reactive vax " << epid.vaccines_used_reactive[it] << " backlog " << vaccination.peopleToVaccinateReactive.size() << endl;
    }
    
    epid.vaccines_used[it] = vaccination.vacc_total_rand + vaccination.vacc_total_react + vaccination.vacc_init_elderly+ vaccination.vacc_init_adults;
}
