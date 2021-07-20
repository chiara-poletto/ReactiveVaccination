/**
    Simulation of an ABM for COVID-19
    @file init.cpp
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

/**
 *
 * remove detection times in all places
 *
 */
void initImprovedCountInCases(Population & city) {
	city.timeDetectedInPlaces.clear(); // empty this
	city.numberDetectedInPlaces.clear();
}

/**
 *
 * put in infected people with different approaches
 *
 */
void initSeed(Population & city,Params &params, default_random_engine & generator, SimulationFiles & sf) {
	int initial_cases = params.init_incid_100000*city.N_use/100000.;
	vector<int> indexForSeeds;
	// CHOOSE INFECTED SEEDS AT RANDOM
	if (params.typeSeeding==0) {
		for (int k=0; k < city.N_use; k++) {
			if (city.status[k]==0) indexForSeeds.push_back(k);
		}
		cout << "initialized " << initial_cases << " exposed at random" <<endl;
		sf.miscellaneous_file << "initialized " << initial_cases << " exposed at random" <<endl;
		shuffle(indexForSeeds.begin(), indexForSeeds.end(),generator);
		for (int k=0; k < initial_cases; k++)
			city.status[indexForSeeds[k]]= 10;
	} else
		// CHOOSE INFECTED SEEDS FROM INFECTEDS IN IMMUNITY FILE
		if (params.typeSeeding==1) {
		sf.miscellaneous_file<< "initialized " << initial_cases << " exposed from list" <<endl;
		cout << "initialized " << initial_cases << " exposed from list" <<endl;
		shuffle(city.exposed_file[0].begin(), city.exposed_file[0].end(),generator);
		int j=0; // index of all cases
		int k=0; // number of exposed
		while (k < initial_cases) {
			int index = city.exposed_file[1][city.exposed_file[0][j]];
			int state = city.exposed_file[2][city.exposed_file[0][j]];
			if (state == 10) {
				city.status[index]= state;
				k++;
			}
			j++;
		}
	} else
		// CHOOSE INFECTED SEEDS FROM INFECTEDS IN IMMUNITY FILE AND PUT OTHER STATUS AS WELL
		if (params.typeSeeding==2) {
		shuffle(city.exposed_file[0].begin(), city.exposed_file[0].end(),generator);
		int j=0; // index of all cases
		int k=0; // number of exposed
		while (k < initial_cases) {
			int index = city.exposed_file[1][city.exposed_file[0][j]];
			int state = city.exposed_file[2][city.exposed_file[0][j]];
			city.status[index]= state;
			if (state == 10) k++;
			j++;
		}
		sf.miscellaneous_file<< "initialized " << j << " cases from file with " << k << " initial_cases" <<endl;
		cout << "initialized " << j << " cases from file with " << k << " initial_cases" <<endl;
	} else
		// CHOOSE INFECTED SEEDS FROM INFECTEDS IN IMMUNITY FILE - ALL
		if (params.typeSeeding==3) {
		int k=0;
		for (int j: city.exposed_file[0]) {
			int index = city.exposed_file[1][j];
			int state = city.exposed_file[2][j];
			city.status[index]= state;
			if (state == 10) k++;
		}
		sf.miscellaneous_file<< "initialized " << city.exposed_file[0].size() << " cases from file with " << k << " initial_cases" <<endl;
		cout <<  "initialized " << city.exposed_file[0].size() << " cases from file with " << k << " initial_cases" <<endl;
	} else
		// CHOOSE INFECTED SEEDS FROM INFECTEDS IN IMMUNITY FILE AND OTHER AS WELLFOR ONE PART, AT RANDOM FOR THE REST
		if (params.typeSeeding==4) {
		shuffle(city.exposed_file[0].begin(), city.exposed_file[0].end(),generator);
		int j=0; // index of all cases
		int k=0; // number of exposed
		vector<int> seedStatus; // to hold the status of people
		while (k < initial_cases) {
			int state = city.exposed_file[2][city.exposed_file[0][j]];
			seedStatus.push_back(state);
			if (state == 10) k++;
			j++;
		}
		// then get people and allocate at random
		for (int k=0; k < city.N_use; k++) {
			if (city.status[k]==0) indexForSeeds.push_back(k);
		}
		shuffle(indexForSeeds.begin(), indexForSeeds.end(),generator);
		for (int k=0; k < (int) seedStatus.size(); k++)
			city.status[indexForSeeds[k]]= seedStatus[k];
		sf.miscellaneous_file<< "obtained " << j << " cases from file with " << k << " exposed cases" <<endl;
		cout << "initialized " << seedStatus.size() << " cases from file with " << initial_cases << " initial_cases and randomized in "<< indexForSeeds.size() << " people" <<endl;
	}  else
		//CHOOSE INFECTED FOR INVASION
		if (params.typeSeeding==5) {
			// we seed with at least R so that one starts an epidemic...
			// params.beta */(1/params.e_rate + 1/params.p1_rate + 1/params.gamm)
			//int nSeed = (int) (params.beta /(1/params.e_rate + 1/params.p1_rate + 1/params.gamm))+1;
			int nSeed=3;
			//we seed only one person at random
			for (int k=0; k < city.N_use; k++) {
				if (city.status[k]==0) indexForSeeds.push_back(k);
			}
			shuffle(indexForSeeds.begin(), indexForSeeds.end(),generator);
			for (int i=0; i < nSeed; i++) {
				city.status[indexForSeeds[i]]= 10;
			}
			cout << "initialized "<<nSeed<< "case at time 0" <<endl;
		}
}

/**
 *
 * intialize population status. taken at random, from file or bit of both
 *
 */
void initImmunity(Population & city, Params &params,default_random_engine& generator, SimulationFiles & sf) {
	 uniform_real_distribution<double> distribution(0,1);

	 if (params.immunity_option!=0) {
		 // If we want some immunity
		 if (params.immunity_random==0) {
			 // Initialisation of immunity from external file
			 city.status = city.status_file; //this is already 0/2
			 int initial_imm= counter(city.status_file,2);
			 double p_imm= initial_imm*1.0/city.N_use;
			 sf.miscellaneous_file << " immunity from file at level" << p_imm*100 <<"%"<<endl;
			 cerr << " immunity from file at level" << p_imm*100 <<"%"<<endl;
		 } else if (params.immunity_random==1 ) {
			 int initial_imm= counter(city.status_file,2);
			 double p_imm= initial_imm*1.0/city.N_use;
			 for (int i{0}; i<city.N_use ;++i)
				 if (distribution(generator)<p_imm)
					 city.status[i]=2;
			 sf.miscellaneous_file << " immunity random at level" << p_imm*100 <<"%"<<endl;
			 cerr << " immunity random at level" << p_imm*100 <<"%"<<endl;
		 } else if (params.immunity_random == 2) {
			 // initialize part from random, part from file
			 vector<int> idxImmuneIndividuals; // indices of the immune
			 vector<int> idxNonImmuneIndividuals; // indices of the non immune
			 int nImmune=0;
			 for (int i=0; i< (int) city.status_file.size(); i++) {
				 if (city.status_file[i] == 0) {
					 idxNonImmuneIndividuals.push_back(i);
				 } else {
					 idxImmuneIndividuals.push_back(i);
					 nImmune++; // immune in initial file
				 }
			 }
			 // random shuffle indices
			 shuffle(idxNonImmuneIndividuals.begin(), idxNonImmuneIndividuals.end(),generator);
			 shuffle(idxImmuneIndividuals.begin(), idxImmuneIndividuals.end(),generator);
			 // fill in :
			 // random part
			 for (int i=0; i<(int) (nImmune * params.pct_rand_init_immunity);i++) city.status[idxNonImmuneIndividuals[i]]=2;
			 // from file
			 for (int i=(int) (nImmune*params.pct_rand_init_immunity); i<nImmune;i++) city.status[idxImmuneIndividuals[i]]=2;
		 }
	 }
	 else
	 {
		 cerr << "Immunity 0 "<<endl;
		 sf.miscellaneous_file << "Immunity 0 "<<endl;
		 for (int i{0}; i<city.N_use ;++i)
				 city.status[i]=0;
	 }
}

/**
 *
 * reset all compartments
 *
 */
void initEpid(Compartments &epid, Params &params) {
	
    //allocate space & fill with  0s
    epid.S.resize(params.max_steps);
    fill(epid.S.begin(), epid.S.end(),0);
    epid.E.resize(params.max_steps);
    fill(epid.E.begin(), epid.E.end(),0);
    epid.P1.resize(params.max_steps);
    fill(epid.P1.begin(), epid.P1.end(),0);
    epid.P2.resize(params.max_steps);
    fill(epid.P2.begin(), epid.P2.end(),0);
    epid.SI.resize(params.max_steps);
    fill(epid.SI.begin(), epid.SI.end(),0);
    epid.AI.resize(params.max_steps);
    fill(epid.AI.begin(), epid.AI.end(),0);
    epid.R.resize(params.max_steps);
    fill(epid.R.begin(), epid.R.end(),0);
    epid.wvS.resize(params.max_steps);
    fill(epid.wvS.begin(), epid.wvS.end(),0);
    epid.vS.resize(params.max_steps);
    fill(epid.vS.begin(), epid.vS.end(),0);
    epid.vE.resize(params.max_steps);
    fill(epid.vE.begin(), epid.vE.end(),0);
    epid.vP1.resize(params.max_steps);
    fill(epid.vP1.begin(), epid.vP1.end(),0);
    epid.vP2.resize(params.max_steps);
    fill(epid.vP2.begin(), epid.vP2.end(),0);
    epid.vSI.resize(params.max_steps);
    fill(epid.vSI.begin(), epid.vSI.end(),0);
    epid.vAI.resize(params.max_steps);
    fill(epid.vAI.begin(), epid.vAI.end(),0);
    epid.vR.resize(params.max_steps);
    fill(epid.vR.begin(), epid.vR.end(),0);
    
    epid.infected.resize(params.max_steps);
    fill(epid.infected.begin(), epid.infected.end(),0);
    epid.WP_infected.resize(params.max_steps);
    fill(epid.WP_infected.begin(), epid.WP_infected.end(),0);
    epid.S_infected.resize(params.max_steps);
    fill(epid.S_infected.begin(), epid.S_infected.end(),0);

    epid.isolated.resize(params.max_steps);
    fill(epid.isolated.begin(), epid.isolated.end(),0);
	epid.infected_and_isolated.resize(params.max_steps);
    fill(epid.infected_and_isolated.begin(), epid.infected_and_isolated.end(),0);
	epid.vaccines_used.resize(params.max_steps);
    fill(epid.vaccines_used.begin(), epid.vaccines_used.end(),0);
	epid.vaccines_used_random.resize(params.max_steps);
    fill(epid.vaccines_used_random.begin(), epid.vaccines_used_random.end(),0);
	epid.vaccines_used_reactive.resize(params.max_steps);
    fill(epid.vaccines_used_reactive.begin(), epid.vaccines_used_reactive.end(),0);
    // vaccinated
	epid.vaccinated_FUL.resize(params.max_steps);
    fill(epid.vaccinated_FUL.begin(), epid.vaccinated_FUL.end(),0);
	epid.vaccinated_MID.resize(params.max_steps);
    fill(epid.vaccinated_MID.begin(), epid.vaccinated_MID.end(),0);
	epid.vaccinated_LOW.resize(params.max_steps);
    fill(epid.vaccinated_LOW.begin(), epid.vaccinated_LOW.end(),0);
    //number detected
    epid.detected.resize(params.max_steps);
    fill(epid.detected.begin(), epid.detected.end(),0);
    //number found by CT
    epid.found_by_ct.resize(params.max_steps);
    fill(epid.found_by_ct.begin(), epid.found_by_ct.end(),0);
    //number found by screening
    epid.found_by_regular_screening.resize(params.max_steps);
    fill(epid.found_by_regular_screening.begin(), epid.found_by_regular_screening.end(),0);
    //number found by screening
    epid.found_by_reactive_screening.resize(params.max_steps);
    fill(epid.found_by_reactive_screening.begin(), epid.found_by_reactive_screening.end(),0);
}

/**
 * reset incidence
 */
void initIncidence(Incidences & incidence, Params &params) {
    incidence.P1.resize(params.N_real,vector<int>(params.max_steps));
    incidence.SI.resize(params.N_real,vector<int>(params.max_steps));
    incidence.P2.resize(params.N_real,vector<int>(params.max_steps));
    incidence.AI.resize(params.N_real,vector<int>(params.max_steps));
    incidence.vP1.resize(params.N_real,vector<int>(params.max_steps));
    incidence.vSI.resize(params.N_real,vector<int>(params.max_steps));
    incidence.vP2.resize(params.N_real,vector<int>(params.max_steps));
    incidence.vAI.resize(params.N_real,vector<int>(params.max_steps));
}

/**
 * reset realisation structuresto0
 *
 */
void initCity(Population & city,Params & params)
{
    // Nstatus set to susceptibles
    city.status.resize(city.N_use);
    fill(city.status.begin(),city.status.end(), 0);
    // Nodes that will be detected after a delay will be passed
    city.will_be_detect.resize(city.N_use);
    fill(city.will_be_detect.begin(),city.will_be_detect.end(), 0);
    // Nodes detected, compliant HH & MCO contacts are also detected
    city.detected.resize(city.N_use);
    fill(city.detected.begin(),city.detected.end(),0);
    // Isolation status of nodes in time (1= isolated; 0= not isolated) first iteration, then individual
    city.isolation.resize(params.max_steps, vector<int> (city.N_use));
    fill(city.isolation[0].begin(), city.isolation[0].end(), 0);
    // Time of symptom onset for clinical cases
    city.onset_time.resize(city.N_use);
    fill(city.onset_time.begin(),city.onset_time.end(),0);
    // isolated
    city.iso_time.resize(city.N_use);
    fill(city.iso_time.begin(),city.iso_time.end(),0);
    // Time when (last) isolation begane
    city.iso_time.resize(city.N_use);
    fill(city.iso_time.begin(),city.iso_time.end(),0);
    

}

/**
 *
 *  	 erase networks, set new links from existing possibilities
 *
 *
 */
void initNetworks(Networks & contacts, int N_use, Places & places, Params & params, default_random_engine& generator) {
	for (int layer=0; layer<5;layer++) {
		contacts.keptContactsInLayers[layer]=0;
		for (int m=0; m <params.number_networks; m++) {
			for (int i=0; i<N_use;i++) {
					(contacts.layers[layer][m][i]).clear();
			}
		}
    }
    uniform_real_distribution<double> distribution(0,1);

	// find people who are working
	vector<int> isTeleWorking;
	isTeleWorking.resize(N_use);
	fill(isTeleWorking.begin(),isTeleWorking.end(),0);
// fill with 1, then cancel after
	for (string idPlace : 	places.WPs) {// loop on WPs
		for (int idWorker : places.WP_nodes[idPlace]) { //loop on people inside
			if (distribution(generator) < params.pct_teleworking) {
				isTeleWorking[idWorker]=1;
			}
		}
	}


    // Loop through all contacts
	for (int ind1{0}; ind1< (int) contacts.node1.size();++ind1)
	{
		// Create a random vector of numbers for each network I want to generate
		vector<double> rand_vec_days = random_vector(generator, params.number_networks);

		for (int m{0};m<params.number_networks;++m)
		{
            // freq is the activation rate, so here I decide if the contact is active
            if (rand_vec_days[m]<contacts.freq[ind1])
            {
            	// place it in the corresponding layer :
            	// - always if household /school
            	// - if LAYER IS WORK/TRAVEL only if none are teleworking
            	// - if LAYER is COMMUNITY with probability params.pct_comty_contact
            	if ((contacts.place[ind1] == LAYER_HH) ||
            			(contacts.place[ind1] == LAYER_S) ||
            			((contacts.place[ind1] == LAYER_WP) && (isTeleWorking[contacts.node1[ind1]]==0) && (isTeleWorking[contacts.node2[ind1]]==0)) ||
            			((contacts.place[ind1] == LAYER_TR) && (isTeleWorking[contacts.node1[ind1]]==0) && (isTeleWorking[contacts.node2[ind1]]==0)) ||
						((contacts.place[ind1] == LAYER_CO) && (distribution(generator) < params.pct_comty_contact))
            			) {
            		(contacts.layers[contacts.place[ind1]][m][contacts.node1[ind1]]).push_back(contacts.node2[ind1]);
            		(contacts.layers[contacts.place[ind1]][m][contacts.node2[ind1]]).push_back(contacts.node1[ind1]);
            		contacts.keptContactsInLayers[contacts.place[ind1]]++;
            	}
            }
        }
    }

    // reset the network number
	contacts.idxNetwork = 0;
	// erase networks
    cout<<"kept Contacts : HH " << contacts.keptContactsInLayers[0]/params.number_networks << "(" << contacts.keptContactsInLayers[0]*100/contacts.allContactsInLayers[0]/params.number_networks <<"%) " <<
    		"WP " << contacts.keptContactsInLayers[1]/params.number_networks << "(" << contacts.keptContactsInLayers[1]*100/contacts.allContactsInLayers[1]/params.number_networks <<"%) " <<
			"S " << contacts.keptContactsInLayers[2]/params.number_networks << "(" << contacts.keptContactsInLayers[2]*100/contacts.allContactsInLayers[2]/params.number_networks <<"%) " <<
			"CO " << contacts.keptContactsInLayers[3]/params.number_networks << "(" << contacts.keptContactsInLayers[3]*100/contacts.allContactsInLayers[3]/params.number_networks <<"%) " <<
			"TR " << contacts.keptContactsInLayers[4]/params.number_networks << "(" << contacts.keptContactsInLayers[4]*100/contacts.allContactsInLayers[4]/params.number_networks <<"%) " <<endl;
}

/**
 *
 *  init Places
 * //for all places, define  whether the place is too young to test; to young to vaccinate based on params
 *
 *
 */
void initPlaces(Places & places, Population & city, Params & params) {
	city.has_WP.resize(city.N_use);
    fill(city.has_WP.begin(), city.has_WP.end(),0);
	for (string idPlace : 	places.WPs) {// loop on WPs
		int too_young_to_test=1; // 1 if true
		int too_young_to_vaccinate=1;
		for (int i : places.WP_nodes[idPlace]) { //loop on people inside
			if (city.Age[i] > params.age_min_testing) too_young_to_test=0;
			if (city.Age[i] > params.age_th) too_young_to_vaccinate=0;
			city.has_WP[i]=1;
		}
		if (too_young_to_test == 1)  {places.too_young_to_test.insert(idPlace);}
		if ((params.regular_screening_WP_inplace == 1) && ((int) places.WP_nodes[idPlace].size() < params.screening_min_size))
		{places.too_small_to_test.insert(idPlace);}
		if (too_young_to_vaccinate == 1) places.too_young_to_vaccinate.insert(idPlace);
	}
	city.has_school.resize(city.N_use);
    fill(city.has_school.begin(), city.has_school.end(),0);
	for (string idPlace : 	places.Ss) {// loop on WPs
		int too_young_to_test=1; // 1 if true
		int too_young_to_vaccinate=1;
		for (int i : places.S_nodes[idPlace]) { //loop on people inside
			if (city.Age[i] > params.age_min_testing) too_young_to_test=0;
			if (city.Age[i] > params.age_th) too_young_to_vaccinate=0;
			city.has_school[i]=1;
		}
		if (too_young_to_test == 1) places.too_young_to_test.insert(idPlace);
		if ((params.regular_screening_S_inplace	== 1) && ((int) places.S_nodes[idPlace].size() < params.screening_min_size))
		{places.too_small_to_test.insert(idPlace);}
		if (too_young_to_vaccinate == 1) places.too_young_to_vaccinate.insert(idPlace);
	}
}

/**
 *
 * allocate memory to store all vectors. improve efficiency.
 * allocate acquainatances (used for contact tracing)
 *
 */
void initDynamicalNetworks(Networks& contacts, int N_use,Params & params) {
	 // first read Network, then compute dynamical networks
	// The network is a 3-dim vector
	// The 3 dimensions have size "number_networks", N_use and variable
	// We have "number_networks" networks, each storing the contacts for the N_use individuals that accurr during a day
	// Networks are generated extracting contacts with probability freq[ind1]
	// we add outer dimension to hold all networks in same structure
	// 0 :Household, 1: school, 2 : workplace, 3 : transport , 4 : other
	// form outermost to innermost : LAYER/NETWORKS/N_USE/INDIV
	// also set acquaintances
// set weights_layers
	contacts.weights_layer[0] = params.wgt_lyr_household;
	contacts.weights_layer[1] = params.wgt_lyr_workspace;
	contacts.weights_layer[2] = params.wgt_lyr_school;
	contacts.weights_layer[3] = params.wgt_lyr_community;
	contacts.weights_layer[4] = params.wgt_lyr_transport;

	// resize for 5 layers
	// this will be put to size = 0 later
	contacts.layers.resize(5, vector<vector<vector<int>>>(params.number_networks,vector<vector<int>>(N_use)));

	// pre-allocate memory
	for (int layer=0; layer<5;layer++) {
		for (int m=0; m <params.number_networks; m++) {
			for (int i=0; i<N_use;i++) {
				if (contacts.maxContacts[i][layer]>0) {
					(contacts.layers[layer][m][i]).reserve(contacts.maxContacts[i][layer]);
				}
			}
		}
    }
    // Define acquaintances (i.e. frequent contacts) for manual contact tracing
    //resize for N_use
    contacts.acquaintances.resize(N_use);
    for (int ind1{0}; ind1< (int) contacts.node1.size();++ind1)
    	if (contacts.freq[ind1]>params.acquai_frequency) {
    		contacts.acquaintances[contacts.node1[ind1]].insert(contacts.node2[ind1]);
    		contacts.acquaintances[contacts.node2[ind1]].insert(contacts.node1[ind1]);
    	}
}

/**
 * convert params from list to struc to avoid runtime dereferencing
 * params grouped by topic
 */
void convertParams(Params & params, MapParams & mapParams) {
	// REMEMBER : UPON CHANGE ALSO ADD LINE IN printParams (InputOutput.cpp)
	//simulation
	params.max_steps=mapParams["max_steps"];
	params.N_real=mapParams["N_real"];
	params.number_networks=mapParams["number_networks"];
	params.refresh_network=mapParams["refresh_network"];

	// print
	params.print_infectors=mapParams["print_infectors"];
	params.print_isolations=mapParams["print_isolations"];

	// initial
	params.immunity_option=mapParams["immunity_option"];
	params.immunity_random=mapParams["immunity_random"];
	params.init_incid_100000=mapParams["init_incid_100000"];
	params.pct_rand_init_immunity = mapParams["pct_rand_init_immunity"];

	// transmission
	params.beta=mapParams["beta"];
	params.e_rate=mapParams["e_rate"];
	params.E_ve_rate=mapParams["E_ve_rate"];
	params.factor_isol=mapParams["factor_isol"];
	params.gamm=mapParams["gamm"];
	params.p1_rate=mapParams["p1_rate"];
	params.p2_rate=mapParams["p2_rate"];
	params.r_a=mapParams["r_a"];
	params.r_V=mapParams["r_V"]; // rate
	params.v_sus_MID=mapParams["v_sus_MID"]; // reduction P(Inf|Vax)/P(Inf|nVax) in primoVax
	params.v_sus_FUL=mapParams["v_sus_FUL"]; // reduction P(Inf|Vax)/P(Inf|nVax) in fullyVax
	params.v_symp=mapParams["v_symp"]; //
	params.v_trans=mapParams["v_trans"]; //
	params.ve_rate=mapParams["ve_rate"]; //
	params.vgamm=mapParams["vgamm"]; //
	params.vp1_rate=mapParams["vp1_rate"]; //
	params.vp2_rate=mapParams["vp2_rate"]; //
	params.w_rate_MID_to_FUL=mapParams["w_rate_MID_to_FUL"]; // passage de no Protection a Mid protection
	params.w_rate_LOW_to_MID=mapParams["w_rate_LOW_to_MID"]; // passage de mid protection a ful protection

	//contacts
	params.acquai_frequency=mapParams["acquai_frequency"];
	params.interv_cumul_cases=mapParams["interv_cumul_cases"];
	params.pct_teleworking = mapParams["pct_teleworking"];
	params.pct_comty_contact = mapParams["pct_comty_contact"];


	//vaccination
	params.coverage_vacc_higher_65=mapParams["coverage_vacc_higher_65"];
	params.coverage_vacc_lower_65=mapParams["coverage_vacc_lower_65"];
	params.age_th=mapParams["age_th"];
	params.frac_elderly_vacc=mapParams["frac_elderly_vacc"];
	params.frac_init_vacc=mapParams["frac_init_vacc"];
	params.reactive_vacc_inplace=mapParams["reactive_vacc_inplace"];
	params.th_cluster_vacc=mapParams["th_cluster_vacc"];
	params.vac_daily_randl_pop=mapParams["vac_daily_randl_pop"];
	params.vac_elderly_inplace=mapParams["vac_elderly_inplace"];
	params.vac_init_inplace=mapParams["vac_init_inplace"];
	params.vac_rand_inplace=mapParams["vac_rand_inplace"];
	params.vacc_daily_react=mapParams["vacc_daily_react"];
	params.vacc_limit_rand_pop=mapParams["vacc_limit_rand_pop"];
	params.vacc_limit_react=mapParams["vacc_limit_react"];
	params.react_type=mapParams["react_type"];
    params.random_type=mapParams["random_type"];
    params.vacc_min_WP_size = mapParams["vacc_min_WP_size"];
    params.age_vax_elderly_th = mapParams["age_vax_elderly_th"];

	//screening
	params.regular_screening_HH_inplace=mapParams["regular_screening_HH_inplace"];
	params.regular_screening_S_inplace=mapParams["regular_screening_S_inplace"];
	params.regular_screening_WP_inplace=mapParams["regular_screening_WP_inplace"];
	params.age_min_testing=mapParams["age_min_testing"]; // minimum age for testing
	params.f_S_screening=mapParams["f_S_screening"]; // frequency of schools participating
	params.f_HH_screening=mapParams["f_HH_screening"]; // frequency of HHs participating
	params.f_WP_screening=mapParams["f_WP_screening"];// frequency of WPs participating
	params.testing_proportion_S=mapParams["testing_proportion_S"]; // individual participation
	params.testing_proportion_WP=mapParams["testing_proportion_WP"]; // individual participation
	params.testing_proportion_HH=mapParams["testing_proportion_HH"]; // individual participation
	params.interval_screening=mapParams["interval_screening"]; // interval screening
	params.reactive_screening_S_inplace=mapParams["reactive_screening_S_inplace"]; //0 or 1
	params.reactive_screening_WP_inplace=mapParams["reactive_screening_WP_inplace"]; //0 or 1
	params.reactive_tests_daily=mapParams["reactive_tests_daily"]; // max number of reactive  tests daily
	params.reactive_tests_limit=mapParams["reactive_tests_limit"]; // max number of reactive tests overall
	params.th_cluster_screen=mapParams["th_cluster_screen"]; // taille cluster for screening
	params.r_react_screen=mapParams["r_react_screen"]; // rate for reactive screening
	params.screening_min_size=mapParams["screening_min_size"]; // rate for reactive screening


	// contact tracing
	params.manualCT_inplace=mapParams["manualCT_inplace"];
	params.p_c_HH=mapParams["p_c_HH"];
	params.p_rA=mapParams["p_rA"];
	params.p_rS=mapParams["p_rS"];
	params.r_d_asymp=mapParams["r_d_asymp"];
	params.r_d_symp=mapParams["r_d_symp"];
	params.r_mCT=mapParams["r_mCT"];
	params.tracing_lookback=mapParams["tracing_lookback"];
	params.will_be_d_asymp=mapParams["will_be_d_asymp"];
	params.will_be_d_symp=mapParams["will_be_d_symp"];

	//after Invasion
	params.p_c_HH_aft_inv=mapParams["p_c_HH_aft_inv"];
	params.p_rA_aft_inv=mapParams["p_rA_aft_inv"];
	params.p_rS_aft_inv=mapParams["p_rS_aft_inv"];
	params.r_d_asymp_aft_inv=mapParams["r_d_asymp_aft_inv"];
	params.r_d_symp_aft_inv=mapParams["r_d_symp_aft_inv"];
	params.will_be_d_asymp_aft_inv=mapParams["will_be_d_asymp_aft_inv"];
	params.will_be_d_symp_aft_inv=mapParams["will_be_d_symp_aft_inv"];
	params.skip_iso_aft_inv = mapParams["skip_iso_aft_inv_aft_inv"];

	//isolation
	params.HHisolation_inplace=mapParams["HHisolation_inplace"];
	params.skip_iso=mapParams["skip_iso"]; // rate of getting out of isolation
	params.dur_isol_sus=mapParams["dur_isol_sus"];
	params.dur_isol_inf=mapParams["dur_isol_inf"];
	params.dur_isol_rec=mapParams["dur_isol_rec"];

	//weight layers
	params.wgt_lyr_household=mapParams["wgt_lyr_household"];
	params.wgt_lyr_workspace=mapParams["wgt_lyr_workspace"];
	params.wgt_lyr_school=mapParams["wgt_lyr_school"];
	params.wgt_lyr_community=mapParams["wgt_lyr_community"];
	params.wgt_lyr_transport=mapParams["wgt_lyr_transport"];

		// invasion
	params.n_start_reactive=mapParams["n_start_reactive"];
	params.incrVaccAccept_inplace=mapParams["incrVaccAccept_inplace"];
}


/**
 * EOF
 */
