/**
    Simulation of an ABM for COVID-19
    @file main.cpp
    @author Pierre-Yves Boelle, Chiara Poletto
    @acknowledgment Jésus A Moreno López for first version fo the code
    @version 1.0 2020-04-16
    @license GPL-3.0-or-later
*/

#define DEBUG

#include <stdio.h>
#include <stdlib.h>
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
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>
#include <utility>
#include <unistd.h>
#include<bits/stdc++.h>

#include "defs.h"
#include "routine_miscellanea.hpp"
#include "inputOutput.h"
#include "init.h"
#include "transmission.h"
#include "vaccination.h"
#include "testing.h"
#include "screening.h"


using namespace std;


/**
*	update the list of detection date in each place, then
*	count the number of cases in places over a time period given by params.interv_cumul_cases
*	@param it : time in simulation
*	@param upd : structure with list of transitions
*	@param screening :
*	@param vaccination
*	@param city :
*/
void countCasesInWorkingPlacesAndSchools(int it, Compartments & upd, Screening & screening, Vaccination & vaccination,
		Population & city, Places & places, Params & params) {
	for (int i: upd.detected) {
		// detec is index of case that is detected now
		// if case was not isolated before
		if (city.isolation[it-1][i]==0) {
			if (places.node_WP.find(i)!=places.node_WP.end()) {//idividuals has a WP
				string wps= places.node_WP[i]; // get the place
				if (city.timeDetectedInPlaces.find(wps) == city.timeDetectedInPlaces.end()) { // not in the map
					city.timeDetectedInPlaces[wps]=vector<int>(0); // add to the map
				}
				city.timeDetectedInPlaces[wps].push_back(it); // push back the time of detection
			} else if (places.node_S.find(i)!=places.node_S.end()) {//individuals has a School
				string wps= places.node_S[i]; // get the place
				if (city.timeDetectedInPlaces.find(wps) == city.timeDetectedInPlaces.end()) { // not in the map
					city.timeDetectedInPlaces[wps]=vector<int>(0); // add to the map
				}
				city.timeDetectedInPlaces[wps].push_back(it); // push back the time of detection
			}
		}
	}
	// remove date of cases that are too distant in time
	int last_valid_time= it - params.interv_cumul_cases;
	for (pair<string, vector<int>> element : city.timeDetectedInPlaces) {
		if (element.second.size()>0) { //there were cases in this place
			// remove times if there are observations and these were from before
			while ((element.second.size()>0) && (element.second[0] < last_valid_time)) {
				element.second.erase(element.second.begin()); // efficient in small vector
			}
			//compute number of cases in place
			city.timeDetectedInPlaces[element.first]=element.second;
			city.numberDetectedInPlaces[element.first]=element.second.size();
			if (city.numberDetectedInPlaces[element.first] > params.th_cluster_screen) {// there are more cases than threshold,
				if (places.too_young_to_test.find(element.first) == places.too_young_to_test.end()) { // and not too young to test
					screening.place_reactive_screening.insert(element.first); //insert
				}
			}
			if (city.numberDetectedInPlaces[element.first] > params.th_cluster_vacc) {// there are more cases than threshold,
				if (places.too_young_to_vaccinate.find(element.first) == places.too_young_to_vaccinate.end()) { // and not too young to vaccinate
					if (vaccination.WPS_is_vaccinated.find(element.first) == vaccination.WPS_is_vaccinated.end()) { // and not already vaccinated
						vaccination.place_reactive_vaccination.insert(element.first); //insert
					}
				}
			}
		}
	}
}

/**
 *
 * compute number of persons with given status
 * then compute number of persons in different boxes
 *
 */
void computePrevalence(int t, Compartments & epid, Population & city, Vaccination & vaccination, Testing & testing) {
	// take status and count
	epid.S[t] = counter(city.status,0);
	epid.E[t] = counter(city.status,10);
	epid.P1[t] = counter(city.status,101);
	epid.SI[t] = counter(city.status,12);
	epid.P2[t] = counter(city.status,102);
	epid.AI[t] = counter(city.status,13);
	epid.R[t] = counter(city.status,2);
	epid.wvS[t] = counter(city.status,909) + counter(city.status,908);
	epid.vS[t] = counter(city.status,90);
	epid.vE[t] = counter(city.status,910);
	epid.vP1[t] = counter(city.status,9101);
	epid.vSI[t] = counter(city.status,912);
	epid.vP2[t] = counter(city.status,9102);
	epid.vAI[t] = counter(city.status,913);
	epid.vR[t] = counter(city.status,92);
	epid.isolated[t] = counter(city.isolation[t],1);

	vector<int> tmp;
	tmp.resize(city.N_use);
	// compute infected in WPs
	for (int i =0; i<city.N_use; i++) {
		tmp[i] = (city.status[i]!=0 && city.status[i]!=90 && city.status[i]!=908 && city.status[i]!=909 && city.status[i]!=2 && city.status[i]!=92)*city.has_WP[i];
	}
	epid.WP_infected[t]=counter(tmp,1);
	// compute infected in schools
	for (int i =0; i<city.N_use; i++) {
		tmp[i] = (city.status[i]!=0 && city.status[i]!=90 && city.status[i]!=908 && city.status[i]!=909 && city.status[i]!=2 && city.status[i]!=92)*city.has_school[i];
	}
	epid.S_infected[t]=counter(tmp,1);
	// compute isolated and infected
	for (int i =0; i<city.N_use; i++) {
		tmp[i] = (city.status[i]!=0 && city.status[i]!=90 && city.status[i]!=908 && city.status[i]!=909 && city.status[i]!=2 && city.status[i]!=92 )*city.isolation[t][i];
	}
	epid.infected_and_isolated[t] = counter(tmp,1);

	// 0 is susceptible, 10 is exposed, 101 pre-symptomatic(1), 12 symptomatic,
	    //102 pre-symptomatic(2), 13 asymptomatic, 2 is recovered.
	    // Waiting for vaccine to act: 908
	    // primoVaccinated : 909
	    // Vaccination status numbers:
	    // 90 is susceptible, 910 is exposed,  9101 pre-symptomatic(1), 912 symptomatic,
	    // 9102 pre-symptomatic(2),  913 asymptomatic, 92 is recovered.
	for (int i =0; i<city.N_use; i++) {
		tmp[i] = (city.status[i]!=0 && city.status[i]!=90 && city.status[i]!=908 && city.status[i]!=909 && city.status[i]!=2 && city.status[i]!=92);
	}
	epid.infected[t] = counter(tmp,1);

	// compute number detected
	epid.detected[t]=counter(city.detected,1);
	// compute number vaccinated
	epid.vaccinated_FUL[t] = counter(vaccination.vaccinated,VAX_FUL);
	epid.vaccinated_MID[t] = counter(vaccination.vaccinated,VAX_MID);
	epid.vaccinated_LOW[t] = counter(vaccination.vaccinated,VAX_LOW);
	// compute number contacttracing_done
	epid.found_by_ct[t] = counter(testing.found_in[TRACING],1);
	// compute number found by regular screening
	epid.found_by_regular_screening[t] = counter(testing.found_in[REGULAR_SCREENING],1);
	// compute number found by reactive screening
	epid.found_by_reactive_screening[t] = counter(testing.found_in[REACTIVE_SCREENING],1);
}


/**
 *
 * Compute incidence as number in update
 *
 */
void computeIncidence(int real, int it, Incidences & incidence, Compartments & upd) {
    //update computed above in transmission
     incidence.P1[real][it] = upd.P1.size();
     incidence.SI[real][it] = upd.SI.size();
     incidence.P2[real][it] = upd.P2.size();
     incidence.AI[real][it] = upd.AI.size();
     incidence.vP1[real][it] = upd.vP1.size();
     incidence.vSI[real][it] = upd.vSI.size();
     incidence.vP2[real][it] = upd.vP2.size();
     incidence.vAI[real][it] = upd.vAI.size();
}

/**
 *
 * actually update transmission after transitions have been decided
 *
 */
void updateCompartments(int real, int it, Compartments & upd, Population & city, Vaccination & vaccination,
		SimulationFiles & simulationFiles) {

	for (int uppvS: upd.vS)
	 {
		 if (city.status[uppvS]==908) {
			 vaccination.vaccinated[uppvS]=VAX_MID;
			 city.status[uppvS]=909;
		 } else if (city.status[uppvS]==909) {
			 vaccination.vaccinated[uppvS]=VAX_FUL;
			 city.status[uppvS]=90;
		 }
	 }
	 for (int upe: upd.E)
		 city.status[upe]=10;
	 for (int upe: upd.vE)
		 city.status[upe]=910;
	 for (int upp1: upd.P1)
		 city.status[upp1]=101;
	 for (int upp2: upd.P2)
		 city.status[upp2]=102;
	 for (int upsi: upd.SI) {
		 city.status[upsi]=12;
		 city.onset_time[upsi]=it;
		 simulationFiles.symptom_list_file << upsi << "," << city.Age[upsi]  << ","<< vaccination.vaccinated[upsi] << "," << it << "," << real << endl;
	 }
	 for (int upai: upd.AI) {
		 city.status[upai]=13;
		 simulationFiles.asymptom_list_file << upai << "," << city.Age[upai]  << "," <<vaccination.vaccinated[upai] << "," << it << "," << real << endl;
	 }
	 for (int upr: upd.R)
		 city.status[upr]=2;
	 for (int upp1: upd.vP1)
		 city.status[upp1]=9101;
	 for (int upp2: upd.vP2)
		 city.status[upp2]=9102;
	 for (int upsi: upd.vSI) {
		 city.status[upsi]=912;
		 city.onset_time[upsi]=it;
		 simulationFiles.symptom_list_file << upsi << "," << city.Age[upsi]  << ","<< vaccination.vaccinated[upsi]  <<"," << it << "," << real << endl;
	 }
	 for (int upai: upd.vAI) {
		 city.status[upai]=913;
		 simulationFiles.asymptom_list_file << upai << "," << city.Age[upai]  << ","<< vaccination.vaccinated[upai]  <<"," << it << "," << real << endl;
	 }
	 for (int upr: upd.vR)
		 city.status[upr]=92;

	 // update detection, isolation, ct_done
	// copy isolation, detection, etc. status from prev timestep
	city.isolation[it] = city.isolation[it-1];
	// update isolated
	 for (int upisol: upd.isolated) {
		 city.isolation[it][upisol] =1;
		 city.iso_time[upisol]=it;
	 }
	 //update detected
	 for (int updetec:upd.detected) {
		 city.detected[updetec]=1;
	 }

}

/**
 *
 * update compartments :
 *  typeComputeImmunity =0 -> regular simulation.
 *  typeComputeImmunity =1 -> constrained simulation - never more than init_incid incidence
 *
 */
void updateCompartmentsWithFiltering(int real, int it, Compartments & upd, Population & city, Vaccination & vaccination,Params & params,
		default_random_engine & generator, SimulationFiles & simulationFiles, int typeComputeImmunity) {

	 // UPDATE TRANSMISSION
	// what we do here is that we limit the number of people going to E so that incidence is always controlled
	// we shuffle first
	shuffle(upd.E.begin(),upd.E.end(),generator);
	// then infect people
	int maxInf = city.N_use;
	// put a cap if required
	if (typeComputeImmunity ==1) maxInf *= params.init_incid_100000/100000;
	if ((int) upd.E.size() < maxInf) maxInf = upd.E.size();
	simulationFiles.miscellaneous_file << "time " << it << " infected " << maxInf << " out of " << upd.E.size() <<endl;
	for (int k=0; k < maxInf; k++)
		 city.status[upd.E[k]]=10;
	 for (int upp1: upd.P1)
		 city.status[upp1]=101;
	 for (int upp2: upd.P2)
		 city.status[upp2]=102;
	 for (int upsi: upd.SI) {
		 city.status[upsi]=12;
		 city.onset_time[upsi]=it;
		 simulationFiles.symptom_list_file << upsi << "," << city.Age[upsi]  << "," << vaccination.vaccinated[upsi] << "," << it << "," << real << endl;
	 }
	 for (int upai: upd.AI) {
		 city.status[upai]=13;
		 simulationFiles.asymptom_list_file << upai << "," << city.Age[upai]  << "," << vaccination.vaccinated[upai] << "," << it << "," << real << endl;
	 }
	 for (int upr: upd.R)
		 city.status[upr]=2;

	 // update detection, isolation, ct_done
	// copy isolation, detection, etc. status from prev timestep
	city.isolation[it] = city.isolation[it-1];
	// update isolated
	 for (int upisol: upd.isolated) {
		 city.isolation[it][upisol] =1;
		 city.iso_time[upisol]=it;
	 }
	 //update detected
	 for (int updetec:upd.detected) {
		 city.detected[updetec]=1;
	 }

}



/**
 *
 * discharge isolated from isolation
 *
 */
void dischargeFromIsolation(int it,Population & city,
        Params & params, default_random_engine & generator) {

     // DISCHARGE FROM ISOLATION

    // Initialisation of immunity at random
    uniform_real_distribution<double> distribution(0,1);
     for (int i{0} ; i<city.N_use ; ++i)
     {
    	 int status= city.status[i];
    	 int isIsolated= city.isolation[it-1][i];
    	 int durIsol = it-city.iso_time[i];
    	 if (isIsolated==1) {
    		 // if infected (not susceptible not recovered)
    		 if (
    				 //infected not symptomatic isolated for > dur_isol_inf
    				 ((status!=0 && status!=90 && status!=909 && status!=908 && status!=2 && status!=912 ) && (durIsol> params.dur_isol_inf)) ||
							 //not infected isolated for > dur_isol_sus
							 ((status==0 || status==90 || status==909 || status == 908) && (durIsol> params.dur_isol_sus)) ||
									 //recovered isolated for > dur_isol_rec
									 ((status==2 || status==902) && (durIsol> params.dur_isol_rec))
    		 )
    		 {
    			 //release from isolation
    			 city.isolation[it][i]=0;
    			 city.iso_time[i]=0;
    		 }
    		 // If a node is isolated and is not a clinical case (either vaccinated or not) it may skip isolation
    		 if (distribution(generator)<params.skip_iso && status!=12 && status!=912)
    		 {
    			 city.isolation[it][i]=0;
    			 city.iso_time[i]=0;
    		 }
    	 }
     }
}

/**
 *
 * run a realisation of simulation with Seeding in city.immun*
 * starts with intialising everyting
 * then loop on steps
 *
 */
int runRealisationWithSeeding (int real, Population & city, Networks & contacts, Places & places,
		Compartments & epid, Incidences & incidence,
		Params & params,
		default_random_engine & generator,SimulationFiles& simulationFiles,int nSim) {

	// Reinitialise compartment time series, setting all to zero
	initEpid(epid, params);

	// Reinitialise nodes' variables related to detection, isolation, contact tracing, vaccination and screening
	initCity(city, params);
    // Itialise immunity - fill in status
    initImmunity(city, params, generator, simulationFiles);

     // Reinitialise vaccination: set vaccine opinion, set vaccines' counting to zero
    // vaccine elderlies if required
	Vaccination vaccination;
    initVaccination(vaccination, city, places, params, generator,simulationFiles);
	
    Testing testing;
    initTesting(testing, city);

    Screening screening;
    initRegularScreening(screening, places, params,generator,simulationFiles);


    // seed - to avoid canceling
	initSeed(city, params, generator,simulationFiles);

	// Compute initial prevalence
    computePrevalence(0, epid, city,vaccination,testing);

    initImprovedCountInCases(city); // set counts to 0
   // START REALISATION

	if (params.incrVaccAccept_inplace ==1 ) {
		params.doIncreasedVaccination=1; // for reinforced vaccination
		cout << "vax acceptance turned to 100% ";
	} else {
		params.doIncreasedVaccination=0; // no reinforced vaccination
		cout << "vax acceptance untouched ";
	}

    for(int t{1} ; t<params.max_steps ; ++t) {
		
        // Create empty variables to update nodes' status
		Compartments update;

		// Put in place vaccination strategy (if any)
		vaccinationStep( real,  t, vaccination, city, places,  epid,  params, generator, simulationFiles);
		
        // Test + Trace
		testingStep(real, t, testing, city, contacts, update, params,  generator, simulationFiles);

		// screening
		screeningStep(real, t, screening, testing, city, contacts, places, update, params,  generator, simulationFiles);

		// Update status of infected people and infect susceptibles
		transmissionStep( real,  t,  city, vaccination, contacts, places,update, params, generator, simulationFiles);
		
		//this must be done after update is filled in
		countCasesInWorkingPlacesAndSchools(t, update, screening, vaccination, city, places, params);
        // apply updats to city
		updateCompartments(real, t, update, city,vaccination,simulationFiles);
	    // get out of isolation - moved back here for updates
	    dischargeFromIsolation(t, city, params, generator);

		//compute current prevalence
		computePrevalence(t, epid, city,vaccination, testing);
		// compute incidence
		computeIncidence (real,t, incidence, update);

		//END OF TIMESTEP

		// Rotate the network
		if (contacts.idxNetwork >= params.number_networks-1){
			contacts.idxNetwork=0;
		} else {
			contacts.idxNetwork++;
		}
		cout << ".";
	}
	cout <<"\n";
	// --------------------------------------------------------------------------------------------------------------
	// END OF SINGLE REALIZATION
	// --------------------------------------------------------------------------------------------------------------
	writeVaccinated(vaccination, city, real, simulationFiles);
	return(1);
}

/**
 *
 * run a realisation of simulation, starting from a few cases (3)
 * starts tracing and vaccination only after a preset number of cases is attained
 *
 */
int runRealisationWithInvasion(int real, Population & city, Networks & contacts, Places & places,
		Compartments & epid, Incidences & incidence,
		Params & params,
		default_random_engine & generator,SimulationFiles& simulationFiles, int nSim) {

	// Reinitialise compartment time series, setting all to zero
	initEpid(epid, params);

	// Reinitialise nodes' variables related to detection, isolation, contact tracing, vaccination and screening
	initCity(city, params);
    // Itialise immunity - fill in status
    initImmunity(city, params, generator, simulationFiles);

     // Reinitialise vaccination: set vaccine opinion, set vaccines' counting to zero
    // vaccine elderlies if required
	Vaccination vaccination;
    initVaccination(vaccination, city, places, params, generator,simulationFiles);

    // we start testing
    Testing testing;
    initTesting(testing, city);

    Screening screening;
    initRegularScreening(screening, places, params,generator,simulationFiles);
    // Compute initial prevalence

    // seed - we do this after to avoid cancelling seeds
	initSeed(city, params, generator,simulationFiles);

	computePrevalence(0, epid, city,vaccination,testing);

    initImprovedCountInCases(city); // set counts to 0
   // START REALISATION

    //reset increased Vaccination
    params.doIncreasedVaccination=0;


    int afterFirstDetection=0;
    int afterThresholdDetection=0;
    int runOK=1; // 1 if epidemic does not go extinct before seeding
    int t_thr = params.max_steps;
    int status_thr=0; // this is 0 if censored, 1 otherwise
    int t_extinct = params.max_steps ;
	int status_extinct=0; // this is 0 if censored, 1 otherwise
	int t_first = params.max_steps;
	int status_first=0;

	// store original parameters that may be changed
	int orig_skip_iso = params.skip_iso ;
	double orig_p_c_HH = params.p_c_HH ;
	double orig_p_rA = params.p_rA ;
	double orig_p_rS = params.p_rS ;
	double orig_r_d_asymp =params.r_d_asymp;
	double orig_r_d_symp=params.r_d_symp ;
	double orig_will_be_d_asymp =params.will_be_d_asymp;
	double orig_will_be_d_symp =params.will_be_d_symp;
	int orig_incrVaccAccept_inplace=params.incrVaccAccept_inplace ;

    for(int t{1} ; t<params.max_steps ; ++t) {

        // Create empty variables to update nodes' status
		Compartments update;

		// we look whether the number of detected is larger than XX
		if (epid.detected[t-1] >= 1  && afterFirstDetection == 0) {
			afterFirstDetection=1;
			cout << " first detection at time " << t-1 << " with " << epid.detected[t-1] <<endl;
			t_first=t-1;
			status_first=1;
		}
		// we look whether the number of detected is larger than XX
		if (epid.detected[t-1] >= params.n_start_reactive  && afterThresholdDetection == 0) {
			afterThresholdDetection =1;
			cout << " threshold met at time " << t-1 << " with " << epid.detected[t-1] <<" detections ";
			if (params.incrVaccAccept_inplace ==1 ) {
				params.doIncreasedVaccination=1; // for reinforced vaccination
				cout << "vax acceptance turned to 100% ";
			} else {
				params.doIncreasedVaccination=0; // no reinforced vaccination
				cout << "vax acceptance untouched ";
			}
			t_thr = t-1;
			status_thr=1;
			cout << "changed parameters"<<endl;
		}

		// change parameters for tracing after detection
		if (afterFirstDetection==1) {
			params.skip_iso = params.skip_iso_aft_inv;
			params.p_c_HH = params.p_c_HH_aft_inv;
			params.p_rA = params.p_rA_aft_inv;
			params.p_rS = params.p_rS_aft_inv;
			params.r_d_asymp = params.r_d_asymp_aft_inv;
			params.r_d_symp = params.r_d_symp_aft_inv;
			params.will_be_d_asymp = params.will_be_d_asymp_aft_inv;
			params.will_be_d_symp = params.will_be_d_symp_aft_inv;

		}

		// Put in place vaccination strategy (if any - only after invasion)
		if (afterThresholdDetection==1) {
			vaccinationStep( real,  t, vaccination, city, places,  epid,  params, generator, simulationFiles);
		} else if (afterThresholdDetection == 0) {
			// store reactive vax
			int par_reac_vac = params.reactive_vacc_inplace;
			// force reactive vac out
			params.reactive_vacc_inplace = 0;
			vaccinationStep( real,  t, vaccination, city, places,  epid,  params, generator, simulationFiles);
			//set to original reac vac
			params.reactive_vacc_inplace = par_reac_vac;
		}

		// we detect in testingStep
        // Test + Trace - we start testing and tracing at the beginning
		testingStep(real, t, testing, city, contacts, update, params,  generator, simulationFiles);

        // Screen
		if (afterThresholdDetection==1) screeningStep(real, t, screening, testing, city, contacts, places,update, params,  generator, simulationFiles);

		// Update status of infected people and infect susceptibles
		transmissionStep( real,  t,  city, vaccination, contacts, places,update, params, generator, simulationFiles);

		//this must be done after update is filled in
		countCasesInWorkingPlacesAndSchools(t, update, screening, vaccination, city, places, params);
        // apply updates to city
		updateCompartments(real, t, update, city,vaccination,simulationFiles);
	    // get out of isolation - moved back here for updates
	    dischargeFromIsolation(t, city, params, generator);

		//compute current prevalence
		computePrevalence(t, epid, city,vaccination, testing);

		// compute incidence
		computeIncidence (nSim,t, incidence, update);

		// check for extinction
		if (epid.infected[t]==0 && epid.infected[t-1]>0) {
			t_extinct=t;
			status_extinct=1;
		}

		//END OF TIMESTEP

		// Rotate the network
		if (contacts.idxNetwork >= params.number_networks-1){
			contacts.idxNetwork=0;
		} else {
			contacts.idxNetwork++;
		}
		cout << ".";
	}
	cout <<"\n";
	// --------------------------------------------------------------------------------------------------------------
	// END OF SINGLE REALIZATION
	// --------------------------------------------------------------------------------------------------------------
	// if runOK, add to epid, incidence and write, otherwise dismiss
	if (runOK==1) {
		simulationFiles.extinction_file << real << "," << nSim <<"," << t_first << "," << status_first << "," << t_thr << "," << status_thr << "," << t_extinct << "," << status_extinct <<endl;
		writeVaccinated(vaccination, city, real, simulationFiles);
		simulationFiles.reactive_vaccination_file<<real<<","<<nSim;
		for (int j=0; j<params.max_steps;j++) {
			simulationFiles.reactive_vaccination_file <<","<<epid.vaccines_used_reactive[j];
		}
		simulationFiles.reactive_vaccination_file<<endl;
	}

	// set parameters back to normal
	params.incrVaccAccept_inplace = orig_incrVaccAccept_inplace	;
	params.skip_iso = orig_skip_iso;
	params.p_c_HH = orig_p_c_HH;
	params.p_rA = orig_p_rA;
	params.p_rS = orig_p_rS;
	params.r_d_asymp = orig_r_d_asymp;
	params.r_d_symp = orig_r_d_symp;
	params.will_be_d_asymp = orig_will_be_d_asymp;
	params.will_be_d_symp = orig_will_be_d_symp;

	return(runOK);
}

/**
 *
 * run the simulation with no intervention and write status files at given intervals of prevalence in all infected
 *
 */
void runRealisationForImmunity(int real, Population & city, Networks & contacts, Places & places,
		Compartments & epid,
		Params & params,
		default_random_engine & generator,SimulationFiles& simulationFiles,
		int typeComputeImmunity) {

		//typeComputeImmunity : 0 for wild epidemic; 1 for constrained epidemic
	cout << "updating networks" << endl;
	initNetworks(contacts, city.N_use,places, params, generator);
	// we want no interventions :
	params.immunity_option=0;
	cout << "set immunity to 0"<<endl;
	params.HHisolation_inplace=0;
	params.manualCT_inplace=0;
	params.vac_rand_inplace=0;
	params.vac_elderly_inplace=0;
	params.regular_screening_WP_inplace=0;
	params.regular_screening_S_inplace=0;
	params.regular_screening_HH_inplace=0;
	params.reactive_screening_S_inplace=0;
	params.reactive_screening_WP_inplace=0;
	params.reactive_vacc_inplace=0;
	cout << "set all interventions to 0"<<endl;
	// Reinitialise compartment time series, setting all to zero
	initEpid(epid, params);

	// Reinitialise nodes' variables related to detection, isolation, contact tracing, vaccination and screening
	initCity(city, params);
    // Itialise immunity - fill in status
    initImmunity(city, params, generator, simulationFiles);
    // seed
	initSeed(city, params, generator,simulationFiles);

     // Reinitialise vaccination: set vaccine opinion, set vaccines' counting to zero
    // vaccine elderlies if required
	Vaccination vaccination;
	//
    initVaccination(vaccination, city, places, params, generator,simulationFiles);

    Testing testing;
    initTesting(testing, city);
    Screening screening;
    initRegularScreening(screening,places,  params,generator,simulationFiles);
    // Compute initial prevalence
    computePrevalence(0, epid, city,vaccination,testing);

    initImprovedCountInCases(city); // set counts to 0
   // START REALISATION

    int next_prev=2; // step in prevalence for output

    int cont=0;
    int t=0;
    while (cont == 0) {
    	t++;
        // Create empty variables to update nodes' status
		Compartments update;
		// Update status of infected people and infect susceptibles
		transmissionStep( real,  t,  city, vaccination, contacts,places, update, params, generator, simulationFiles);

        // apply updats to city
		updateCompartmentsWithFiltering(real, t, update, city,vaccination,params,generator, simulationFiles, typeComputeImmunity);
	    // get out of isolation - moved back here for updates

		//compute current prevalence
		computePrevalence(t, epid, city,vaccination, testing);
		// compute incidence
//		computeIncidence (real,t, incidence, update);

		// write prevalence file
		int prev = (epid.P1[t]+epid.P2[t]+epid.SI[t]+epid.AI[t]+epid.R[t])*100/city.N_use;
		if (prev>= next_prev) {
			writeImmunity(t,city ,prev);
			next_prev=prev+2;
			if (next_prev == 30) cont=1;
		}
		//END OF TIMESTEP

		// Rotate the network
		if (contacts.idxNetwork >= params.number_networks-1){
			contacts.idxNetwork=0;
		} else {
			contacts.idxNetwork++;
		}
		cout << ".";
	}
	cout <<"\n";
	// --------------------------------------------------------------------------------------------------------------
	// END OF SINGLE REALIZATION
	// --------------------------------------------------------------------------------------------------------------
}


/**
 *
 * add compartments to sum over all realisations
 *
 */
void addToAllEpid(Compartments & allEpids, Compartments & epid, Params & params) {
	// add to all compartments result of current  simulation

	 for (int i{0} ; i<=params.max_steps ; ++i)
	        {
	            allEpids.S[i] += epid.S[i];
	            allEpids.E[i]+= epid.E[i];
	            allEpids.P1[i]+= epid.P1[i];
	            allEpids.P2[i]+= epid.P2[i];
	            allEpids.SI[i]+= epid.SI[i];
	            allEpids.AI[i]+= epid.AI[i];
	            allEpids.R[i]+= epid.R[i];
	            allEpids.wvS[i]+= epid.wvS[i];
	            allEpids.vS[i]+= epid.vS[i];
	            allEpids.vE[i]+= epid.vE[i];
	            allEpids.vP1[i]+= epid.vP1[i];
	            allEpids.vP2[i]+= epid.vP2[i];
	            allEpids.vSI[i]+= epid.vSI[i];
	            allEpids.vAI[i]+= epid.vAI[i];
	            allEpids.vR[i]+= epid.vR[i];
	            // average isolations
	            allEpids.WP_infected[i]+= epid.WP_infected[i];
	            allEpids.S_infected[i]+= epid.S_infected[i];
	            allEpids.isolated[i]+= epid.isolated[i];
	            allEpids.infected_and_isolated[i]+= epid.infected_and_isolated[i];
	            allEpids.vaccines_used[i]+= epid.vaccines_used[i];
	            allEpids.vaccines_used_random[i]+= epid.vaccines_used_random[i];
	            allEpids.vaccines_used_reactive[i]+= epid.vaccines_used_reactive[i];
	            allEpids.vaccinated_FUL[i]+= epid.vaccinated_FUL[i];
	            allEpids.vaccinated_MID[i]+= epid.vaccinated_MID[i];
	            allEpids.vaccinated_LOW[i]+= epid.vaccinated_LOW[i];
	            allEpids.detected[i]+= epid.detected[i];
	            allEpids.found_by_ct[i]+= epid.found_by_ct[i];
	            allEpids.found_by_regular_screening[i]+= epid.found_by_regular_screening[i];
	            allEpids.found_by_reactive_screening[i]+= epid.found_by_reactive_screening[i];
	        }
}


#define OPTIONS "hs:n:f:i:j:x:"

/**
 *
 * main program
 *
 */
int main(int argc, char* argv[]) {

    int opt;
    int seed=0;
    int scenario=0;
    string fichierInput="filenames.txt";

    int simulInvasion=0;
    int typeSeeding=0;
    int computeImmunity=0;
    int typeComputeImmunity=0;
    int n_start_reactive=0;
    // Retrieve the options:
    while ( (opt = getopt(argc, argv, OPTIONS)) != -1 ) {  // for each option...
    	switch ( opt ) {
    	case 'h':
    		cout << argv[0] << endl;
    		cout <<"\t-h : print this help"<<endl;
    		cout << "\t-n N runs scenario N"<<endl;
    		cout <<	"\t-s SEED sets random seed to SEED (0 for time)"<<endl;
			cout << "\t-f FILENAME for input (default filenames.txt)"<<endl;
			cout <<"\t-i 0 for computing (natural) initial immunity/ -i 1 constrained to max incidence"<<endl;
			cout << "\t-j 0 (random seed) / 1 (seed only exposed) / 2 (seed exposed and others) / 3 (seed all exposed file) / 4 (same as 2 but random)" <<endl;
			exit(1);
    	case 's':
    		seed = atoi(optarg);
    		break;
    	case 'n':
    		scenario=atoi(optarg);
    		cout << "run scenario " << scenario << endl;
    		break;
    	case 'f' :
    		fichierInput=optarg;
    		cout << "reading filenames in " << fichierInput << endl;
    		break;
    	case 'i' :
    		typeComputeImmunity = atoi(optarg);
    		computeImmunity=1;
    		cout << "computing immunity type " <<typeComputeImmunity <<endl; ;
    		break;
    	case 'j' :
    		typeSeeding = atoi(optarg);
    		cout << "seeding is type " << typeSeeding <<endl;
    		break;
    	case 'x':
    		simulInvasion=1;
    		cout << "simulating invasion : ";
    		n_start_reactive = atoi(optarg);
    		cout << "start reactive vaccination after " << n_start_reactive << " cases"<<endl;
    		break;
		case '?' :
    		cout << "unknown option" <<endl;
    		exit(1);
    		break;
    	}
    }


    if (scenario==0 && computeImmunity==0) {
    	cout << "scenario number " << scenario << " is not OK\n";
    	exit(1);
    } else if (scenario==0 && computeImmunity ==1){
    	cout << "need a no_intervention scenario" <<endl;
    	cout << "run 'grep -n no_intervention input/params.csv' "<<endl;
    	exit(1);
    }

    time_t seed_gen;
    if (seed == 0) {
    	seed_gen = time(NULL); //lit le seed
    	cout << "use random seed : " << (int) seed_gen<<"\n";
    } else {
    	seed_gen = (time_t) seed;
    	cout << "use fixed seed : " << (int) seed_gen<<"\n";
    }
    default_random_engine generator(seed_gen); //Random number itialisaton
    uniform_real_distribution<double> distribution(0,1); //Intialisation of the uniform distribution between 0 and 1

    	//pointer to realisation function
    int (* runRealisation) (int, Population & , Networks & , Places &,
    		Compartments & , Incidences & ,
    		Params & ,
    		default_random_engine & ,SimulationFiles&, int );


    // READ PARAMETERS FROM OUTSIDE
    // filenames are in fichierInput and read here
    // this allow using long names
    Filenames filenames;
    readFileNames(fichierInput, filenames);
    int lineNumberSought= scenario;
    string scenario_tag;
    MapParams mapParams;
    Params params;
    Population city;
    Networks contacts;
    Places places;
    
    // READ PARAMS
    readParams(filenames, lineNumberSought, scenario_tag, mapParams);
    // parameters to move to file
    mapParams.insert({"pct_rand_init_immunity",1./3.});
    mapParams.insert({"n_start_reactive",n_start_reactive});
    // change to list
    convertParams( params, mapParams);
    params.typeSeeding = typeSeeding;
    params.computeImmunity = computeImmunity;

    //
    if (simulInvasion == 1) {
       	runRealisation = runRealisationWithInvasion;
       	params.doIncreasedVaccination=0;
       	if (params.incrVaccAccept_inplace==1) {
       		cout << "params.doIncreasedVaccination will be set to 1 after threshold"<<endl;
       	} else {
       		cout << "params.doIncreasedVaccination will be set to 0 after threshold"<<endl;
       	}
   		typeSeeding=5;
   		if (params.n_start_reactive<1) {
   			params.n_start_reactive=1;
   			cout << "n_start_reactive set to 1"<<endl;
   		}
       } else {
       	runRealisation = runRealisationWithSeeding;
       	if (params.incrVaccAccept_inplace==1) {
       		cout << "params.doIncreasedVaccination will be set to 1"<<endl;
       	} else {
       		cout << "params.doIncreasedVaccination will be set to 0"<<endl;
       	}
       }
    printParams(params);
//to check
    // to store :
    Params origParams = params;

    // READ AND INITIALISE POPULATION PROPERTIES - AGE SUSCEPTIBILITY ASYMPTO IMMUNITY
    int numFirstIndv; //number of first indiv
    int numLastIndv; //number of last indiv
    readPopulation(filenames, city, &numFirstIndv, &numLastIndv, params);
    
    // READ THE NETWORK FILE
    // uses number first/last read above fills in nodes, place,freq
	readNetworkFile(filenames, numFirstIndv,numLastIndv,contacts, city.N_use) ;
    
    // READ HOUSEHOLD, SCHOOL AND WORKPLACE INFORMATION
    readHouseholdSchoolWorkplaces(filenames,places);
    
    //  INITIALISATION
    initPlaces(places, city, params); // define places that are too young to vaccinate; too young to test

    // TO STORE RESULTS of all epidemics (summed over realisations)
    Compartments allEpids; //to store and average
    initEpid(allEpids,params);
    
    // to hold value of incidence over all realisations
    Incidences incidence;
    initIncidence(incidence,params);
    
    // resize network and compute acquaintances - building dynamic network defered to main loop
    initDynamicalNetworks(contacts,city.N_use, params);
    
    //open FILES to be used during simulations : these files are only closed at exit/saves IO
    SimulationFiles simulationFiles;
    openSimulationFiles(simulationFiles, scenario_tag);
    
    // To store final number of recovered in each realization
    vector<double> R_final(params.N_real); // this could be somewhere else
    
    //compute files for initial immunity : simulate epidemic then save status and exit
    if (computeImmunity==1) {
    	Compartments epid;
    	runRealisationForImmunity(0, city,contacts, places, epid, params, generator,simulationFiles, typeComputeImmunity);
    	exit(1);
    }

    vector<int> idxSim;
    //  START LOOPING OVER SIMULATIONS
    int nSimOK=0; //count simulations that starts (not extinct before )
    for (int real{0}; real<5*params.N_real; ++real) //max we make 5 times the required
    {
 //   	params = origParams;
        if (nSimOK==params.N_real) break;

    	int simOK=0;
        // UPDATE DYNAMICAL NETWORK - DONE ON FIRST LOOP, THEN EVERY refresh_network SIMS
    	if ((real==0) || (nSimOK>0 && !(nSimOK % (int) params.refresh_network))) {
    		cout << "updating networks" << endl;
    		initNetworks(contacts, city.N_use,places, params, generator);
    	}
        cout << "realisation number " << nSimOK << " - attempted " << real << endl;
        
        // OUTPUT OF ONE REALISATION
        Compartments epid;
        //simOK is always 1
        simOK=runRealisation(real, city,contacts, places, epid, incidence, params, generator,simulationFiles,nSimOK);
        
        if (simOK==1) {
        	idxSim.push_back(real);
        	nSimOK++;
        	// SUMMARIZE REALISATION
        	addToAllEpid(allEpids, epid, params);
        	R_final[nSimOK] = ((epid.R[epid.R.size()-1]+epid.vR[epid.vR.size()-1])*1.0)/city.N_use;
        } else {
        	if (((nSimOK==0) && (real > params.refresh_network)) || ((nSimOK>0) && (real/nSimOK > params.refresh_network))) {
        		cout << "initialization is not working properly: too many extinctions. exiting."<<endl;
        		break;
        	}
        }
    }
    // WRITE ALL FILES
    writeAverageFiles(city.N_use, allEpids, incidence, R_final, scenario_tag,  params);
    // CLOSE ALL FILES
    closeSimulationFiles(simulationFiles);
    // write the idx of Simulations
    writeIdxSim(idxSim,scenario_tag);

    exit(0);
}

