/**
    Simulation of an ABM for COVID-19
    @file transmission.cpp
    @author Pierre-Yves Boelle, Chiara Poletto
    @acknowledgment Jesus Moreno for earlier version
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


/*
 * compute probabilities not infectious according to characteristics
 */
double computepNotInf(double pInfStub, int statusContact, int noneIsolated, int isVaccinated, int inHousehold,
		Params & params) {
	// this function computes P(Infection) for all layers
	// pInfStub : baseline prob for layer
	// statusContact : infectious status of alter
	// noneIsolated is true is isol(ego) =0 and isol(alter)=0
	// isVaccinated : for ego this can be VAX_NOT/VAX_LOW/VAX_MID/VAX_FUL depending on vaccine status
	// in Household : true for layer 0, false otherwise
	// params vector of params

	double pInf=-1.0; //negative by default

	if (!(inHousehold) && !(noneIsolated)) return(pInf); // outside household, no transmission for isolated persons

    if (statusContact==102  || statusContact==13) {
        // if neighbour is pre-asymptomatic or asymptomatic, include re-scale in transmissibility because of asymptomatic route "r_a"
        if (noneIsolated) {// if neighbour and susceptible are not isolated
        	if (isVaccinated == VAX_FUL) {   // if susceptible is vaccinated, add factor v_sus
            	pInf = 1.0-pInfStub*params.r_a*params.v_sus_FUL;
            } else if (isVaccinated == VAX_MID) {
              pInf = 1.0-pInfStub*params.r_a*params.v_sus_MID;
            } else {
            	// if susceptible is not vaccinated, no factor v_sus
            	pInf = 1.0-pInfStub*params.r_a;
            }
        } else {  // if either of them are isolated, due to being the household layer, there can be a contact and a possible infection but there is an isolation factor "factor_isol" added
            if (isVaccinated == VAX_FUL)  {
            	pInf=1.0-pInfStub*params.r_a*params.v_sus_FUL*params.factor_isol;
            } else if (isVaccinated == VAX_MID)  {
            	pInf=1.0-pInfStub*params.r_a*params.v_sus_MID*params.factor_isol;
            } else {
            	pInf=1.0-pInfStub*params.r_a*params.factor_isol;
            }
        }
     } else if (statusContact==9102  || statusContact==913) {
    	// same as before but if contacts are pre-asymptomatics or asymptomatics vaccinated, add factor v_trans in transmission
        if (noneIsolated) {
            if (isVaccinated == VAX_FUL){
            	pInf=1.0-pInfStub*params.r_a*params.v_sus_FUL*params.v_trans;
            } else  if (isVaccinated == VAX_MID){
            	pInf=1.0-pInfStub*params.r_a*params.v_sus_MID*params.v_trans;
            } else {
            	pInf=1.0-pInfStub*params.r_a*params.v_trans;
            }
        } else {
            if (isVaccinated == VAX_FUL) {
            	pInf=1.0-pInfStub*params.r_a*params.v_sus_FUL*params.factor_isol*params.v_trans;
            } else  if (isVaccinated == VAX_MID) {
            	pInf=1.0-pInfStub*params.r_a*params.v_sus_MID*params.factor_isol*params.v_trans;
            } else {
            	pInf=1.0-pInfStub*params.r_a*params.factor_isol*params.v_trans;
            }
        }
    } else if (statusContact==101 || (statusContact==12 && inHousehold)) {
        // if neighbour is pre-symptomatic or symptomatic
        if (noneIsolated) {
        	if (isVaccinated == VAX_FUL) {
        		pInf = 1.0 - pInfStub*params.v_sus_FUL;
        	} else if (isVaccinated == VAX_MID) {
        		pInf = 1.0 - pInfStub*params.v_sus_MID;
        	} else{
        		pInf=1.0-pInfStub;
        	}
        } else {
        	if (isVaccinated == VAX_FUL) {
        		pInf=1.0- pInfStub*params.v_sus_FUL*params.factor_isol;
        	} else	if (isVaccinated == VAX_MID) {
        		pInf=1.0- pInfStub*params.v_sus_MID*params.factor_isol;
        	} else{
        		pInf=1.0 - pInfStub*params.factor_isol;
        	}
        }
    } else if (statusContact==9101 || (statusContact==912 && inHousehold)) {
    	// vaccinated + presymptomatic or symptomatic
        if (noneIsolated) {
            if (isVaccinated==VAX_FUL){
            	pInf= 1.0-pInfStub*params.v_sus_FUL*params.v_trans;
            } else if (isVaccinated==VAX_MID){
            	pInf= 1.0-pInfStub*params.v_sus_MID*params.v_trans;
            } else {
            	pInf=1.0-pInfStub*params.v_trans;
            }
        } else {
            if (isVaccinated == VAX_FUL) {
            	pInf=1.0 - pInfStub*params.v_sus_FUL*params.factor_isol*params.v_trans;
            } else if (isVaccinated == VAX_MID) {
            	pInf=1.0 - pInfStub*params.v_sus_MID*params.factor_isol*params.v_trans;
            } else{
            	pInf = 1.0-pInfStub*params.factor_isol*params.v_trans;
            }
        }
    }
    return(pInf);
}

/*
 *
 *
 */
int addInfectorsInLayer(int i,
		vector<int> & contactLayerAtNetItAtEgo,//list of js in contact
		Population & city,
		Vaccination & vaccination,
		int layer,//1 if true
		double pInfStub,//common to layer
		int it,
		int real,
		vector<double> & p_ID, //probInfection
		vector<vector<int>> & infectious_neighbours, /*infectious neighbours*/
		Params & params) {

	int inHousehold;
	int nContacts = contactLayerAtNetItAtEgo.size();
	if (nContacts == 0) {return(0);}
	else {
		if (layer==LAYER_HH) {inHousehold=1;} else{inHousehold=0;};
		// this takes into account VAX_NOT/...
		int isVaccinated = vaccination.vaccinated[i];
//
		int isIsolated = (city.isolation[it-1][i]==1)?1:0;
        for (int j{0}; j< nContacts ; ++j)
        {
            int statusContact = (int) (city.status[contactLayerAtNetItAtEgo[j]]);
            int noneIsolated = (1 - (int) (city.isolation[it-1][contactLayerAtNetItAtEgo[j]]))*(1-isIsolated) ;
            double pNotInf=computepNotInf(pInfStub, statusContact, noneIsolated, isVaccinated, inHousehold, params);
            if (pNotInf>0.0) {
                vector<int> neighbours1;
                neighbours1.push_back(contactLayerAtNetItAtEgo[j]);
                neighbours1.push_back(city.Age[contactLayerAtNetItAtEgo[j]]);
                neighbours1.push_back(city.status[contactLayerAtNetItAtEgo[j]]);
                neighbours1.push_back(vaccination.vaccinated[contactLayerAtNetItAtEgo[j]]); // newly added
                neighbours1.push_back(layer);
                neighbours1.push_back(i);
                neighbours1.push_back(city.Age[i]);
                neighbours1.push_back(vaccination.vaccinated[i]); // newly added
                neighbours1.push_back(it);
                neighbours1.push_back(real);

            	p_ID.push_back(pNotInf);
            	infectious_neighbours.push_back(neighbours1);
            }
        }
    }
	return(1);
}

/*
 *
 */
int updateNotSusceptible(int i,//of ego
		int status, //status
		double a_propor,
		double va_propor,
		Compartments & upd,
		Params &params,
		default_random_engine& generator) {
    
    uniform_real_distribution<double> distribution(0,1);
    // Look at exposed and determine if they leave exposed phase
    if (status==10 && distribution(generator)<params.e_rate)
    {
        // determine if they become subclinical
        if (distribution(generator)<a_propor)
            upd.P2.push_back(i);
        else
            upd.P1.push_back(i);
    }
    else if (status==910 && distribution(generator)<params.ve_rate) {
        // Look at vaccinated exposed and determine if they leave exposed phase
        if (distribution(generator)<va_propor)
            upd.vP2.push_back(i);
        else
            upd.vP1.push_back(i);
    }
    else if (status==101 && distribution(generator)<params.p1_rate) {
        // Look at presymptomatic 1 and determine if the leave to be symptomatic
        upd.SI.push_back(i);
    }
    else if (status==9101 && distribution(generator)<params.vp1_rate) {
        // Look at vaccinated presymptomatic 1 and determine if the leave to be symptomatic
        upd.vSI.push_back(i);
    }
    else if (status==102 && distribution(generator)<params.p2_rate) {
        // L	ook at presymptomatic 2 and determine if the leave to be asymptomatic
        upd.AI.push_back(i);
    }
    else if (status==9102 && distribution(generator)<params.vp2_rate) {
        //Look at vaccinated presymptomatic 2 and determine if the leave to be asymptomatic
        upd.vAI.push_back(i);
    }
    else if (status==12 || status==13) {
        // Look at symptomatic or asymptomatic to determine if they recover
        if (distribution(generator)<params.gamm)
            upd.R.push_back(i);
    }
    else if (status==912 || status==913) {
        // Look at vaccinated symptomatic or asymptomatic to determine if they recover
        if (distribution(generator)<params.vgamm)
            upd.vR.push_back(i);
    } else 	if (status==909) {
		if (distribution(generator)<params.w_rate_MID_to_FUL) {
	    	//PYB - moved to updateNotSusceptible
			// Look at vaccinated waiting to determine if they become vaccinated
	    		upd.vS.push_back(i);
		}
	} else if (status==908) {
		if(distribution(generator)<params.w_rate_LOW_to_MID) {
	    	//PYB - moved to updateNotSusceptible
			// Look at vaccinated waiting to determine if they become vaccinated
	    		upd.vS.push_back(i);
		}
	}

	return(1);
}

/*
 *
 */
void updateSusceptible(int i,
		Population &city,
		Vaccination & vaccination,
		Networks & contacts,
		Places & places,
		int it, int real,
		Compartments & upd,
		Params &params,
		default_random_engine& generator,
		SimulationFiles & simulationFiles) {

	double age_w = city.age_w[i];
	int isVaccinated = vaccination.vaccinated[i];
	uniform_real_distribution<double> distribution(0,1);

	vector<vector<int>> infectious_neighbours;
	vector<double> p_ID;
	// probability of infection p_ID = beta_overall*beta_node*beta_layer
	for (int layer=0; layer<5; layer ++) {
		addInfectorsInLayer(i,contacts.layers[layer][contacts.idxNetwork][i],
				city, vaccination, layer, params.beta*contacts.weights_layer[layer]*age_w,
				it,real,
			p_ID, infectious_neighbours,
			params);
	}
	//  if risk of contagion > 0 decide if node i gets infected
	if (p_ID.size()>0) {
		double prob_ID{1.0};
		for (double pr: p_ID)
		{prob_ID = prob_ID*(pr);}

		double r_final{distribution(generator)};

		if (r_final<1.0-prob_ID) { // actually infected
			if (isVaccinated==VAX_NOT) {
				upd.E.push_back(i);
			} else if (isVaccinated==VAX_LOW) {
				upd.E.push_back(i);
//				vaccination.vaccinated[i] = VAX_FUL; - not so that we keep info on status at infection
			} else if (isVaccinated==VAX_MID) {
				if (distribution(generator)<params.E_ve_rate) {
					upd.E.push_back(i);
				} else {
					upd.vE.push_back(i);
				}
//				vaccination.vaccinated[i] = VAX_FUL;  - not so that we keep info on status at infection
			} else if (isVaccinated==VAX_FUL) {
				upd.vE.push_back(i);
			}
		
            if (params.print_infectors==1) {
                //Determine the infector from sampling of p_ID
                vector<int> infector{0,0,0,0,0,0,0};
                vector<double> cumulative;
                cumulative.push_back(0);
                double norm{0.0};
                for (double p: p_ID) norm=norm+p;
                double cumul{0.0};
                for (double p: p_ID) {
                    cumul=cumul+p;
                    cumulative.push_back(cumul/norm);
                }
                double r{distribution(generator)};
                for (int k=1 ; k < (int) cumulative.size() ; ++k) {
                    if (r>cumulative[k-1] && r<cumulative[k]) { infector= infectious_neighbours[k-1]; break;}
                }
                // Print information about infector and infected
                // PYB changed to print vax information (position 3 and 7)
                simulationFiles.infector_file <<
                		infector[0] << "," << //ID infector
                		infector[1] << "," << //age infector
						infector[2] << "," << //status infector
						infector[3] << "," << //vax infector
						infector[4] << "," << //layer
						infector[5] << "," << //ID infected
						infector[6] << "," << //age infected
						infector[7] << "," ; //vax infected
//add print place of infection
                		if (infector[4]==LAYER_HH) simulationFiles.infector_file << places.node_HH[i]<<",";
                		else if (infector[4]==LAYER_WP) simulationFiles.infector_file << places.node_WP[i]<<",";
                		else if (infector[4]==LAYER_S) simulationFiles.infector_file << places.node_S[i]<<",";
                   		else if (infector[4]==LAYER_CO) simulationFiles.infector_file << "community" <<",";
                   		else if (infector[4]==LAYER_TR) simulationFiles.infector_file << "transport" <<",";

                    	simulationFiles.infector_file << infector[8] << "," << //iteration
						infector[9] << endl; //simulation
            }
        }
	}
}

// TRANSMISSION

void transmissionStep( int real, int it,  Population & city, Vaccination & vaccination, Networks & contacts, Places & places,
		Compartments &upd,Params & params, default_random_engine & generator, SimulationFiles & simulationFiles) {
	 
	 for (int i{0} ; i<city.N_use ; ++i)
	 {
		 // note 909/908 goes both ways
		 if (city.status[i]!=0 && city.status[i]!=90)
			 updateNotSusceptible(i, city.status[i], city.a_propor[i], city.va_propor[i],
					 upd,
					 params, generator);

		 if (city.status[i]==0 || city.status[i]==90 || city.status[i]==909 || city.status[i]==908)
			 // Look at susceptible, vaccinated susceptible and vaccinated waiting
			 updateSusceptible(i,city,vaccination, contacts, places, it, real,
					 upd, params,generator,simulationFiles);
	 }
}


