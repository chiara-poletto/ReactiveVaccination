
/**
    Simulation of an ABM for COVID-19
    @file defs.h
    @author Pierre-Yves Boelle, Chiara Poletto
    @version 1.0 2020-04-16
    @license GPL-3.0-or-later
*/
#ifndef DEFS_H
#define DEFS_H

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cerrno>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <utility>
#include <time.h>
#include <string>
#include <sstream>
#include <random>

#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <set>
#include <list>
#include <unordered_set>
#include <tuple>
#include <map>

using namespace std;

//--------------------------------------
//SOME STRUCTURES FOR AN EASIER LIFE
//--------------------------------------

// function to compute probability of infection
typedef std::map<std::string,double> MapParams;

typedef struct FILENAMES {
	string paramfile;
	string agefile;
	string networkfile;
	string HH_list;
	string WP_list;
	string S_list;
	string immunityfile;
	string universityfile;
	string susceptibilityfile;
	string subclinicalfile;
	string subclinicalVaccinatedfile;
	string invasionfile;
} Filenames;


#define LAYER_HH 0
#define LAYER_WP 1
#define LAYER_S 2
#define LAYER_CO 3
#define LAYER_TR 4

typedef struct NETWORKS {
	//pointeur vers la version courante du ntwork
	int idxNetwork;
	//these are changing with time
	//layer/num/ego/alter
	vector<vector<vector<vector<int>>>> layers;
	vector<set<int>> acquaintances;
	//these are fixed
	// contacts are described as node1/node2/place/freq
	vector<int> node1;
	vector<int> node2;
	vector<int> place;
	vector<double> freq;
	vector<vector<int>> maxContacts;
	int maxOfMaxContacts;
	int allContactsInLayers[5];
	int keptContactsInLayers[5];
	double weights_layer[5];  // household, work, schools, community, transport
} Networks;

//
// this could be simplified to
//  typedef 	vector<ofstream> SimulationFiles;
//
//

typedef struct FILES {
	// To print a file with details of infectors / infected pairs
    ofstream infector_file;
    // To print a file with details of isolations
    ofstream isolation_file;
    // To print a file with time of symptom onset and other informations for clinical and subclinical nodes
    ofstream symptom_list_file;
	ofstream asymptom_list_file;
    // TO print detailed information on vaccination
	ofstream outv_detailed_WPS_file;
	ofstream outv_detailed_HH_file;
	// To print detailed information on the age of the vaccinated
	ofstream out_vaccinated_age_file;
	// To print number of vaccines used
	ofstream vaccine_number_file;
	// To print places in regular screening
	ofstream 	out_regular_screening_file;
	// To print places in reactive screening
	ofstream out_reactive_screening_file;
	// To print miscellaneous info
	ofstream miscellaneous_file;
	// to print extinction info
	ofstream extinction_file;
	// to print reactive vaccination numbers
	ofstream reactive_vaccination_file;
} SimulationFiles;

typedef struct COMPARTMENTS {
	/// COMPARTMENTS are either for time series WHEN USED as epid
	// also used to store individual numbers in update
	vector<int> S;
    vector<int> E;
    vector<int> P1;
    vector<int> P2;
    vector<int> SI;
    vector<int> AI;
    vector<int> R;
    vector<int> wvS;
    vector<int> vS;
    vector<int> vE;
    vector<int> vP1;
    vector<int> vP2;
    vector<int> vSI;
    vector<int> vAI;
    vector<int> vR;
    // To monitor infection in workplaces
    vector<int> WP_infected;
    // To monitor infection in schools
    vector<int> S_infected;
    // To monitor infection overall
    vector<int> infected;
    // To monitor vaccines used in time
    vector<int> vaccines_used;
    vector<int> vaccines_used_reactive;
    vector<int> vaccines_used_random;
    // To monitor vaccinated
    vector<int> vaccinated_FUL;
    vector<int> vaccinated_MID;
    vector<int> vaccinated_LOW;
	// number isolated
    vector<int> isolated;
    //number infected and isolated
    vector<int> infected_and_isolated;
    //number detected
    vector<int> detected;
    //number found by CT
    vector<int> found_by_ct;
    //number found by CT
    vector<int> found_by_regular_screening;
    //number found by CT
    vector<int> found_by_reactive_screening;
} Compartments;
//---------------------------------------

typedef struct INCIDENCES {
	vector<vector<int>> P1;
	vector<vector<int>> SI;
	vector<vector<int>> P2;
	vector<vector<int>> AI;
	vector<vector<int>> vP1;
	vector<vector<int>> vSI;
	vector<vector<int>> vP2;
	vector<vector<int>> vAI;
} Incidences;

typedef struct PLACES {
    // Map individual --> his/her place
    // Map place --> list of individuals
    // List of places
    // Household
    map<int, string> node_HH;
    map<string, vector<int>> HH_nodes;
    set<string> HHs;
    //School
	map<int, string> node_S;
    map<string, vector<int>> S_nodes;
    set<string> Ss;
    //Workplace
	map<int, string> node_WP;
    map<string, vector<int>> WP_nodes;
    set<string> WPs;
    //
    set<string> too_small_to_test; ///these are places where everybody is too young to test or too small
    set<string> too_young_to_test; ///these are places where everybody is too young to test or too small
    set<string> too_small_to_vaccinate; ///these are places where everybody is too young to test
    set<string> too_young_to_vaccinate; ///these are places where everybody is too young to test
    set<string> universities;
} Places;

typedef struct POPULATION {
	//nombre de personnes
    int N_use;
    ////////////////////////FIXED
    // Age
    vector<int> Age;
    //susceptibility
    vector<double> age_w;
    //asymptomatic
    vector<double> a_propor;
    //asymptomatic vaccinated
    vector<double> va_propor;
    ////////////////////////DYNAMIC
    //Current status
    vector<int> status;
    // 0 is susceptible, 10 is exposed, 101 pre-symptomatic(1), 12 symptomatic,
    //102 pre-symptomatic(2), 13 asymptomatic, 2 is recovered.
    // Waiting for vaccine to act: 908
    // primoVaccinated : 909
    // Vaccination status numbers:
    // 90 is susceptible, 910 is exposed,  9101 pre-symptomatic(1), 912 symptomatic,
    // 9102 pre-symptomatic(2),  913 asymptomatic, 92 is recovered.
    //status read in file immunity
    vector<int> status_file;
    //status read in file immunity
    vector<vector<int>> exposed_file;
    // Nodes detected, compliant HH & MCO contacts are also detected
    vector<int> detected;
	// Nodes that will be detected after a delay will be passed
    vector<int> will_be_detect;
    // Isolation status of nodes in time (1= isolated; 0= isolated) :: this is max_time * N_use not N_use * max_time
    vector<vector<int>> isolation;
    // Time of symptom onset for clinical cases
    vector<int> onset_time;
    // Time when isolation begane
    vector<int> iso_time;
    // structure to keep the time at which people were detected in places
    map<string, vector<int>> timeDetectedInPlaces;
    // structure to keep count of cases in places cumulated over the last params.interv_cumul_cases
    map<string, int> numberDetectedInPlaces;
    // has WP
    vector<int> has_WP;
    //has School
    vector<int> has_school;
} Population;

typedef struct PARAMS {
	// ANY CHANGE to PARAMS MUST BE :
	// -here
	// - in convertParams (init.cpp)
	// - in printParams (inputOutput.cpp)

	//simulation
	int max_steps;
	int N_real;
	int number_networks;
	int refresh_network;

	// print output
	int print_infectors;
	int print_isolations;

	//initial setup
	int immunity_option;
	double immunity_random;
	double init_incid_100000;
	int typeSeeding;
	int computeImmunity;
	double pct_rand_init_immunity; // pct for initial random immunization

	// transmission
	double beta;
	double e_rate;
	double E_ve_rate; // P(vax) tretaed as unvaccinatedon infection
	double factor_isol;
	double p1_rate;
	double p2_rate;
	double r_a;
	double v_sus_MID; //rate of susceptibility in VAX_MID
	double v_sus_FUL; //rate of susceptibility in VAX_FUL
	double v_symp;
	double v_trans;
	double ve_rate;
	double gamm;
	double vgamm;
	double vp1_rate;
	double vp2_rate;

	// contacts
	double acquai_frequency;
	double pct_teleworking;
	double pct_comty_contact;

	// cumulatinc cases to count how many in place
	int interv_cumul_cases;

	//vaccination
	double w_rate_MID_to_FUL; //from 9 to 1 - primoVac to Vac
	double w_rate_LOW_to_MID; //from 8 to 9 - susc to primoVac
	double coverage_vacc_higher_65; //willing to vaccinate percentage 0-100
	double coverage_vacc_lower_65; // willing to vaccinate percentage 0-100
	int age_th; //vaccination threshold
	int vac_elderly_inplace; //1 if true
	int vac_init_inplace; //1 if true
	double frac_elderly_vacc; //
	double frac_init_vacc;
	double r_V;
	int react_type;
    int random_type;
	int vac_rand_inplace;
	int vac_daily_randl_pop;
	int vacc_limit_rand_pop;
	int reactive_vacc_inplace;
	int vacc_daily_react;
	int vacc_limit_react;
	int th_cluster_vacc;
	int vacc_min_WP_size; // min size for random vacc WPs
	int age_vax_elderly_th; // age threshold for initial vaccination of the old

	//screening
	int regular_screening_HH_inplace; //0 or 1
	int regular_screening_S_inplace; //0 or 1
	int regular_screening_WP_inplace; //0 or 1
	int age_min_testing; // min age for testing (int)
	double f_HH_screening; // agrees to test
	double f_S_screening; // agrees to test
	double f_WP_screening; // agrees to test
	double testing_proportion_S; // overall proportion
	double testing_proportion_WP; // overall proportion
	double testing_proportion_HH; // overall proportion
	int reactive_screening_S_inplace; //0 or 1
	int reactive_screening_WP_inplace; //0 or 1
	int interval_screening; // interval min for reactive screening
	int reactive_tests_daily; // for reactive screening
	int reactive_tests_limit; // for reactive screening
	int th_cluster_screen; // cluster for reactive testing
	double r_react_screen; // rate of reactive screening
	int screening_min_size; // min size for regular screening

	//contact tracing
	int manualCT_inplace;
	int tracing_lookback;
	double r_mCT;

	//isolation
	int HHisolation_inplace;
	int dur_isol_inf; // infected
	int dur_isol_rec; // recovered
	int dur_isol_sus; // susceptible
	int skip_iso; // pct skip islation
	double p_c_HH; // pct isolation in household
	double p_rA; // pct detection contact A
	double p_rS; // pct detection contact S
	double r_d_asymp;  // rate detectection asymp
	double r_d_symp; // rate detection symp
	double will_be_d_asymp; // overall pct detection as
	double will_be_d_symp; // overall pct detection s

	// after invasion
	int incrVaccAccept_inplace; // 1 if 100% accept vax, 0 if opinion matters
	int skip_iso_aft_inv; // pct skip islation
	double p_c_HH_aft_inv; // pct isolation in household
	double p_rA_aft_inv; // pct detection contact A
	double p_rS_aft_inv; // pct detection contact S
	double r_d_asymp_aft_inv;  // rate detectection asymp
	double r_d_symp_aft_inv; // rate detection symp
	double will_be_d_asymp_aft_inv; // overall pct detection as
	double will_be_d_symp_aft_inv; // overall pct detection s

	// weights layer
	double wgt_lyr_household;
	double wgt_lyr_workspace;
	double wgt_lyr_school;
	double wgt_lyr_community;
	double wgt_lyr_transport;

	// reactive vaccination on invasion
	int n_start_reactive;
	int doIncreasedVaccination;
} Params;


#define VAX_NOT 0
#define VAX_LOW 8
#define VAX_MID 9
#define VAX_FUL 1

typedef struct VACCINATION {
	// to hold all information related to vaccination
	// population information
    // Define vaccine opinion
    // 0 = ok to vaccinate; 1 = do not want to be vaccinated ==> never vaccinate
    vector<int> vaccination_opinion;
    // Nodes (among the detected) for which step vaccination in HH was triggered
    vector<int> vaccination_triggered;
    // List of Workplaces/Schools where the vaccination campaign was done
    set<string> WPS_is_vaccinated;
    // List of Workplaces where the vaccination campaign was not done
    set<string> WP_is_not_vaccinated;
    // List of Schools where the vaccination campaign was not done
    set<string> S_is_not_vaccinated;
    // People who have been vaccinated (0= no vacc; 1= vacc; 9= vacc waiting)
    vector<int> vaccinated;
    // List of Workplaces/Schools where reactive testing must be done
    set<string> place_reactive_vaccination;
    // list of people in the population that can be vaccinated
    vector<int> allNotVaccinated;
    //list of people in WPs that can be vaccinated
//    list<int> WPPeople;
    // list of people in Schools that can be vaccinated
//    list<int> SPeople;
    vector<int> peopleToVaccinateReactive; // for random vaccination
    //
    vector<int> peopleToVaccinateRandom; // for random vaccination
    //random vaccination
    int vacc_total_rand;
    //reactive vaccination
    int vacc_total_react;
    //elderly
    int vacc_init_elderly;
    //adults
    int vacc_init_adults;

} Vaccination;


#define TESTING_MODE 4
#define HOUSEHOLD 0
#define TRACING 1
#define REGULAR_SCREENING 2  // regular
#define REACTIVE_SCREENING 3 // reactive

typedef struct TESTING {
// to hold information related to case detection by testing/tracing
    // Nodes (once detected) for which manual contact tracing was done
    vector<int> ct_done;
    // count those actually found in
    vector<vector<int>> found_in;
    int test_in[TESTING_MODE];
} Testing;


typedef struct SCREENING {
// to hold information related to case detection by screening
    // List of Workplaces/Schools where testing campaign was done
    set<string> WPs_regular_screening;
    // List of Workplaces/Schools where testing campaign was done
    set<string> Ss_regular_screening;
    // List of Workplaces/Schools where testing campaign was done
    set<string> HHs_regular_screening;
    // this is the time at which they start testing  so that it can be regular
    map<string,int> date_regular_screening;
    // List of Workplaces/Schools where reactive testing must be done
    set<string> place_reactive_screening;
    // this is the time at which they were tested so dont repeat too shortly
    map<string,int> date_reactive_screening;
} Screening;


#endif
