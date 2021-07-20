/**
    Simulation of an ABM for COVID-19
    @file inputOutput.cpp
    @author Pierre-Yves Boelle, Chiara Poletto
    @acknowledgment Jesus Moreno for first version
    @version 1.0 2020-04-16
    @license GPL-3.0-or-later
*/

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
#include <sstream>
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

#if defined(_WIN32)
#include <direct.h>
#endif

#include "defs.h"
#include "routine_miscellanea.hpp"

using namespace std;


// --------------------------------------------------------------------------------------------------------------
//  PATH OF INPUT FILES
// --------------------------------------------------------------------------------------------------------------
void readFileNames(string fichierInput, Filenames & filenames)
{
	vector<string> paramNameRequired={"params","age","network","households","workplaces","schools","immunity","susceptibility","subclinical","subclinicalVaccinated"};
	vector<string> paramNameOptional={"universities","invasion"};
	string param;
    ifstream ip;
	ip.open(fichierInput);
	if(!ip.is_open()) std::cout << "ERROR: Could not open " << fichierInput << '\n';

	// read filenames in a map
	map<string,string> fromFile;
	string destination;
	string value;
	string line;
	while (getline(ip,line)) {
    	istringstream sstr(line);
    	// Read in an item
    	getline(sstr, destination, ',');
    	getline(sstr, value);
    	fromFile.insert({destination,value});
    }
	ip.close();

	for (string item: paramNameRequired) {
		if(fromFile.find(item) == fromFile.end()) {
			cout <<"cannot find "<< item << " filename in " << fichierInput <<endl;
			exit(1);
		}
	}
	filenames.paramfile = fromFile["params"];
	cout << "params in " <<filenames.paramfile <<endl;
	filenames.agefile = fromFile["age"];
	cout << "age in " <<filenames.agefile<<endl;
	filenames.HH_list = fromFile["households"];
	cout << "households in " <<filenames.HH_list<<endl;
	filenames.S_list = fromFile["schools"];
	cout << "schools in " <<filenames.S_list<<endl;
	filenames.WP_list = fromFile["workplaces"];
	cout << "workplaces in " <<filenames.WP_list<<endl;
	filenames.immunityfile = fromFile["immunity"];
	cout << "immunity in " <<filenames.immunityfile<<endl;
	filenames.networkfile = fromFile["network"];
	cout << "network in " <<filenames.networkfile<<endl;
	filenames.subclinicalVaccinatedfile = fromFile["subclinicalVaccinated"];
	cout << "subclinicalVaccinated in " <<filenames.subclinicalVaccinatedfile<<endl;
	filenames.subclinicalfile = fromFile["subclinical"];
	cout << "subclinical in " <<filenames.subclinicalfile<<endl;
	filenames.susceptibilityfile = fromFile["susceptibility"];
	cout << "susceptibility in " <<filenames.susceptibilityfile<<endl;

	for ( string item: paramNameOptional) {
		if(fromFile.find(item) == fromFile.end()) {
			cout <<"did not find "<< item << " filename in " << fichierInput <<endl;
		}
	}
	filenames.universityfile = fromFile["universities"];
	cout << "universities in " <<filenames.universityfile<<endl;
	filenames.invasionfile = fromFile["invasion"];
	cout << "invasion in " <<filenames.invasionfile<<endl;

}

/**
 *
 * open simulation files. kept open during the whole sim.
 */
void openSimulationFiles(SimulationFiles & simulationFiles, string &scenario_tag)
{
    // For output files
    string param1 = "_"+ scenario_tag;
    string filename;

	string chemin = "output/" + scenario_tag;
	cout <<"printing results to directory"<<chemin<<endl;
	// Creating scenario_tag directory
	#if defined(_WIN32)
		if (_mkdir(chemin.c_str()) == -1)
	#else
		if (mkdir(chemin.c_str(), 0777) == -1)
	#endif
		cerr << strerror(errno) << endl;
	else
		cout << "Directory "<<chemin <<" created" << endl;

    // Print information about infector and infected
    filename="output/"+scenario_tag+"/infectors"+param1+".csv";
    simulationFiles.infector_file.open(filename);
    (simulationFiles.infector_file) << "infector_node_ID" << "," << "age_infector" << "," << "status_infector" << "," << "vax_infector" << "," << "layer" << "," << "infected_ID" << "," << "age_infected" << ","<< "vax_infected" <<"," <<"ID_place" <<"," << "timestep" << "," << "run" << endl;
   
    // Print information about isolated
    filename = "output/"+scenario_tag+"/isolations"+param1+".csv";
    simulationFiles.isolation_file.open(filename);
    (simulationFiles.isolation_file) << "index_node_ID" << "," << "age_index" << "," << "status" << "," << "sec_isol_ID" << "," << "age_sec" << "," << "status" << "," << "timestep" << "," << "run" << endl;
    
    // Print information about clinical cases
    string symptom_list_name ="output/"+scenario_tag+"/clinical"+param1+".csv";
    simulationFiles.symptom_list_file.open(symptom_list_name);
    (simulationFiles.symptom_list_file) << "index_node_ID" << "," << "age_index" << "," << "vaccination_status" << "," << "timestep" << "," << "run" << endl;
    
    // Print information about subclinical case
    string asymptom_list_name="output/"+scenario_tag+"/subclinical"+param1+".csv";
    simulationFiles.asymptom_list_file.open(asymptom_list_name);
    (simulationFiles.asymptom_list_file) << "index_node_ID" << "," << "age_index" << "," << "vaccination_status" << "," << "timestep" << "," << "run" << endl;
    
    // Print detailed information about vaccination on WP and S
    string outv_detailed_WPS_name="output/"+scenario_tag+"/vaccination_detailed_WPS"+param1+".csv";
    simulationFiles.outv_detailed_WPS_file.open(outv_detailed_WPS_name);
    (simulationFiles.outv_detailed_WPS_file) << "workplace/school ID" << "," << "workplace/school size" << "," << "detected/isolated" << "," << "already vaccinated" << "," << "clinical" << "," << "young" << "," << "vaccinated" << "," << "vacc not suscept" << "," << "vacc suscept" << "," <<  "not willing to vacc" << "," << "timestep" << "," << "run"  << endl;
    
    // Print detailed information about vaccination on HH
    string outv_detailed_HH_name="output/"+scenario_tag+"/vaccination_detailed_HH"+param1+".csv";
    simulationFiles.outv_detailed_HH_file.open(outv_detailed_HH_name);
    (simulationFiles.outv_detailed_HH_file) << "HH vaccinated" << "," << "detected/isolated" << "," << "already vaccinated" << "," << "clinical" << "," << "young" << "," << "vaccinated" << "," << "vacc not suscept" << "," << "vacc suscept" << "," <<  "not willing to vacc" << "," << "timestep" << "," << "run"  << endl;
    
    // Print list of people who got vaccinated
    string out_vaccinated_age_name="output/"+scenario_tag+"/vaccinated_age"+param1+".csv";
    simulationFiles.out_vaccinated_age_file.open(out_vaccinated_age_name);
    (simulationFiles.out_vaccinated_age_file) << "Index ID" << "," << "Age" << "," << "Vaccination status" << "," << "run" << endl;
    
    // Print total number of vaccines used
    string out_vaccine_number_name="output/" + scenario_tag+"/vacc_total"+param1+".csv";
    simulationFiles.vaccine_number_file.open(out_vaccine_number_name);
    simulationFiles.vaccine_number_file << "vacc_rand,vacc_react,vacc_elderly,vacc_adults,vacc_total" << endl;

    // Print places where regular screening occurs
    string out_regular_screening="output/"+scenario_tag+"/regular_screening"+param1+".csv";
    simulationFiles.out_regular_screening_file.open(out_regular_screening);

    // Print places where regular screening occurs
    string out_reactive_screening="output/"+scenario_tag+"/reactive_screening"+param1+".csv";
    simulationFiles.out_reactive_screening_file.open(out_reactive_screening);

    // Print whatever
    string out_miscellaneous="output/"+scenario_tag+"/miscellaneous"+param1+".csv";
    simulationFiles.miscellaneous_file.open(out_miscellaneous);

    // Print extinction
    string out_extinction="output/"+scenario_tag+"/extinction"+param1+".csv";
    simulationFiles.extinction_file.open(out_extinction);
    simulationFiles.extinction_file << "real,nSim,t_first,status_first,t_thr,status_thr,t_extinct,status_extinct" <<endl;

    //Print reactive vaccination numbers
    string out_reactive_vaccination="output/"+scenario_tag+"/reactive_vaccination"+param1+".csv";
    simulationFiles.reactive_vaccination_file.open(out_reactive_vaccination);

}

/*
 * close simulation files before exit
 */
void closeSimulationFiles(SimulationFiles & simulationFiles) {
	simulationFiles.infector_file.close();
	    // To print a file with details of isolations
	simulationFiles.isolation_file.close();
	    // To print a file with time of symptom onset and other informations for clinical and subclinical nodes
	simulationFiles.symptom_list_file.close();
	simulationFiles.asymptom_list_file.close();
	    // TO print detailed information on vaccination
	simulationFiles.outv_detailed_WPS_file.close();
	simulationFiles.outv_detailed_HH_file.close();
		// To print detailed information on the age of the vaccinated
	simulationFiles.out_vaccinated_age_file.close();
		// To print number of vaccines used
	simulationFiles.vaccine_number_file.close();

    simulationFiles.out_regular_screening_file.close();

    simulationFiles.out_reactive_screening_file.close();

    simulationFiles.miscellaneous_file.close();
    simulationFiles.extinction_file.close();
    simulationFiles.reactive_vaccination_file.close();

}//

/*
 * read param files, make map
 */
void readParams( Filenames & filenames,   int lineNumberSought, string & scenario_tag, MapParams &params) {
	//read parameters from designated file
	vector<string> list_params_names;
	vector<string> list_values;
	string line, csvItem, csvName;
	ifstream myfile;
	int lineNumber = 0;
	int lineNumberName = 1;

	myfile.open(filenames.paramfile);

	if(myfile.is_open()) {
		while (getline(myfile,line)) {
			lineNumber++;
			if(lineNumber == lineNumberName) {
				istringstream myline(line);
				while(getline(myline, csvName, ',')) {
					//cout << csvName << endl;
					list_params_names.push_back(csvName);
				}
			}
			if(lineNumber == lineNumberSought) {
				//cout << line << endl; ;
				istringstream myline(line);
				while(getline(myline, csvItem, ',')) {
					//cout << csvItem << endl;
					list_values.push_back(csvItem);
				}
			}
		}
		myfile.close();
	} else {
		cout << "could not read params in " << filenames.paramfile << "\n";
	}

	scenario_tag = list_values[0];
	cout << scenario_tag << endl;

	vector<double> list_params_values;
	for(int i=1; i<(int) list_params_names.size(); ++i){
		//cout << "Names = " << list_params_names[i] << endl;
		//cout << "Values = " << list_values[i] << endl;
		list_params_values.push_back(stod(list_values[i])); //THIS LINE DOES NOT WORK
	}
	list_params_names.erase (list_params_names.begin());

	//Map creation contaning parameters names and values
	for (size_t i = 0; i < list_params_values.size(); ++i)
		(params)[list_params_names[i]] = list_params_values[i];

//	for (const auto& x :  params) {
//		cout << x.first << ": " << x.second << "\n";
//	}

}



/*
 * read immunity file. starts with V2
 */
void readImmunity(Filenames & filenames, vector<int> & immune, vector<vector<int>> & exposed, Params & params) {
    ifstream ip4;

    if (params.immunity_option!=0)
    {
    	cout << "Immunity level " << params.immunity_option << endl;
    	int im_option = (int)params.immunity_option;
    	string fullFilename = filenames.immunityfile+to_string(im_option)+".csv";
    	cout << "reading Immunity from " << fullFilename<<endl;
        ip4.open(fullFilename);

        if(!ip4.is_open()) std::cout << "ERROR: could not read from" << fullFilename <<'\n';

        string index, status1, status2;
        getline(ip4,index,',');
        getline(ip4,status1,',');
        getline(ip4,status2,'\n');

        if (index != "V2") {
        	cout << "Old version of immunity files\nRun vaccination.exe -i to regenerate immunity files\n";
        	exit(1);
        }
        int n=0;
        while(ip4.good())
        {
            getline(ip4,index,',');
            getline(ip4,status1,',');
            getline(ip4,status2,'\n');
            immune.push_back(stoi(status1));
            if(status2 != "0") {
            	exposed[0].push_back(n); // vector of exposed
            	exposed[1].push_back(stoi(index)); // vector of exposed
            	exposed[2].push_back(stoi(status2)); // vector of exposed
            	n++;
            }
        }
        ip4.close();
    }
}

/*
 * read invasion file. file of index of cases
 */
void readInvasion(Filenames & filenames, vector<vector<int>> & invasion) {
    ifstream ip4;

    string fullFilename = filenames.invasionfile;
    cout << "reading invasion from " << fullFilename<<endl;
    ip4.open(fullFilename);

    if(!ip4.is_open()) std::cout << "ERROR: could not read from" << fullFilename <<'\n';

    //file is a list of indices in the order of apparition
    int n=0;
    string line;
    while (getline(ip4,line)) {
    	string item;
    	istringstream sstr(line);
    	// Read in an item
    	while (getline(sstr, item, ','))
    	{
    		invasion[n].push_back(stoi(item));
    	}
    	n++;
    }
    ip4.close();
}


/*
 * read susceptibility, asymptomatic, then age
 */
void readPopulation(Filenames & filenames, Population & city,  int*numFirstIndiv,int * numLastIndiv,
		 Params& params) {

	ifstream ip3;

	// format is age.lo-age.high susceptibility
	string agelo;
	string ageup;
	string value;

	// read susceptibility file
	//	vector<double> susceptibility={0.4,0.38,0.79,0.86,0.8,0.82,0.88,0.74,0.74,0.74,0.74,0.74,0.74};
	vector<double> susceptibility; // one value per age
	ip3.open(filenames.susceptibilityfile);
	if(!ip3.is_open()) std::cout << "ERROR: could not read from " << filenames.susceptibilityfile  << '\n';
	while(ip3.good()) {
		getline(ip3,agelo,',');
		getline(ip3,ageup,',');
		getline(ip3,value,'\n');
		if (agelo=="") break;
		for (int age=stoi(agelo); age <=stoi(ageup); age++)
			susceptibility.push_back(stof(value));
	}
	ip3.close();

	// read subclinical file
	//	vector<double> propSubclinical={0.29,0.21,0.27,0.33,0.4,0.49,0.63,0.69,0.69,0.69,0.69,0.69,0.69#};
	vector<double> propSubclinical; // one value per age
	ip3.open(filenames.subclinicalfile);
	if(!ip3.is_open()) std::cout << "ERROR: could not read from " << filenames.subclinicalfile  << '\n';
	// format is age.lo-age.high susceptibility
	while(ip3.good()) {
		getline(ip3,agelo,',');
		getline(ip3,ageup,',');
		getline(ip3,value,'\n');
		if (agelo=="") break;
		for (int age=stoi(agelo); age <=stoi(ageup); age++)
			propSubclinical.push_back(stof(value));
	}
	ip3.close();

	// read subclinicalVaccinated file
	//	vector<double> propSubclinicalVaccinated={0.29,0.21,0.27,0.33,0.4,0.49,0.63,0.69,0.69,0.69,0.69,0.69,0.69};
	vector<double> propSubclinicalVaccinated; // one value per age
	ip3.open(filenames.subclinicalVaccinatedfile);
	if(!ip3.is_open()) std::cout << "ERROR: could not read from " << filenames.subclinicalVaccinatedfile  << '\n';
	// format is age.lo-age.high susceptibility
	while(ip3.good()) {
		getline(ip3,agelo,',');
		getline(ip3,ageup,',');
		getline(ip3,value,'\n');
		if (agelo=="") break;
		for (int age=stoi(agelo); age <=stoi(ageup); age++)
			propSubclinicalVaccinated.push_back(stof(value));
	}
	ip3.close();


		ip3.open(filenames.agefile);

	if(!ip3.is_open()) std::cout << "ERROR: could not read from " << filenames.agefile  << '\n';

	string indv4, age1;
	vector<int> indv;
	vector<string> comm_name;

// 	intialise here
	for (int i=0; i< (int) propSubclinical.size(); i++) {
		propSubclinical[i] = 1.0-propSubclinical[i];
		propSubclinicalVaccinated[i] = 1.0- propSubclinicalVaccinated[i]*params.v_symp;
	}

	city.N_use=0;
	// For the header
	getline(ip3,indv4,',');
	getline(ip3,age1,'\n');
	while(ip3.good()) {
		getline(ip3,indv4,',');
		getline(ip3,age1,'\n');
		indv.push_back(stoi(indv4)-1);
		int ageN=stoi(age1);
		city.Age.push_back(ageN); //read age
		city.age_w.push_back(susceptibility[ageN]); // susceptibility
		city.a_propor.push_back(propSubclinical[ageN]); // asymptomatic
		city.va_propor.push_back(propSubclinicalVaccinated[ageN]);  //asympto vaccinated
		++city.N_use;
	}
	ip3.close();

	if (params.computeImmunity==0) {
		city.exposed_file.resize(3);
		readImmunity(filenames, city.status_file, city.exposed_file, params);
	} else {
		cout << "immunity set to 0" << endl;
	}

	// Print the number of individuals in simulation
	cout << "N_use " << city.N_use << endl;
	// Displays file ID of first individual used in simulation
	cout << "indv[0]" << indv[0] << endl;
	//add 1 because count start from 1
	*numFirstIndiv=indv[0]+1;
	*numLastIndiv=indv[indv.size()-1]+1;
}

/*
 * read network files
 */
void readNetworkFile(Filenames & filenames, int numFirstIndiv, int numLastIndiv, Networks & contacts, int N_use) {
	// files are red once and kept
    ifstream ip;
    ip.open(filenames.networkfile);
    if(!ip.is_open()) std::cout << "ERROR: could not read from " << filenames.networkfile << '\n';

    string indv1, indv2, place1, freqnum1;

    //compute number of contacts by people
    contacts.maxContacts.resize(N_use,vector<int>(5));
    contacts.maxContacts.assign(N_use,vector<int>(5,0));

    for (int i =0; i<5; i++) contacts.allContactsInLayers[i] = 0;
    //header
    getline(ip,indv1,',');
    getline(ip,indv2,',');
    getline(ip,place1,',');
    getline(ip,freqnum1,'\n');
    while(ip.good())
    {
        getline(ip,indv1,',');
        getline(ip,indv2,',');
        getline(ip,place1,',');
        getline(ip,freqnum1,'\n');
        if (stoi(indv1)<=numLastIndiv && stoi(indv1)>=numFirstIndiv && stoi(indv2)<=numLastIndiv && stoi(indv2)>=numFirstIndiv)
        {
        	int nIdv1=stoi(indv1)-numFirstIndiv;
        	int nIdv2=stoi(indv2)-numFirstIndiv;
        	int place=stoi(place1);
            contacts.node1.push_back(nIdv1);
            contacts.node2.push_back(nIdv2);
            contacts.place.push_back(place);
            contacts.freq.push_back(stod(freqnum1));
            contacts.maxContacts[nIdv1][place]+=1;
            contacts.maxContacts[nIdv2][place]+=1;
            contacts.allContactsInLayers[place]++;
        }
    }
    ip.close();
    cout<<"all Contacts : HH " << contacts.allContactsInLayers[0] << " WP " << contacts.allContactsInLayers[1] <<
    		" S " <<contacts.allContactsInLayers[2] <<
    		" CO " << contacts.allContactsInLayers[3] << " TR " << contacts.allContactsInLayers[4] <<endl ;
    //	 look for maximum Contacts
    contacts.maxOfMaxContacts = 0;
    for (int i=0; i < N_use; i++) {
    	for (int m=0; m<5; m++) {
    		if (contacts.maxContacts[i][m]>contacts.maxOfMaxContacts) contacts.maxOfMaxContacts =  contacts.maxContacts[i][m];
    	}
    }
    cout << "max contacts :" << contacts.maxOfMaxContacts <<endl;
}



/*
 * read details of places
 */
void readPlace(string &fileName, map<int, string> & node_place, map<string, vector<int>> &place_nodes, set<string>&places) {

	string indv5, unused, idPlace, wps;
	ifstream ip;
	ip.open(fileName);
	//
	if(!ip.is_open()) std::cout << "ERROR: unable to open file " << fileName << " \n";
	getline(ip,indv5,',');
	getline(ip,unused,',');
	getline(ip,idPlace,'\n');
	while(ip.good()) {
		getline(ip,indv5,',');
		getline(ip,unused,',');
		getline(ip,idPlace,'\n');
		//cout << indv5 << ", " << hh << endl;
		if(indv5!=""){
			node_place[stoi(indv5)-1]= idPlace;
			place_nodes[idPlace].push_back(stoi(indv5)-1);
			places.insert(idPlace);
		}
	}
	ip.close();
}

/*
 * read WPs/Ss
 * change universities to WPs
 */
void readHouseholdSchoolWorkplaces(Filenames & filenames, Places&places) {
	readPlace(filenames.HH_list, places.node_HH, places.HH_nodes,places.HHs);
	readPlace(filenames.S_list, places.node_S, places.S_nodes, places.Ss);
	readPlace(filenames.WP_list, places.node_WP, places.WP_nodes, places.WPs);

	// if required change universities fo workplaces
	if (filenames.universityfile!="") {
		// read in university names
		ifstream ip;
		ip.open(filenames.universityfile);
		if(!ip.is_open()) std::cout << "ERROR: unable to open " << filenames.universityfile << " \n";
		string idUniv;
		string tmp1;
		string tmp2;
		while(ip.good()) {
			getline(ip,tmp1,',');
			getline(ip,tmp2,',');
			getline(ip,idUniv,'\n');
			places.universities.insert(idUniv);
		}
		cout <<"read "<<places.universities.size()<< "universities"<<endl;
		for (string univ: places.universities) {
			while (places.S_nodes[univ].size()>0) {
				places.WP_nodes[univ].push_back(places.S_nodes[univ].back());
				places.S_nodes[univ].pop_back();
			}
			places.Ss.erase(univ);
			places.WPs.insert(univ);
			places.S_nodes.erase(univ);
		}
	}
}

/*
 * write average files
 */
void writeAverageFiles(int N_use, Compartments & allEpids, Incidences & incidence, vector<double> R_final, string &scenario_tag, Params & params) {

    //PRINT OUT STATISTICS

	vector<double> S_av; S_av.resize(params.max_steps);
	vector<double> E_av; E_av.resize(params.max_steps);
	vector<double> P1_av; P1_av.resize(params.max_steps);
	vector<double> P2_av; P2_av.resize(params.max_steps);
	vector<double> SI_av; SI_av.resize(params.max_steps);
	vector<double> AI_av; AI_av.resize(params.max_steps);
	vector<double> R_av; R_av.resize(params.max_steps);
	vector<double> wvS_av; wvS_av.resize(params.max_steps);
	vector<double> vS_av; vS_av.resize(params.max_steps);
	vector<double> vE_av; vE_av.resize(params.max_steps);
	vector<double> vP1_av; vP1_av.resize(params.max_steps);
	vector<double> vP2_av; vP2_av.resize(params.max_steps);
	vector<double> vSI_av; vSI_av.resize(params.max_steps);
	vector<double> vAI_av; vAI_av.resize(params.max_steps);
	vector<double> vR_av; vR_av.resize(params.max_steps);
	vector<double> isol_av; isol_av.resize(params.max_steps);
	vector<double> isol_and_infected_av; isol_and_infected_av.resize(params.max_steps);
	vector<double> vaccinated_FUL_av; vaccinated_FUL_av.resize(params.max_steps);
	vector<double> vaccinated_MID_av; vaccinated_MID_av.resize(params.max_steps);
	vector<double> vaccinated_LOW_av; vaccinated_LOW_av.resize(params.max_steps);
	vector<double> vaccines_used_av; vaccines_used_av.resize(params.max_steps);
	vector<double> vaccines_used_random_av; vaccines_used_random_av.resize(params.max_steps);
	vector<double> vaccines_used_reactive_av; vaccines_used_reactive_av.resize(params.max_steps);
	vector<double> detected_av; detected_av.resize(params.max_steps);
	vector<double> found_by_ct_av; found_by_ct_av.resize(params.max_steps);
	vector<double> found_by_reg_av; found_by_reg_av.resize(params.max_steps);
	vector<double> found_by_reac_av; found_by_reac_av.resize(params.max_steps);
	vector<double> WP_infected; WP_infected.resize(params.max_steps);
	vector<double> S_infected; S_infected.resize(params.max_steps);
	string param1 = "_"+ scenario_tag;

    // Print out averages of: fraction of pop in each compartment, isolated, isolated and infected
    for (int i{0} ; i<params.max_steps ; ++i)
    {
        S_av[i] = allEpids.S[i]/params.N_real;
        E_av[i] = allEpids.E[i]/params.N_real;
        P1_av[i] = allEpids.P1[i]/params.N_real;
        P2_av[i] = allEpids.P2[i]/params.N_real;
        SI_av[i] = allEpids.SI[i]/params.N_real;
        AI_av[i] = allEpids.AI[i]/params.N_real;
        R_av[i] = allEpids.R[i]/params.N_real;
        wvS_av[i] = allEpids.wvS[i]/params.N_real;
        vS_av[i] = allEpids.vS[i]/params.N_real;
        vE_av[i] = allEpids.vE[i]/params.N_real;
        vP1_av[i] = allEpids.vP1[i]/params.N_real;
        vP2_av[i] = allEpids.vP2[i]/params.N_real;
        vSI_av[i] = allEpids.vSI[i]/params.N_real;
        vAI_av[i] = allEpids.vAI[i]/params.N_real;
        vR_av[i] = allEpids.vR[i]/params.N_real;

        WP_infected[i] = allEpids.WP_infected[i]/params.N_real;
        S_infected[i] = allEpids.S_infected[i]/params.N_real;

        isol_av[i] = allEpids.isolated[i]/params.N_real;
        isol_and_infected_av[i] = allEpids.infected_and_isolated[i]/params.N_real;

        vaccinated_FUL_av[i] = allEpids.vaccinated_FUL[i]/params.N_real;
        vaccinated_MID_av[i] = allEpids.vaccinated_MID[i]/params.N_real;
        vaccinated_LOW_av[i] = allEpids.vaccinated_LOW[i]/params.N_real;
        vaccines_used_av[i] = allEpids.vaccines_used[i]/params.N_real;
        vaccines_used_random_av[i] = allEpids.vaccines_used_random[i]/params.N_real;
        vaccines_used_reactive_av[i] = allEpids.vaccines_used_reactive[i]/params.N_real;

        detected_av[i] = allEpids.detected[i]/params.N_real;
        found_by_ct_av[i] = allEpids.found_by_ct[i]/params.N_real;
        found_by_reg_av[i] = allEpids.found_by_regular_screening[i]/params.N_real;
        found_by_reac_av[i] = allEpids.found_by_reactive_screening[i]/params.N_real;
    }

    string final_name="output/"+scenario_tag+"/average_SIR_series" + param1;
    save_csv_double_all(final_name,S_av,vS_av,wvS_av,E_av,vE_av,P1_av,vP1_av,SI_av,vSI_av,P2_av,vP2_av,AI_av,vAI_av,R_av,vR_av,N_use);
    string av_percent_isol_st="output/"+scenario_tag+"/av_percent_isolated_series"+param1;
    save_csv_double_vec(av_percent_isol_st,isol_av);
    string inf_a_iso_st="output/"+scenario_tag+"/infected_and_isolated"+param1;
    save_csv_double_vec(inf_a_iso_st,isol_and_infected_av);

    string place_infected="output/"+scenario_tag+"/place_infected"+param1;
    string colnames="WP,S";
    save_csv_double_vec2(place_infected,colnames,WP_infected,S_infected);
    string percent_detected="output/"+scenario_tag+"/detected"+param1;
    save_csv_double_vec(percent_detected,detected_av);
    string percent_vaccinated="output/"+scenario_tag+"/vaccinated"+param1;
    string v1="FUL";
    string v2="MID";
    string v3="LOW";

    save_csv_double_vec3(percent_vaccinated,v1,v2,v3, vaccinated_FUL_av, vaccinated_MID_av, vaccinated_LOW_av);
    string vaccines_used="output/"+scenario_tag+"/vaccines_used"+param1;
    v1="all";
    v2="reactive";
    v3="random";

    save_csv_double_vec3(vaccines_used,v1,v2,v3,vaccines_used_av,vaccines_used_reactive_av,vaccines_used_random_av);



    string percent_found_by_ct="output/"+scenario_tag+"/found_by_ct"+param1;
    save_csv_double_vec(percent_found_by_ct,found_by_ct_av);
    string percent_found_by_reg_s="output/"+scenario_tag+"/found_by_regular_scr"+param1;
    save_csv_double_vec(percent_found_by_reg_s,found_by_reg_av);
    string percent_found_by_reac_s="output/"+scenario_tag+"/found_by_reactive_scr"+param1;
    save_csv_double_vec(percent_found_by_reac_s,found_by_reac_av);

    // Print out final number of Recovered for each run
    string final_R="output/"+scenario_tag+"/R_final" + param1;
    save_csv_double_vec(final_R, R_final);


    // PRINT OUT INCIDENCES RUN BY RUN
    string name ="output/"+scenario_tag+"/incP1"+param1+".csv";
    ofstream csvfile(name);
    if (csvfile.is_open())
    {   csvfile << "time";

        for (int n=0; n<params.N_real ; n++)
        {   string real_head = to_string(n);
            csvfile << "," << real_head ;}
        csvfile << endl;

        for (int j=0; j<params.max_steps ; j++)
        {   csvfile << j;
            for (int i=0; i<params.N_real ; i++)
            {csvfile << "," << incidence.P1[i][j];}
            csvfile << endl; }
    	csvfile.close();
    } else
    {cout << "File "<<name <<" could not be opened";}

    string name2 ="output/"+scenario_tag+"/incP2"+param1+".csv";
    ofstream csvfile2(name2);
    if (csvfile2.is_open())
    {   csvfile2 << "time";

        for (int n=0; n<params.N_real ; n++)
        {   string real_head = to_string(n);
            csvfile2 << "," << real_head ;}
        csvfile2 << endl;

        for (int j=0; j<params.max_steps ; j++)
        {   csvfile2 << j;
            for (int i=0; i<params.N_real ; i++)
            {csvfile2 << "," << incidence.P2[i][j];}
            csvfile2 << endl; }
    	csvfile2.close();
    }
    else
    {cout << "File "<< name2 <<" could not be opened";}

    string name3 ="output/"+scenario_tag+"/incSI"+param1+".csv";
    ofstream csvfile3(name3);
    if (csvfile3.is_open())
    {   csvfile3 << "time";

        for (int n=0; n<params.N_real ; n++)
        {   string real_head = to_string(n);
            csvfile3 << "," << real_head ;}
        csvfile3 << endl;

        for (int j=0; j<params.max_steps ; j++)
        {   csvfile3 << j;
            for (int i=0; i<params.N_real ; i++)
            {csvfile3 << "," << incidence.SI[i][j];}
            csvfile3 << endl; }
    	csvfile3.close();
    }
    else
    {cout << "File "<< name3<<" could not be opened";}

    string name4 ="output/"+scenario_tag+"/incAI"+param1+".csv";
    ofstream csvfile4(name4);
    if (csvfile4.is_open())
    {   csvfile4 << "time";

        for (int n=0; n<params.N_real ; n++)
        {   string real_head = to_string(n);
            csvfile4 << "," << real_head ;}
        csvfile4 << endl;

        for (int j=0; j<params.max_steps ; j++)
        {   csvfile4 << j;
            for (int i=0; i<params.N_real ; i++)
            {csvfile4 << "," << incidence.AI[i][j];}
            csvfile4 << endl; }
    	csvfile4.close();
    }
    else
    {cout << "File "<< name4<<" could not be opened";}

    string namev ="output/"+scenario_tag+"/incvP1"+param1+".csv";
    ofstream csvfilev(namev);
    if (csvfilev.is_open())
    {   csvfilev << "time";

        for (int n=0; n<params.N_real ; n++)
        {   string real_head = to_string(n);
            csvfilev << "," << real_head ;}
        csvfilev << endl;

        for (int j=0; j<params.max_steps ; j++)
        {   csvfilev << j;
            for (int i=0; i<params.N_real ; i++)
            {csvfilev << "," << incidence.vP1[i][j];}
            csvfilev << endl; }
    	csvfilev.close();
    }
    else
    {cout << "File "<< namev<<" could not be opened";}

    string name2v="output/"+scenario_tag+"/incvP2"+param1+".csv";
    ofstream csvfile2v(name2v);
    if (csvfile2v.is_open())
    {   csvfile2v << "time";

        for (int n=0; n<params.N_real ; n++)
        {   string real_head = to_string(n);
            csvfile2v << "," << real_head ;}
        csvfile2v << endl;

        for (int j=0; j<params.max_steps ; j++)
        {   csvfile2v << j;
            for (int i=0; i<params.N_real ; i++)
            {csvfile2v << "," << incidence.vP2[i][j];}
            csvfile2v << endl; }
    	csvfile2v.close();
    }
    else
    {cout << "File "<<name2v <<" could not be opened";}

    string name3v ="output/"+scenario_tag+"/incvSI"+param1+".csv";
    ofstream csvfile3v(name3v);
    if (csvfile3v.is_open())
    {   csvfile3v << "time";

        for (int n=0; n<params.N_real ; n++)
        {   string real_head = to_string(n);
            csvfile3v << "," << real_head ;}
        csvfile3v << endl;

        for (int j=0; j<params.max_steps ; j++)
        {   csvfile3v << j;
            for (int i=0; i<params.N_real ; i++)
            {csvfile3v << "," << incidence.vSI[i][j];}
            csvfile3v << endl; }
    	csvfile3v.close();
    }
    else
    {cout << "File "<<name3v <<" could not be opened";}

    string name4v ="output/"+scenario_tag+"/incvAI"+param1+".csv";
    ofstream csvfile4v(name4v);
    if (csvfile4v.is_open())
    {   csvfile4v << "time";

        for (int n=0; n<params.N_real ; n++)
        {   string real_head = to_string(n);
            csvfile4v << "," << real_head ;}
        csvfile4v << endl;

        for (int j=0; j<params.max_steps ; j++)
        {   csvfile4v << j;
            for (int i=0; i<params.N_real ; i++)
            {csvfile4v << "," << incidence.vAI[i][j];}
            csvfile4v << endl; }
    	csvfile4v.close();
    }
    else
    {cout << "File "<<name4v <<" not be opened";}

}

void writeVaccinated(Vaccination & vaccination, Population & city, int real, SimulationFiles & simulationFiles) {

	simulationFiles.vaccine_number_file << vaccination.vacc_total_rand <<"," <<vaccination.vacc_total_react <<"," <<vaccination.vacc_init_elderly <<"," <<vaccination.vacc_init_adults <<"," << vaccination.vacc_total_rand + vaccination.vacc_total_react + vaccination.vacc_init_elderly + vaccination.vacc_init_adults  << endl;

	for (int i{0}; i< (int) vaccination.vaccinated.size(); ++i){
        if (vaccination.vaccinated[i]!=VAX_NOT){
            simulationFiles.out_vaccinated_age_file << i << "," << city.Age[i] << "," << vaccination.vaccinated[i] << "," << real << endl;
        }
    }
}

void writeIdxSim(vector<int> idxSim, string scenario_tag) {
    string final_name="output/"+scenario_tag+"/idxSim";
    save_csv_int_vec(final_name, idxSim);
}

void writeImmunity(int t, Population & city, int prev) {
// write the status of persons for a range of values of immunity in the 0 - 30 range.
	// Creating scenario_tag directory
	string chemin;
	chemin ="output/immunity";

	#if defined(_WIN32)
		if (_mkdir(chemin.c_str()) == -1)
	#else
		if (mkdir(chemin.c_str(), 0777) == -1)
	#endif
			cout<<"created "<< chemin;

	string filename;
	filename= "output/immunity/status_file_" + to_string(prev) + ".csv";
	ofstream outfile;
	outfile.open(filename);
	outfile << "V2,immune,infected - t=" << t <<endl;
	for (int i=0; i<city.N_use; i++) {
		outfile << i;
		if (city.status[i]==0)
			outfile << ",0,0" <<endl;
		else if(city.status[i]==2)
			outfile << ",2,0" <<endl;
		else if (city.status[i]==10)
			outfile << ",0,"<<city.status[i] <<endl;
		else
			outfile << ",2,"<<city.status[i] <<endl;
	}
	outfile.close();
	cout << "wrote "<<filename<< " t= " << t << "\n";
}

/*
 * print parameters for check
 */
void printParams(Params & params) {
	cout <<"Parameters after conversion:"<<endl;
	cout <<"acquai_frequency: " << params.acquai_frequency << endl;
	cout <<"age_min_testing: " << params.age_min_testing << endl;
	cout <<"age_th: " << params.age_th << endl;
	cout <<"age_vax_elderly_th: " << params.age_vax_elderly_th << endl;
	cout <<"beta: " << params.beta << endl;
	cout <<"computeImmunity: " << params.computeImmunity << endl;
	cout <<"coverage_vacc_higher_65: " << params.coverage_vacc_higher_65 << endl;
	cout <<"coverage_vacc_lower_65: " << params.coverage_vacc_lower_65 << endl;
	cout <<"dur_isol_inf: " << params.dur_isol_inf << endl;
	cout <<"dur_isol_rec: " << params.dur_isol_rec << endl;
	cout <<"dur_isol_sus: " << params.dur_isol_sus << endl;
	cout <<"e_rate: " << params.e_rate << endl;
	cout <<"E_ve_rate: " << params.E_ve_rate << endl;
	cout <<"f_HH_screening: " << params.f_HH_screening << endl;
	cout <<"f_S_screening: " << params.f_S_screening << endl;
	cout <<"f_WP_screening: " << params.f_WP_screening << endl;
	cout <<"factor_isol: " << params.factor_isol << endl;
	cout <<"frac_elderly_vacc: " << params.frac_elderly_vacc << endl;
	cout <<"frac_init_vacc: " << params.frac_init_vacc << endl;
	cout <<"gamm: " << params.gamm << endl;
	cout <<"HHisolation_inplace: " << params.HHisolation_inplace << endl;
	cout <<"immunity_option: " << params.immunity_option << endl;
	cout <<"immunity_random: " << params.immunity_random << endl;
	cout <<"init_incid_100000: " << params.init_incid_100000 << endl;
	cout <<"interv_cumul_cases: " << params.interv_cumul_cases << endl;
	cout <<"interval_screening: " << params.interval_screening << endl;
	cout <<"manualCT_inplace: " << params.manualCT_inplace << endl;
	cout <<"max_steps: " << params.max_steps << endl;
	cout <<"N_real: " << params.N_real << endl;
	cout <<"n_start_reactive: " << params.n_start_reactive << endl;
	cout <<"number_networks: " << params.number_networks << endl;
	cout <<"p_c_HH: " << params.p_c_HH << endl;
	cout <<"p_rA: " << params.p_rA << endl;
	cout <<"p_rS: " << params.p_rS << endl;
	cout <<"p1_rate: " << params.p1_rate << endl;
	cout <<"p2_rate: " << params.p2_rate << endl;
	cout <<"pct_rand_init_immunity: " << params.pct_rand_init_immunity << endl;
	cout <<"pct_teleworking: " << params.pct_teleworking << endl;
	cout <<"pct_comty_contact: " << params.pct_comty_contact << endl;
	cout <<"print_infectors: " << params.print_infectors << endl;
	cout <<"print_isolations: " << params.print_isolations << endl;
	cout <<"r_a: " << params.r_a << endl;
	cout <<"r_d_asymp: " << params.r_d_asymp << endl;
	cout <<"r_d_symp: " << params.r_d_symp << endl;
	cout <<"r_mCT: " << params.r_mCT << endl;
	cout <<"r_react_screen: " << params.r_react_screen << endl;
	cout <<"r_V: " << params.r_V << endl;
	cout <<"random_type: " << params.random_type << endl;
	cout <<"react_type: " << params.react_type << endl;
	cout <<"reactive_screening_S_inplace: " << params.reactive_screening_S_inplace << endl;
	cout <<"reactive_screening_WP_inplace: " << params.reactive_screening_WP_inplace << endl;
	cout <<"reactive_tests_daily: " << params.reactive_tests_daily << endl;
	cout <<"reactive_tests_limit: " << params.reactive_tests_limit << endl;
	cout <<"reactive_vacc_inplace: " << params.reactive_vacc_inplace << endl;
	cout <<"refresh_network: " << params.refresh_network << endl;
	cout <<"regular_screening_HH_inplace: " << params.regular_screening_HH_inplace << endl;
	cout <<"regular_screening_S_inplace: " << params.regular_screening_S_inplace << endl;
	cout <<"regular_screening_WP_inplace: " << params.regular_screening_WP_inplace << endl;
	cout <<"screening_min_size: " << params.screening_min_size << endl;
	cout <<"skip_iso: " << params.skip_iso << endl;
	cout <<"testing_proportion_HH: " << params.testing_proportion_HH << endl;
	cout <<"testing_proportion_S: " << params.testing_proportion_S << endl;
	cout <<"testing_proportion_WP: " << params.testing_proportion_WP << endl;
	cout <<"th_cluster_screen: " << params.th_cluster_screen << endl;
	cout <<"th_cluster_vacc: " << params.th_cluster_vacc << endl;
	cout <<"tracing_lookback: " << params.tracing_lookback << endl;
	cout <<"typeSeeding: " << params.typeSeeding << endl;
	cout <<"v_sus_FUL: " << params.v_sus_FUL << endl;
	cout <<"v_sus_MID: " << params.v_sus_MID << endl;
	cout <<"v_symp: " << params.v_symp << endl;
	cout <<"v_trans: " << params.v_trans << endl;
	cout <<"vac_daily_randl_pop: " << params.vac_daily_randl_pop << endl;
	cout <<"vac_elderly_inplace: " << params.vac_elderly_inplace << endl;
	cout <<"vac_init_inplace: " << params.vac_init_inplace << endl;
	cout <<"vac_rand_inplace: " << params.vac_rand_inplace << endl;
	cout <<"vacc_daily_react: " << params.vacc_daily_react << endl;
	cout <<"vacc_limit_rand_pop: " << params.vacc_limit_rand_pop << endl;
	cout <<"vacc_limit_react: " << params.vacc_limit_react << endl;
	cout <<"vacc_min_WP_size: " << params.vacc_min_WP_size << endl;
	cout <<"ve_rate: " << params.ve_rate << endl;
	cout <<"vgamm: " << params.vgamm << endl;
	cout <<"vp1_rate: " << params.vp1_rate << endl;
	cout <<"vp2_rate: " << params.vp2_rate << endl;
	cout <<"w_rate_LOW_to_MID: " << params.w_rate_LOW_to_MID << endl;
	cout <<"w_rate_MID_to_FUL: " << params.w_rate_MID_to_FUL << endl;
	cout <<"wgt_lyr_household: " << params.wgt_lyr_household << endl;
	cout <<"wgt_lyr_workspace: " << params.wgt_lyr_workspace << endl;
	cout <<"wgt_lyr_school: " << params.wgt_lyr_school << endl;
	cout <<"wgt_lyr_community: " << params.wgt_lyr_community << endl;
	cout <<"wgt_lyr_transport: " << params.wgt_lyr_transport << endl;
	cout <<"will_be_d_asymp: " << params.will_be_d_asymp << endl;
	cout <<"will_be_d_symp: " << params.will_be_d_symp << endl;
	cout << "after Invasion:"<<endl;
	cout << "p_c_HH_aft_inv: " << params.p_c_HH_aft_inv << endl;
	cout << "p_rA_aft_inv: " << params.p_rA_aft_inv << endl;
	cout << "p_rS_aft_inv: " <<  params.p_rS_aft_inv << endl;
	cout << "r_d_asymp_aft_inv: " << params.r_d_asymp_aft_inv << endl;
	cout << "r_d_symp_aft_inv: " << params.r_d_symp_aft_inv << endl;
	cout << "will_be_d_asymp_aft_inv: " <<  params.will_be_d_asymp_aft_inv << endl;
	cout << "will_be_d_symp_aft_inv: "  << params.will_be_d_symp_aft_inv << endl;
	cout << "skip_iso_aft_inv: " << params.skip_iso_aft_inv << endl;
	cout << "incrVaccAccept_inplace: " << params.incrVaccAccept_inplace<<endl;
	cout << "doIncreasedVaccination: " << params.doIncreasedVaccination<<endl;
	cout <<"End of parameters"<<endl;

}
