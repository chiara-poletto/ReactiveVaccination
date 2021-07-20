/**
    Simulation of an ABM for COVID-19
    @file inputOutput.h
    @author Pierre-Yves Boelle, Chiara Poletto
    @acknowledgment Jesus Moreno for first version
    @version 1.0 2020-04-16
    @license GPL-3.0-or-later
*/
#ifndef INPUTOUTPUT_H
#define INPUTOUTPUT_H

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

using namespace std;

void readFileNames(string fichierInput, Filenames & filenames);

void openSimulationFiles(SimulationFiles & simulationFiles, string &scenario_tag) ;

void closeSimulationFiles(SimulationFiles & simulationFiles);

void readParams(Filenames &filenames,  int lineNumberSought, string & scenario_tag, MapParams &params);

void readPopulation(Filenames &filenames,Population & city,  int*numFirstIndiv,int * numLastIndiv,
		 Params& params) ;

void readNetworkFile(Filenames & filenames, int numFirstIndiv, int numLastIndiv, Networks & contacts, int N_use);

void readHouseholdSchoolWorkplaces(Filenames &filenames,Places&places) ;

void writeAverageFiles(int N_use, Compartments & allEpids, Incidences & incidence, vector<double> R_final, string &scenario_tag, Params & params);

void writeIdxSim(vector<int> idxSim, string scenario_tag) ;

void writeVaccinated(Vaccination & vaccination, Population & city, int real, SimulationFiles & simulationFiles) ;

void writeImmunity(int, Population & city, int);

void printParams(Params & params);

#endif
