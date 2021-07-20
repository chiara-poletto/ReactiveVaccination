/**
    Simulation of an ABM for COVID-19
    @file init.h
    @author Pierre-Yves Boelle, Chiara Poletto
    @version 1.0 2020-04-16
    @license GPL-3.0-or-later
*/

#ifndef INIT_H
#define INIT_H

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

using namespace std;

#include "defs.h"

void initImprovedCountInCases(Population & city);

void initRealisation_Status(Population & city, Params &params,default_random_engine& generator) ;

void initEpid(Compartments &epid, Params &params) ;

void initIncidence(Incidences & incidence, Params &params);

void initImmunity(Population & city, Params &params,default_random_engine& generator, SimulationFiles & sf);

void initSeed(Population & city,Params &params, default_random_engine & generator, SimulationFiles & sf);

void initCity(Population & city,Params & params);

void initNetworks(Networks & contacts, int N_use, Places & places, Params & params, default_random_engine& generator);

void initPlaces(Places & places, Population & city, Params & params);

void initDynamicalNetworks(Networks& contacts,int N_use, Params & params);

void convertParams(Params & params, MapParams & mapParams) ;

#endif
