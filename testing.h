/**
    Simulation of an ABM for COVID-19
    @file testing.h
    @author Pierre-Yves Boelle, Chiara Poletto
    @version 1.0 2020-04-16
    @license GPL-3.0-or-later
*/

#ifndef TRACING_H
#define TRACING_H



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

void initTesting(Testing & testing, Population & city);

int testThisPerson(int j, int loc, Population & city, Testing & testing, Compartments & upd);

void isolateContactsInHousehold (int i, Testing & testing, Population & city, Networks & contacts, Compartments & upd, Params & params,
		default_random_engine & generator, SimulationFiles & simulationFiles,
		int real, int it);

void initRegularScreening(Testing & testing, Places & places, Params & params,default_random_engine generator,SimulationFiles & sf);

void testingStep(int real, int it, Testing & testing, Population & city, Networks &contacts,
		Compartments &upd, Params & params, default_random_engine & generator,SimulationFiles & simulationFiles);

#endif
