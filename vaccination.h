/**
    Simulation of an ABM for COVID-19
    @file vaccination.h
    @author Pierre-Yves Boelle, Chiara Poletto
    @version 1.0 2020-04-16
    @license GPL-3.0-or-later
*/

#ifndef VACCINATION_H
#define VACCINATION_H


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

void initVaccination(Vaccination & vaccination, Population & city, Places & places,
		Params &params,default_random_engine& generator, SimulationFiles & sf);

void vaccinationStep(int real, int it, Vaccination & vaccination,  Population & city, Places &places, Compartments & epid,
		Params & params, default_random_engine & generator, SimulationFiles & simulationFiles);

void vaccinate_elderly(default_random_engine& generator, Population & city,
		int frac_elderly_vacc);
#endif
