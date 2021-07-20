/**
    Simulation of an ABM for COVID-19
    @file routine_miscellanea.cpp
    @author Pierre-Yves Boelle, Chiara Poletto
    @acknowledgment Jesus Moreno for first version
    @version 1.0 2020-04-16
    @license GPL-3.0-or-later
*/

#ifndef functions_hpp
#define functions_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <random>
#include <array>
#include <string>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <sstream>
#include <cmath>
#include <set>
#include <map>
#include <algorithm>

using namespace std;

// ----------------------------------------------------------------------------------------------------------------

int counter(vector<int> (&status_vec),int status1);

void print_all_int(vector<int> (&vec));

void print_all_double(vector<double> (&vec));

void save_csv_int_vec(string name1, vector<int> (&x1) );

void save_csv_double_vec(string name1, vector<double> (&x1) );

void save_csv_double_vec2(string name1, string name2, vector<double> (&x1), vector<double> (&x2) );

void save_csv_double_vec3(string name, string name1, string name2, string name3,vector<double> (&x1),vector<double> (&x2),vector<double> (&x3)  );

void save_csv_double_multi(string name1, vector<double> (&x1),vector<double> (&x2),vector<double> (&x3),vector<double> (&x4),vector<double> (&x5),vector<double> (&x6),vector<double> (&x7),int N1 );

void save_csv_double_all(string name1, vector<double> (&x1),vector<double> (&y1),vector<double> (&z1),vector<double> (&x2),vector<double> (&y2),vector<double> (&x3),vector<double> (&y3),vector<double> (&x4),vector<double> (&y4),vector<double> (&x5),vector<double> (&y5),vector<double> (&x6),vector<double> (&y6),vector<double> (&x7),vector<double> (&y7),int N1 );

void save_csv_int_multi(string name1, vector<int> (&x1),vector<int> (&x2),vector<int> (&x3),vector<int> (&x4),vector<int> (&x5),vector<int> (&x6),vector<int> (&x7),int N1 );

int above_counter(vector<int> (&vector1),int index1);

vector<double> random_vector(default_random_engine& generator, int size);

template<typename K, typename V>
void print_map(map<K,V> const &m)
{
    for (auto const& pair: m) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
}

template<class bidiiter>
bidiiter random_unique(bidiiter begin, bidiiter end, size_t num_random) {
    size_t left = std::distance(begin, end);
    while (num_random--) {
        bidiiter r = begin;
        std::advance(r, rand()%left);
        std::swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}

#endif /* functions_hpp */
