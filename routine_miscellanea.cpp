//
/**
    Simulation of an ABM for COVID-19
    @file routine_miscellanea.cpp
    @author Pierre-Yves Boelle, Chiara Poletto
    @acknowledgment Jesus Moreno for first version
    @version 1.0 2020-04-16
    @license GPL-3.0-or-later
*/

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
#include <algorithm>

#include "defs.h"
#include "routine_miscellanea.hpp"
using namespace std;

// ----------------------------------------------------------------------------------------------------------------
string getFileName(std::string filePath, char separator='/')
{
    // Get last dot position
    size_t sepPos = filePath.rfind(separator);
    if(sepPos != string::npos)
    {
        return filePath.substr(sepPos + 1, filePath.size());
    }
    return "";
}

// Function to produce a vector of random numbers
vector<double> random_vector(default_random_engine& generator, int size)
{
    uniform_real_distribution<double> distribution(0,1);
    vector<double> vec;
    vec.resize(size);

    for (int i=0; i< (int) vec.size();++i)
    {
        vec[i] = distribution(generator);
    }

    return vec;
}



// Function to count people with status == status1 -------------------------

int counter(vector<int> & status_vec, int status1)
{
	return count(status_vec.cbegin(), status_vec.cend(),status1);
}


// Function to print all elements of vector of type int --------------------

void print_all_int(vector<int> (&vec))
{for (int v: vec )
    {cout << v << " ";}
    cout.flush();
    cout << '\n';}

// Function to print all elements of vector of type double -----------------

void print_all_double(vector<double> (&vec))
{for (double v: vec )
    {cout << v << " ";}
    cout.flush();
    
    cout << '\n';}


// Function to save a single vector to CSV when vector is int --------------------------------------

void save_csv_int_vec(string name1, vector<int> (&x1) )
{
    
    string name =name1 + ".csv";
    string name2 = getFileName(name1);


ofstream csvfile(name);
if (csvfile.is_open())
{
    // it is possible that the vectors are not of the same size
    // if you get the maximum, then you may need to fill in empty fields
    //     for some of the columns
    
    csvfile << "Nonsense" << "," << name2 << endl;
    
    int num_of_rows = x1.size();
    
    for (int i = 0; i < num_of_rows; i++)
    {
        csvfile << i << "," << x1[i] << endl;
    }
    csvfile.close();
}
else
{
    cout << "File could not be opened";
}
}

// Function to save a single vector to CSV when vector is double --------------------------------------

void save_csv_double_vec(string name1, vector<double> (&x1) )
{
    
    string name =name1 + ".csv";
    string name2 = getFileName(name1);
ofstream csvfile(name);
if (csvfile.is_open())
{
    // it is possible that the vectors are not of the same size
    // if you get the maximum, then you may need to fill in empty fields
    //     for some of the columns
    
    csvfile << "Nonsense" << "," << name2 << endl;
    
    int num_of_rows = x1.size();
    
    for (int i = 0; i < num_of_rows; i++)
    {
        csvfile << i << "," << x1[i] << endl;
    }
    csvfile.close();
}
else
{
    cout << "File could not be opened";
}
}

// Function to save 2 double vectors to CSV when vector is double --------------------------------------

void save_csv_double_vec2(string name1, string name2,vector<double> (&x1),vector<double> (&x2) )
{

    string name =name1 + ".csv";

ofstream csvfile(name);
if (csvfile.is_open())
{
    // it is possible that the vectors are not of the same size
    // if you get the maximum, then you may need to fill in empty fields
    //     for some of the columns

    csvfile << "Nonsense" << "," << name2 << endl;

    int num_of_rows = x1.size();

    for (int i = 0; i < num_of_rows; i++)
    {
        csvfile << i << "," << x1[i] << "," << x2[i]<<endl;
    }
    csvfile.close();
}
else
{
    cout << "File could not be opened";
}
}

// Function to save 2 double vectors to CSV when vector is double --------------------------------------

void save_csv_double_vec3(string name, string name1, string name2, string name3,vector<double> (&x1),vector<double> (&x2),vector<double> (&x3)  )
{

    name =name + ".csv";

ofstream csvfile(name);
if (csvfile.is_open())
{
    // it is possible that the vectors are not of the same size
    // if you get the maximum, then you may need to fill in empty fields
    //     for some of the columns

    csvfile << "Nonsense" << "," << name1 <<","<< name2<<","<<name3 << endl;

    int num_of_rows = x1.size();

    for (int i = 0; i < num_of_rows; i++)
    {
        csvfile << i << "," << x1[i] << "," << x2[i]<< "," << x3[i] << endl;
    }
    csvfile.close();
}
else
{
    cout << "File could not be opened";
}
}

// Function to print out data to CSV when vectors are double --------------------------------------

void save_csv_double_multi(string name1, vector<double> (&x1),vector<double> (&x2),vector<double> (&x3),vector<double> (&x4),vector<double> (&x5),vector<double> (&x6),vector<double> (&x7),int N1 )
{
    
    string name =name1 + ".csv";
    
ofstream csvfile(name);
if (csvfile.is_open())
{
    // it is possible that the vectors are not of the same size
    // if you get the maximum, then you may need to fill in empty fields
    //     for some of the columns
    int num_of_rows = min({ x1.size(), x2.size(), x3.size(),x4.size(),x5.size(),x6.size(),x7.size() });
    
    csvfile << "Nonsense" << "," << "Time" << "," << "Susceptibles" << "," << "Exposed" << "," << "Presymptomatic1" << "," << "Symptomatics" << "," << "Presymptomatic2" << "," << "Asymptomatics" << "," << "Recovered" << "," << "Total" << endl;
    
    for (int i = 0; i < num_of_rows; i++)
    {
        csvfile << i << "," << i << "," << x1[i]/N1 << "," << x2[i]/N1 << "," << x3[i]/N1 << "," << x4[i]/N1 << "," << x5[i]/N1 << "," << x6[i]/N1 << "," << x7[i]/N1 << "," << x2[i]/N1 + x3[i]/N1 + x4[i]/N1 + x5[i]/N1 + x6[i]/N1 <<  endl;
    }
    csvfile.close();
}
else
{
    cout << "File could not be opened";
}
}

// Function to print out data to CSV when vectors are double --------------------------------------

void save_csv_double_all(string name1, vector<double> (&x1),vector<double> (&y1),vector<double> (&z1),vector<double> (&x2),vector<double> (&y2),vector<double> (&x3),vector<double> (&y3),vector<double> (&x4),vector<double> (&y4),vector<double> (&x5),vector<double> (&y5),vector<double> (&x6),vector<double> (&y6),vector<double> (&x7),vector<double> (&y7),int N1 )
{
    
    string name =name1 + ".csv";
    
ofstream csvfile(name);
if (csvfile.is_open())
{
    // it is possible that the vectors are not of the same size
    // if you get the maximum, then you may need to fill in empty fields
    //     for some of the columns
    int num_of_rows = min({ x1.size(), x2.size(), x3.size(),x4.size(),x5.size(),x6.size(),x7.size() });
    
    csvfile << "Nonsense" << "," << "Time" << "," << "S" << "," << "vS" << "," << "wvS" << "," << "E" << "," << "vE" << "," << "P1" << "," << "vP1" << "," << "Symp" << "," << "vSymp" << "," << "P2" << "," << "vP2" << "," << "A"  << "," << "vA" << "," << "R"  << "," << "vR" << "," << "Total" << endl;
    
    for (int i = 0; i < num_of_rows; i++)
    {
        csvfile << i << "," << i << "," << x1[i]/N1  << "," << y1[i]/N1 << "," << z1[i]/N1 << "," << x2[i]/N1 << "," << y2[i]/N1 << "," << x3[i]/N1  << "," << y3[i]/N1 << "," << x4[i]/N1  << "," << y4[i]/N1 << "," << x5[i]/N1  << "," << y5[i]/N1  << "," << x6[i]/N1 << "," << y6[i]/N1 << "," << x7[i]/N1 << "," << y7[i]/N1 << "," << x2[i]/N1 + x3[i]/N1 + x4[i]/N1 + x5[i]/N1 + x6[i]/N1 + y2[i]/N1 + y3[i]/N1 + y4[i]/N1 + y5[i]/N1 + y6[i]/N1 <<  endl;
    }
    csvfile.close();
}
else
{
    cout << "File could not be opened";
}
}

// Function to print out data to CSV when vectors are int --------------------------------------

void save_csv_int_multi(string name1, vector<int> (&x1),vector<int> (&x2),vector<int> (&x3),vector<int> (&x4),vector<int> (&x5),vector<int> (&x6),vector<int> (&x7),int N1 )
{
    
    string name =name1 + ".csv";
    
ofstream csvfile(name);
if (csvfile.is_open())
{
    // it is possible that the vectors are not of the same size
    // if you get the maximum, then you may need to fill in empty fields
    //     for some of the columns
    int num_of_rows = min({ x1.size(), x2.size(), x3.size(),x4.size(),x5.size(),x6.size(),x7.size() });
    
    csvfile << "Nonsense" << "," << "Time" << "," << "Susceptibles" << "," << "Exposed" << "," << "Presymptomatic1" << "," << "Symptomatics" << "," << "Presymptomatic2" << "," << "Asymptomatics" << "," << "Recovered" << "," << "Total" << endl;
    
    for (int i = 0; i < num_of_rows; i++)
    {
        csvfile << i << "," << i << "," << 1.0*x1[i]/N1 << "," << 1.0*x2[i]/N1 << "," << 1.0*x3[i]/N1 << "," << 1.0*x4[i]/N1 << "," << 1.0*x5[i]/N1 << "," << 1.0*x6[i]/N1 << "," << 1.0*x7[i]/N1 << "," << 1.0*x2[i]/N1 + 1.0*x3[i]/N1 + 1.0*x4[i]/N1 + 1.0*x5[i]/N1 + 1.0*x6[i]/N1 <<  endl;
    }
    csvfile.close();

}
else
{
    cout << "File could not be opened";
}
}

// Function to count integer vector elements above a certain index.

int above_counter(vector<int> (&vector1),int index1)
{
    int counts{0};
    
    for (int v:vector1)
    {
        if (v>=index1)
        {++counts;}
    }
    
    return counts;
}



