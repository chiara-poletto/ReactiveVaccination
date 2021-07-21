# ReactiveVaccination
README

This repository contains file for the individual based model used in the companion paper listed below :
Reactive vaccination in Schools and Workplaces for COVID-19
by Faucher, B et al.

############################
#INSTALLATION
############################

The code must be compiled with g++. Tested with g++ 5.2.0
The code uses standard c++ libraries

source code (.cpp,.h,Makefile) in src folder
filenames_test.txt in src folder
tests.tgz in src/input/test folder

make all # takes <60 secs

executable file is vaccination.exe

############################
#DEMO
############################

# get help
./vaccination.exe -h 

-h : print this help

-n N runs scenario N

-s SEED sets random seed to SEED (0 for time)

-f FILENAME for input (default filenames.txt)

-i 0 for computing (natural) initial immunity/ -i 1 constrained to max incidence

-j 0 (random seed) / 1 (seed only exposed) / 2 (seed exposed and others) / 3 (seed all exposed file) / 4 (same as 2 but random)

# parameter file - list of all parameters
input/params_test.csv

#run demo
mkdir src/output
./vaccination.exe -n2 -f filenames_test.txt

#output
files in output/CT_vrandomtype4_vdaily50_b013_im20_c74482cf_test
run time <30 s

##############################
#
##############################
 
