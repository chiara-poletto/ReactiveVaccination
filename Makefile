# Makefile

CXX=g++
CXX_DEBUG_FLAGS=-g --std=c++11 -Wall -Wextra
CXXFLAGS=--std=c++11 -Wall -Wextra 
CXXFLAGSbasic= -O3 -Wall -Wextra  -lm 

OBJECT= routine_miscellanea.o inputOutput.o init.o transmission.o vaccination.o testing.o screening.o main.o 

all: exe

exe : clean-exe $(OBJECT)
	$(CXX) $(CXXFLAGS) -o vaccination.exe $(OBJECT) 

.PHONY: debug
debug: CXXFLAGS=$(CXX_DEBUG_FLAGS)
debug: clean-gdb $(OBJECT)
	$(CXX) $(CXXFLAGS) -o vaccination.gdb $(OBJECT) 

clean:
	rm -f *.o *.exe *.gdb
	
clean-gdb:
	rm -f *.o *.gdb

clean-exe:
	rm -f *.o *.exe

