#
# General Makefile for root-gravlens analysis examples
#
USER_SRCS = $(wildcard *.cc)
#
# Executable names come from the .cc sources in this directory.
# Replace the wildcard expression with .cc file list if you do
# not want to compile all .cc files in this directory
#
EXE = $(patsubst %.cc,%, $(wildcard *.cc))

#
#############################################################

CLASS_HEADERS = ProcessedRFPower.h

## You should not need to change anything below this line ###

.PHONY: all depend clean
## Getting the root flags
CXXFLAGS    = $(shell root-config --cflags) -std=c++11
LDFLAGS     = $(shell root-config --ldflags)
LDFLAGS    += $(shell root-config --libs)

## ANITA II flags
CXXFLAGS	 +=-I$(ANITA_UTIL_INSTALL_DIR)/include
LDFLAGS	 +=-L$(ANITA_UTIL_INSTALL_DIR)/lib -lAnitaEvent -lRootFftwWrapper -lMathMore -lfftw3 -lBrotterTools

all: $(EXE)

%: %.cc
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) -lMinuit 

#############################################################
# gcc can generate the dependency list

depend: Make-depend

Make-depend: $(USER_SRCS)
	$(CPP) $(CXXFLAGS) -MM $^ > $@
clean:
	- rm -f *.o  *.so *.ps core Make-depend $(EXE)

-include Make-depend
