
ROOT=../../

# directories to look for include files
INCDIRS=$(ROOT)

# directories to look for source files
SRCDIRS=$(ROOT)/common

# libraries, eg LIBS=lapack blas
LIBS=

# compiler and compiler flags.  By default we use lots of debugging
CXX=g++
CXXFLAGS=-std=c++14 -Wall -Wextra -Wconversion -Wno-sign-conversion -Wno-conversion \
          -Wno-shorten-64-to-32 -g -D_GLIBCXX_DEBUG

CXXFLAGS=-std=c++14 -O2

# compiler flags for optimization.  Uncomment this line to enable.
# Don't forget to 'make clean' before recompiling.
# CXXFLAGS=-std=c++11 -O3 -march=native -DNDEBUG

# list of executables (so that 'make clean' can remove them)
executable = interpolate

# dependencies for our executable programs.  The first target listed
# will be the default if 'make' is run with no arguments

# example : example.o randutil.o

# link executable 'example' with lapack and blas libraries,
# alternative to global LIBS above
# example_LIBS = lapack blas

all : $(executable)

#
# everything bellow here is boilerplate incantations,
# do not modify!
#

# turn off some implicit rules that we don't want
% : %.cpp
% : %.cc
%.o : %.cpp
%.o : %.cc

# get all cpp files in VPATH
cxxsrcglob = *.cpp $(patsubst %,%/*.cpp, $(subst :, ,$(VPATH)))
cxx_source := $(notdir $(wildcard $(cxxsrcglob)))
cxx_depend = $(cxx_source:%.cpp=%.d)

# these targets don't require (or want) dependency information
nodep-targets = clean

# don't delete intermediate .d files
.SECONDARY: $(cxx_depend)

# phony targets; run rule unconditionally
.PHONY: clean all

VPATH=$(SRCDIRS)

# convert the list of include directories x y z ... into
# -Ix -Iy -Iz ....
INCLUDE=$(INCDIRS:%=-I%)

# add the include directories to the compiler flags
CXXFLAGS+=$(INCLUDE)

%.o : %.cpp %.d
	$(CXX) $(CXXFLAGS) -MM $< -o $*.d
	@cp -f $*.d $*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp
	$(CXX) $(CXXFLAGS) -c $< -o $@

% : %.o
	$(CXX) $(CXXFLAGS) $^  $(LIB:%=-l%) $($@_LIBS:%=-l%) -o $@

%.d :
	\

clean :
	rm -f *.o *.d $(executable)

# include the dependency information
ifeq (,$(filter $(nodep-targets),$(MAKECMDGOALS)))
-include $(cxx_depend)
endif
