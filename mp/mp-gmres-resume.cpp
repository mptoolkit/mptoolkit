// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "common/conflist.h"
#include "mp-algorithms/solver-gmres.h"
#include "mp-algorithms/stateslist.h"
#include "mp-algorithms/dmrgloop.h"
#include "mp-algorithms/resume.h"
#include <iostream>

int main(int argc, char** argv)
{
   return GenericResume<DMRGLoop<SolverGmres> >(argc, argv, "mp-gmres-resume");
}
