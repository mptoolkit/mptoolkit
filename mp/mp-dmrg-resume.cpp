// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "common/conflist.h"
#include "quantumnumbers/u1.h"
#include "quantumnumbers/su2.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "mp-algorithms/dmrg.h"
#include "mp-algorithms/stateslist.h"
#include "mp-algorithms/dmrgloop.h"
#include "mp-algorithms/resume.h"
#include "mp/copyright.h"
#include <iostream>

int main(int argc, char** argv)
{
   return GenericResume<DMRGLoop<DMRG> >(argc, argv, "mp-dmrg-resume");
}
