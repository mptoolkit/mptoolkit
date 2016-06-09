// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-dmrg-resume.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

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
