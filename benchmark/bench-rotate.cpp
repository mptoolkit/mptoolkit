// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// benchmark/bench-rotate.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"

void DoTest(MPWavefunction& x)
{
   while (x.RightSize() > 1)
      x.RotateRight();

   while (x.LeftSize() > 1)
      x.RotateLeft();
}

int main(int argc, char** argv)
{
   if (argc != 2)
   {
      std::cerr << "usage: benchrotate <wavefunction>\n";
      return 1;
   }

   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(argv[1], 655360, true);

   DoTest(*Psi.mutate());
}
