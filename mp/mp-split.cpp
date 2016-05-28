// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-split.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc != 4)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-split <psi> <real-part> <imag-part>\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(argv[1], CacheSize, true);

   MPWavefunction Psi = *PsiPtr;
   MPWavefunction PsiBar = conj(Psi);

   {
      pvalue_ptr<MPWavefunction> PsiReal = new MPWavefunction(0.5 * (Psi + PsiBar));
      //   }
      //   {
      pvalue_ptr<MPWavefunction> PsiImag = 
         new MPWavefunction(complex(0.0,-0.5) * (Psi - PsiBar));
      pheap::ExportHeap(argv[3], PsiImag);
      pheap::ExportHeap(argv[2], PsiReal);
   }
}
