// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-load-0.7.3.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "matrixproduct/mpwavefunction-compat.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"

int main(int argc, char** argv)
{
   if (argc != 3)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-load-0.7.3 <input-0.7.3-wavefunction> <output-0.7.4-wavefunction>\n";
      return 1;
   }

   std::string InPsi = argv[1];
   std::string OutPsi = argv[2];

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
   pheap::Initialize(OutPsi, 1, PageSize, CacheSize);

   pvalue_ptr<WavefunctionCompat> Psi = pheap::ImportHeap(InPsi);
   pvalue_ptr<LinearWavefunction> PsiNew = new LinearWavefunction(Psi->Psi);
   Psi = pvalue_ptr<WavefunctionCompat>();
   pheap::ShutdownPersistent(PsiNew);
}
