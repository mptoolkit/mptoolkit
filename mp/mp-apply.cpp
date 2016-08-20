// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-apply.cpp
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"
#include "interface/operator-parser.h"

int main(int argc, char** argv)
{
   if (argc < 4 || argc > 6)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-apply <operator> <input-wavefunction> <output-wavefunction> [states] [<quantumnumber>]\n";
      return 1;
   }

   std::string InPsi = argv[2];
   std::string OutPsi = argv[3];

   pvalue_ptr<MPWavefunction> Psi;

   long CacheSize = getenv_or_default("MP_CACHESIZE", DEFAULT_PAGE_CACHE_SIZE);
   if (InPsi == OutPsi)
      Psi = pheap::OpenPersistent(InPsi.c_str(), CacheSize);
   else
   {
      int PageSize = getenv_or_default("MP_PAGESIZE", DEFAULT_PAGE_SIZE);
      pheap::Initialize(OutPsi, 1, PageSize, CacheSize);
      Psi = pheap::ImportHeap(InPsi.c_str());
   }

   MPOperator Op = ParseOperator(argv[1]);

   if (argc >= 5)
   {
      WARNING("mp-apply does not currently obey the [states] option.");
   }

   int NStates = DefaultMaxStates;
   if (argc >= 5)
      NStates = boost::lexical_cast<int>(argv[4]);

   QuantumNumbers::QuantumNumber Q;
   if (argc >= 6)
   {
      Q = QuantumNumbers::QuantumNumber(Op.GetSymmetryList(), argv[5]);
   }
   else
   {
      QuantumNumbers::QuantumNumberList QList = transform_targets(Op.TransformsAs(), Psi->TransformsAs());
      if (QList.size() != 1)
      {
         std::cerr << "quantum number of operator-wavefunction product is ambiguous,\n"
                   << "  choices are ";
         std::copy(QList.begin(), QList.end(), std::ostream_iterator<QuantumNumber>(std::cerr, " "));
         std::cerr << '\n';
         exit(1);
      }
      Q = QList[0];
   }

   *Psi.mutate() = prod(Op, *Psi, Q); // , NStates);

   pheap::ShutdownPersistent(Psi);
}
