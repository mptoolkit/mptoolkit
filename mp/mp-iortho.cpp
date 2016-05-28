// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-iortho.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "matrixproduct/infinitewavefunction.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"

int main(int argc, char** argv)
{
   if (argc != 2)
   {
      std::cout << "usage: mp-iortho <wavefunction>\n";
      return 1;
   }

   std::string FName = argv[1];

   pvalue_ptr<InfiniteWavefunction> Psi = pheap::OpenPersistent(FName, mp_pheap::CacheSize(), true);

   std::cout.precision(14);
   std::cout << (1.0 - orthogonality_fidelity(*Psi)) << '\n';

   pheap::Shutdown();
}
