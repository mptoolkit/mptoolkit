// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-opinfo.cpp
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
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/proccontrol.h"

int main(int argc, char** argv)
{
   if (argc != 2)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-opinfo <operator>\n";
      return 1;
   }

   int const Verbosity = 0;

   std::string TempFile = getenv_or_default("MP_BINPATH", std::string());
   int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   if (Verbosity >= 2)
   {
      std::cerr << "Using page size " << PageSize 
                << ", cache size " << CacheSize << '\n';
   }
   int TempFileDesc = ProcControl::CreateUniqueFile(TempFile);
   CHECK(TempFileDesc != -1);  // if this happens, then we failed to create a temporary file
   pheap::Initialize(TempFile, 1, PageSize, CacheSize);

   MPOperator Op = ParseOperator(argv[1]);
   std::cout.precision(13);
   std::cout << Op << '\n';

   pheap::Shutdown();
}
