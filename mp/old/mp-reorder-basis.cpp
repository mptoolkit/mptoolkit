// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-reorder-basis.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>
#include "interface/inittemp.h"
#include "interface/operator-parser.h"
#include "common/terminal.h"
#include "common/environment.h"
#include <boost/program_options.hpp>
#include "common/stringutil.h"  // for Split

namespace prog_opt = boost::program_options;

LinearWavefunction ReorderLocalBasis(LinearWavefunction const& Psi,
                                     std::list<int> const& NewOrder)
{
   LinearWavefunction Result(Psi.GetSymmetryList());
   LinearWavefunction::const_iterator I = Psi.end();
   while (I != Psi.begin())
   {
      --I;
      Result.push_front(ReorderLocalBasis(*I, NewOrder));
   }
   return Result;
}

template <typename T>
struct DoLexicalCast
{
   typedef T result_type;
   typedef std::string argument_type;

   result_type operator()(argument_type const& x) const
   {
      return boost::lexical_cast<result_type>(x);
   }
};

int main(int argc, char** argv)
{
   try
   {
      if (argc != 3)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-reorder-basis <new-ordering> <psi>\n"
            "the new-ordering is a comma-delimited permutation of numbers 0..(N-1).\n";
         return 1;
      }

      std::string Ordering = argv[1];
      std::string FName = argv[2];
      std::list<std::string> Strings;
      Split(Ordering, ',', std::back_inserter(Strings));
      std::list<int> NewOrder;
      std::transform(Strings.begin(), Strings.end(),
                     std::back_inserter(NewOrder), DoLexicalCast<int>());
      pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(FName, mp_pheap::CacheSize());

      Psi = pvalue_ptr<MPWavefunction>(new MPWavefunction(ReorderLocalBasis(*Psi, NewOrder)));
      pheap::ShutdownPersistent(Psi);
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      return 1;
   }
}
