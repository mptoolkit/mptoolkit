// -*- C++ -*- $Id$

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

namespace prog_opt = boost::program_options;

LinearWavefunction ReorderSymmetry(LinearWavefunction const& Psi, SymmetryList const& NewSL)
{
   LinearWavefunction Result(NewSL);
   LinearWavefunction::const_iterator I = Psi.end();
   while (I != Psi.begin())
   {
      --I;
      Result.push_front(CoerceSymmetryList(*I, NewSL));
   }
   return Result;
}

int main(int argc, char** argv)
{
   try
   {
      if (argc != 3)
      {
	 print_copyright(std::cerr);
	 std::cerr << "usage: mp-reorder-symmetry <new-symmetry-list> <psi>\n";
	 return 1;
      }

      std::string NewSLName = argv[1];
      std::string FName = argv[2];
      pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(FName, mp_pheap::CacheSize());

      SymmetryList NewSL(NewSLName);

      Psi = pvalue_ptr<MPWavefunction>(new MPWavefunction(ReorderSymmetry(*Psi, NewSL)));
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
