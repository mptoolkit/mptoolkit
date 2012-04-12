// -*- C++ -*- $Id: mp-reorder-symmetry.cpp 802 2007-12-08 06:55:30Z ianmcc $

#include "mps/infinitewavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>
#include "interface/inittemp.h"
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

InfiniteWavefunction ReorderSymmetry(InfiniteWavefunction const& Psi, SymmetryList const& NewSL)
{
   InfiniteWavefunction Result;
   Result.Psi = ReorderSymmetry(Psi.Psi, NewSL);
   Result.C_old = Psi.C_old;
   CoerceSymmetryList(Result.C_old, NewSL);
   Result.C_right = Psi.C_right;
   CoerceSymmetryList(Result.C_right, NewSL);
   Result.QShift = Psi.QShift;
   CoerceSymmetryList(Result.QShift, NewSL);
   Result.Attr = Psi.Attr;
   return Result;
}

int main(int argc, char** argv)
{
   try
   {
      if (argc != 3)
      {
	 print_copyright(std::cerr);
	 std::cerr << "usage: mp-ireorder-symmetry <new-symmetry-list> <psi>\n";
	 return 1;
      }

      std::string NewSLName = argv[1];
      std::string FName = argv[2];
      pvalue_ptr<InfiniteWavefunction> Psi = pheap::OpenPersistent(FName, mp_pheap::CacheSize());

      SymmetryList NewSL(NewSLName);

      Psi = pvalue_ptr<InfiniteWavefunction>(new InfiniteWavefunction(ReorderSymmetry(*Psi, NewSL)));
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
