// -*- C++ -*- $Id: mp-rename-symmetry.cpp 802 2007-12-08 06:55:30Z ianmcc $

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

LinearWavefunction RenameSymmetry(LinearWavefunction const& Psi, SymmetryList const& NewSL)
{
   LinearWavefunction Result(NewSL);
   LinearWavefunction::const_iterator I = Psi.end();
   while (I != Psi.begin())
   {
      --I;
      Result.push_front(RenameSymmetry(*I, NewSL));
   }
   return Result;
}

InfiniteWavefunction RenameSymmetry(InfiniteWavefunction const& Psi, SymmetryList const& NewSL)
{
   InfiniteWavefunction Result(Psi);
   Result.Psi = RenameSymmetry(Result.Psi, NewSL);
   Result.C_old = RenameSymmetry(Result.C_old, NewSL);
   Result.C_right = RenameSymmetry(Result.C_right, NewSL);
   Result.QShift = RenameSymmetry(Result.QShift, NewSL);
   return Result;
}

int main(int argc, char** argv)
{
   try
   {
      if (argc != 4)
      {
	 print_copyright(std::cerr);
	 std::cerr << "usage: mp-irename-symmetry <old> <new> <psi>\n";
	 return 1;
      }

      std::string OldQ = argv[1];
      std::string NewQ = argv[2];
      std::string FName = argv[3];
      pvalue_ptr<InfiniteWavefunction> Psi = pheap::OpenPersistent(FName, mp_pheap::CacheSize());

      int const s = Psi->GetSymmetryList().WhichSymmetry(OldQ);
      if (s == -1)
      {
	 std::cerr << "mp-irename-symmetry: warning: \"" + OldQ 
	    + "\" does not name a quantum number.\n";
	 pheap::ShutdownPersistent(Psi);
	 return 0; // this isn't really an error
      }

      // construct the new symmetry list
      SymmetryList NewSL = ChangeSymmetryName(Psi->Psi.GetSymmetryList(), OldQ, NewQ);

      Psi = pvalue_ptr<InfiniteWavefunction>(new InfiniteWavefunction(RenameSymmetry(*Psi, NewSL)));
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
