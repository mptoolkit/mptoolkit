// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-change-lattice.cpp
//
// Copyright (C) 2015-2023 Ian McCulloch <ian@qusim.net>
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

#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcell-parser.h"
#include "common/statistics.h"

namespace prog_opt = boost::program_options;


int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string Lat1Str;
      std::string Lat2Str;
      std::string InputFile;
      std::string OutputFile;
      bool Force = false;
      bool AssumeOrthogonal = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("force,f", prog_opt::bool_switch(&Force),
          "allow overwriting the output file, if it already exists")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("lat1", prog_opt::value(&Lat1Str), "lat1")
         ("lat2", prog_opt::value(&Lat2Str), "lat2")
         ("psi1", prog_opt::value(&InputFile), "psi1")
         ("psi2", prog_opt::value(&OutputFile), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("lat1", 1);
      p.add("lat2", 1);
      p.add("psi1", 1);
      p.add("psi2", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("psi2") < 1)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <current_lattice> <new_lattice> <input-psi> <output-psi>\n";
         std::cerr << desc << '\n';
         std::cerr << "This tool changes the local Hilbert space of an MPS to a new lattice.\n"
                   << "It does this by matching the state names in the local Hilbert space to the new lattice.\n"
                   << "A typical use for this is for Bosonic models, to adjust the maximum number of bosons per site,\n"
                   << "but it can also be used to effect a Gutzwiller projection, or similar projections.\n"
                   << "If there are some states in the current Hilbert space that are not present in the new lattice, then\n"
                   << "those states are projected out, leaving the wavefunction with reduced norm.\n"
                   << "If there are states in the new lattice that are not in the current Hilbert space, they are set to zero.\n";
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (Verbose > 0)
         std::cout << "Loading wavefunction..." << std::endl;

      pvalue_ptr<MPWavefunction> PsiPtr;
      if (InputFile == OutputFile)
         PsiPtr = pheap::OpenPersistent(InputFile.c_str(), mp_pheap::CacheSize());
      else
      {
         pheap::Initialize(OutputFile, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);
         PsiPtr = pheap::ImportHeap(InputFile);
      }

      FiniteWavefunctionLeft Psi = PsiPtr->get<FiniteWavefunctionLeft>();

      pvalue_ptr<InfiniteLattice> Lattice1Ptr = pheap::ImportHeap(Lat1Str);
      InfiniteLattice Lattice1 = *Lattice1Ptr;

      pvalue_ptr<InfiniteLattice> Lattice2Ptr = pheap::ImportHeap(Lat2Str);
      InfiniteLattice Lattice2 = *Lattice2Ptr;

      if (Verbose > 0)
         std::cout << "Mapping wavefunction" << std::flush;
      if (Verbose > 1)
         std::cout << std::endl;
      std::list<StateComponent> NewPsi;
      UnitCell::const_iterator CurrentSite = Lattice1.GetUnitCell().begin();
      UnitCell::const_iterator NewSite = Lattice2.GetUnitCell().begin();
      int s = 0;
      for (auto I = Psi.begin(); I != Psi.end(); ++I)
      {
         if (Verbose == 1)
         {
            std::cout << '.';
         }
         if (Verbose > 1)
         {
            std::cout << "Mapping wavefunction site " << s << std::endl;
         }
         if (Verbose > 2)
         {
            std::cout << "Current Hilbert space dimension: " << (CurrentSite->Basis().size()) << std::endl;
            std::cout << "New Hilbert space dimension: " << (NewSite->Basis().size()) << std::endl;
         }
         int c = 0; // count of the number of states that we've mapped into the new basis
         StateComponent A(NewSite->Basis(), I->Basis1(), I->Basis2());
         for (int i = 0; i < NewSite->Basis().size(); ++i)
         {
            std::string Label = NewSite->Basis().Label(i);
            int j = CurrentSite->Basis().LookupOrNeg(Label);
            if (j >= 0)
            {
               if (Verbose > 3)
               {
                  std::cout << "Mapping state " << Label << " into the new basis." << std::endl;
               }
               ++c;
               A[i] = (*I)[j];
            }
            else if (Verbose > 3)
            {
               std::cout << "New basis state " << Label << " does not exist in the current basis; will be set to zero." << std::endl;
            }
         }
         if (c == 0)
         {
            std::cout << "mp-change-lattice: fatal: new lattice at site " << s << " has no states in common in the local Hilbert spaces!\n";
            exit(1);
         }
         if (Verbose > 2)
         {
            std::cout << "Mapped " << c << " basis states." << std::endl;
         }
         ++s;
         if (++CurrentSite == Lattice1.GetUnitCell().end())
            CurrentSite = Lattice1.GetUnitCell().begin();
         if (++NewSite == Lattice2.GetUnitCell().end())
         {
            NewSite = Lattice2.GetUnitCell().begin();
         }
         NewPsi.push_back(std::move(A));
      }
      if (Verbose == 1)
         std::cout << std::endl;
      if (Verbose > 0)
      {
         std::cout << "Finished mapping wavefunction." << std::endl;
         std::cout << "Orthogonalizing wavefunction..." << std::endl;
      }

      PsiPtr.mutate()->Wavefunction() = FiniteWavefunctionLeft::Construct(LinearWavefunction::FromContainer(NewPsi.begin(), NewPsi.end()));

      PsiPtr.mutate()->AppendHistoryCommand(EscapeCommandline(argc, argv));
      PsiPtr.mutate()->SetDefaultAttributes();

      if (Verbose > 0)
         std::cout << "Finished." << std::endl;

      pheap::ShutdownPersistent(PsiPtr);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      pheap::Cleanup();
      return 1;
   }
}
