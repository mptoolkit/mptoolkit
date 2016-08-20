// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-trunc.cpp
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
#include "matrixproduct/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

// a version of truncate, where we do not apply the truncation operator immediately,
// but instead save it for later and apply all of the truncations together at the end.
// This has the property that the truncations commute.
void truncate_global(LinearWavefunction& Psi, StatesInfo const& SInfo, bool ShowStates)
{
   LinearWavefunction::iterator I = Psi.begin();
   MatrixOperator M = MatrixOperator::make_identity(I->Basis1());
   int BondNr = 1;
   std::vector<MatrixOperator> Truncators;
   while (I != Psi.end())
   {
      *I = prod(M, *I);
      M = ExpandBasis2(*I);
      DensityMatrix<MatrixOperator> DM(scalar_prod(M, herm(M)));
      TruncationInfo Info;
      MatrixOperator U = DM.ConstructTruncator(DM.begin(), TruncateFixTruncationErrorAbsolute(DM.begin(),
                                                                                              DM.end(),
                                                                                              SInfo,
                                                                                              Info));
      if (ShowStates)
         std::cerr << "bond=" << BondNr
                   << ", states=" << Info.KeptStates()
                   << ", trunc=" << Info.TruncationError()
                   << ", largest_discarded_evalue=" << Info.LargestDiscardedEigenvalue()
                   << '\n';

      // Now we want to save the truncator that keeps all the states
      MatrixOperator UFull = DM.ConstructTruncator(DM.begin(),
                                                   DM.begin() + M.Basis2().total_dimension());

      *I = prod(*I, herm(UFull));
      M = UFull*M;

      U = U*herm(UFull);
      Truncators.push_back(U);

      ++I;
      ++BondNr;
   }

   // Now we apply the truncators
   if (ShowStates)
      std::cerr << "Applying truncators...\n";

   I = Psi.begin();
   std::vector<MatrixOperator>::const_iterator TI = Truncators.begin();
   *I = prod(*I, herm(*TI));
   ++I;
   while (I != Psi.end())
   {
      *I = prod(*TI, *I);
      ++TI;
      *I = prod(*I, herm(*TI));
      ++I;
   }
   M = (*TI)*M;

   DEBUG_CHECK(++TI == Truncators.end());

   if (ShowStates)
      std::cerr << "Orthogonalizing...\n";

   Psi = inject_right_old_interface(Psi, M);
   I = Psi.begin();
   *I = prod(M, *I);

   if (ShowStates)
      std::cerr << "Done.\n";
}

int main(int argc, char** argv)
{
   try
   {
      int MinStates = 1;
      int MaxStates = 100000;
      double TruncCutoff = 0;
      double EigenCutoff = -1;
      std::string WavefunctionFile;
      bool Verbose = false;
      bool Normalize = false;
      bool Global = false;
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("wavefunction,w", prog_opt::value(&WavefunctionFile),
          "input wavefunction (required)")
         ("out,o", prog_opt::value<std::string>(),
          "output wavefunction (overwrites the input if this is not specified)")
         ("min-states", prog_opt::value<int>(&MinStates),
          ("Minimum number of states to keep [default "
           +boost::lexical_cast<std::string>(MinStates)+"]").c_str())
         ("max-states,m", prog_opt::value<int>(&MaxStates),
          ("Maximum number of states to keep [default "
           +boost::lexical_cast<std::string>(MaxStates)+"]").c_str())
         ("trunc,r", prog_opt::value(&TruncCutoff),
          ("Cutoff truncation error per site [default "
           +boost::lexical_cast<std::string>(TruncCutoff)+"]").c_str())
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
          ("Cutoff threshold for density matrix eigenvalues [default "
           +boost::lexical_cast<std::string>(EigenCutoff)+"]").c_str())
         ("global", prog_opt::bool_switch(&Global),
          "Truncate globally - determine all truncation operators before applying any of them")
         ("normalize", prog_opt::bool_switch(&Normalize),
          "Normalize the output wavefunction to be the same as the input")
         ("verbose,v", prog_opt::bool_switch(&Verbose),
          "show additional information")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-trunc [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<MPWavefunction> Psi;
      if (vm.count("out"))
      {
         int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
         pheap::Initialize(vm["out"].as<std::string>(), 1, PageSize, CacheSize);
         Psi =  pheap::ImportHeap(WavefunctionFile);
      }
      else
      {
         Psi = pheap::OpenPersistent(WavefunctionFile, CacheSize);
      }

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = MaxStates;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;

      if (Verbose)
         std::cerr << SInfo << '\n';

      MPWavefunction p = *Psi;
      double Norm = norm_2(p);

      if (Global)
         truncate_global(p, SInfo, Verbose);
      else
         truncate(p, SInfo, Verbose);

      if (Normalize)
         p *= Norm / norm_2(p);
      Psi = new MPWavefunction(p);

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
