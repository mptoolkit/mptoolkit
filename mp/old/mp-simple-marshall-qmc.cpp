// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-simple-marshall-qmc.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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
#include "mp-algorithms/random_wavefunc.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "interface/inittemp.h"
#include "common/environment.h"
#include <common/prog_options.h>
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include <tuple>

namespace prog_opt = boost::program_options;

typedef std::complex<double> complex;

// attempt to flip basis states at two sites of the configuration such that
// it stays in the same symmetry sector.
bool TryFlipPair(WavefunctionDesc& Config, std::vector<BasisList> const& Basis, int s1, int s2)
{
   // firstly, find a new configuration for site s1
   int NewState = rand() % Basis[s1].size();
   Projection CurrentDelta = difference(Config.Height[s1], Config.Height[s1+1]);
   ProjectionList PL = enumerate_projections(Basis[s1][NewState]);
   Projection NewDelta = PL[rand() % PL.size()];

   // now our action depends on whether s2 > s1 or s2 < 1
   if (s2 == s1)
   {
      // the corner case, but perhaps it is possible
      if (NewDelta == CurrentDelta)
      {
         Config.State[s1] = NewState;
         return true;
      }
      return false;
   }
   // else

   // try and flip s2 such that the change is opposite that of s1
   Projection DeltaChange = sum(NewDelta, negate(CurrentDelta));
   Projection s2CurrentDelta = difference(Config.Height[s2], Config.Height[s2+1]);
   Projection s2WantDelta = sum(s2CurrentDelta, negate(DeltaChange)); // this is the delta we require for s2
   //   TRACE(DeltaChange)(s2CurrentDelta)(s2WantDelta);
   // try and find a basis state in s2 that permits s2WantDelta
   std::vector<int> s2PossibleNewStates;
   for (unsigned i = 0; i < Basis[s2].size(); ++i)
   {
      //      TRACE(Basis[s2][i])(is_projection(Basis[s2][i], s2WantDelta));
      if (is_projection(Basis[s2][i], s2WantDelta))
         s2PossibleNewStates.push_back(i);
   }
   if (s2PossibleNewStates.empty())
      return false; // no possible shift of s2

   // choose one of the states
   int s2NewState = s2PossibleNewStates[rand() % s2PossibleNewStates.size()];

   //   TRACE(Basis[s1][Config.State[s1]])(Basis[s2][Config.State[s2]]);
   //   TRACE(Basis[s1][NewState])(Basis[s2][s2NewState]);

   // now do the actual flips
   std::vector<QuantumNumber> NewHeight = Config.Height;
   if (s1 < s2)
   {
      for (int i = s1+1; i <= s2; ++i)
      {
         //         TRACE("changing s1 < s2")(s1)(s2)(i)(NewHeight[i])(negate(DeltaChange))(is_possible(NewHeight[i], negate(DeltaChange)));
         if (!is_possible(NewHeight[i], negate(DeltaChange))) return false;
         NewHeight[i] = change(NewHeight[i], negate(DeltaChange));
      }
   }
   else
   {
      // s2 < s1
      for (int i = s2+1; i <= s1; ++i)
      {
         //         TRACE("changing s1 > s2")(s1)(s2)(i)(NewHeight[i])(DeltaChange);
         if (!is_possible(NewHeight[i], DeltaChange)) return false;
         NewHeight[i] = change(NewHeight[i], DeltaChange);
      }
   }

   // if we get here, the flip is allowed.  Update the configuration.
   std::swap(Config.Height, NewHeight);
   Config.State[s1] = NewState;
   Config.State[s2] = s2NewState;

   //   Config.CheckValid(Basis);

   return true;
}

double
CalculateMarshallSign(LinearWavefunction const& Psi, WavefunctionDesc& State, int NumFlips, bool Verbose = false)
{
   long x = 0;
   std::vector<BasisList> LocalBasis = ExtractLocalBasis(Psi);
   complex Current = Amplitude(Psi, State);
   int TotalFlips = 0;
   for (int i = 0; i < NumFlips; ++i)
   {
      if (Verbose && ((i+1) % (NumFlips / 10) == 0))
         std::cerr << (i+1) << '\n';

      for (unsigned n = 0; n < Psi.size(); ++n)
      {
         WavefunctionDesc Trial = State;
         int s1 = rand() % Psi.size();
         int s2 = rand() % Psi.size();
         bool Flipped = TryFlipPair(Trial, LocalBasis, s1, s2);
         //         Trial.CheckValid(LocalBasis);
         if (Flipped)
         {
            complex New = Amplitude(Psi, Trial);
            if (norm_2_sq(Current) == 0 || double(rand()) / RAND_MAX < norm_2_sq(New) / norm_2_sq(Current))
            {
               State = Trial;
               Current = New;
               ++TotalFlips;
            }
         }
      }
      int Next = Current.real() > 0 ? 1 : -1;
      x = x + Next;
   }
   if (Verbose)
      std::cerr << "Total flips = " << TotalFlips << '\n';

   return double(x) / NumFlips;
}

double CalculateMarshallSign(LinearWavefunction const& Psi, int BurnIn, int NumFlips, bool Verbose)
{
   std::vector<BasisList> LocalBasis = ExtractLocalBasis(Psi);
   QuantumNumber Target = Psi.TransformsAs();
   if (Verbose)
      std::cerr << "Constructing initial configuration...\n";
   WavefunctionDesc Config = CreateRandomConfiguration(LocalBasis, Target, 3.0);
   if (Verbose)
      std::cerr << "Burning in...\n";
   CalculateMarshallSign(Psi, Config, BurnIn);
   if (Verbose)
      std::cerr << "Calculating statistics...\n";
   double x = CalculateMarshallSign(Psi, Config, NumFlips, Verbose);
   if (Verbose)
      std::cerr << "Done.\n";
   return x;
}

int main(int argc, char** argv)
{
   try
   {
      bool Verbose = false;
      std::string PsiStr;
      int BurnIn = 100;
      int Sweeps = 100;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("burnin,b", prog_opt::value(&BurnIn),
          FormatDefault("Burn in iterations", BurnIn).c_str())
         ("iter,i", prog_opt::value(&Sweeps),
          FormatDefault("Number of iterations for the statistics", Sweeps).c_str())
         ("seed,s", prog_opt::value<unsigned int>(),
          ("Random seed [range 0.."+boost::lexical_cast<std::string>(RAND_MAX)+"]").c_str())
         ("verbose,v", prog_opt::bool_switch(&Verbose),
          "extra debug output")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value<std::string>(&PsiStr), "psi1")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("psi") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-simple-marshall-qmc [options] <psi>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));

      unsigned int RandSeed = vm.count("seed") ? (vm["seed"].as<unsigned long>() % RAND_MAX) : (ext::get_unique() % RAND_MAX);
      srand(RandSeed);
      if (Verbose)
         std::cerr << "Random seed = " << RandSeed << '\n';

      pvalue_ptr<MPWavefunction> P = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);
      LinearWavefunction Psi = *P;

      double x = CalculateMarshallSign(Psi, BurnIn, Sweeps, Verbose);
      double StdError = std::sqrt((1.0 - x*x) / Sweeps);      // this ignores the autocorrelation
      std::cout << x << " +/- " << StdError << '\n';

      pheap::Shutdown();

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
