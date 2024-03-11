// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-evolve-bonds.cpp
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
#include "matrixproduct/mpoperator.h"
#include "matrixproduct/density.h"
#include "quantumnumbers/all_symmetries.h"
#include "mp-algorithms/local-evolution.h"
#include "tensor/tensor_exponential.h"
#include "pheap/pheap.h"
#include "mp-algorithms/dmrg.h"
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include <iostream>
#include "common/environment.h"

#include "models/spin.h"
#include "models/spin-u1.h"
#include "models/spin-su2.h"
#include "models/spinlessfermion-u1.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string Hamiltonian = "";
      int MinStates = 1;
      int MaxStates = 100000;
      double TruncCutoff = 0;
      double EigenCutoff = -1;
      bool Relative = false;
      half_int Spin = 0.5;
      double Lambda = 1.0;
      double Jx = 1.0;
      double Jy = 1.0;
      double Jz = 1.0;
      double V = 0.0;
      std::vector<std::complex<double> > Timesteps;
      bool EvolveAll = false;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value(&Hamiltonian),
          "operator to use for the Hamiltonian, not an expression!"
          " Valid choices: itf, xyz, xxz-u1, xxx-su2")
         ("wavefunction,w", prog_opt::value<std::string>(),
          "wavefunction to apply DMRG (required)")
         ("max-states,m", prog_opt::value<int>(&MaxStates),
          ("Maximum number of states to keep [default "
          +boost::lexical_cast<std::string>(MaxStates)+"]").c_str())
         ("min-states", prog_opt::value<int>(&MinStates), ("Minimum number of states to keep [default "
          +boost::lexical_cast<std::string>(MinStates)+"]").c_str())
         ("trunc,r", prog_opt::value<double>(&TruncCutoff),
          "Truncation error cutoff [default 0]")
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
          ("Cutoff threshold for density matrix eigenvalues (alternative to truncation error) [default "
           +boost::lexical_cast<std::string>(EigenCutoff)+"]").c_str())
         ("relative", prog_opt::bool_switch(&Relative), "truncate relative weights, rather than"
          " absolute (for imaginary time evolution you surely want this)")
         ("all,a", prog_opt::bool_switch(&EvolveAll), "evolve all bonds at once [PRA 60, 1956 (1999)]")
         ("spin", prog_opt::value(&Spin), "spin (for xxx,xxz,xyz hamiltonians)")
         ("jx", prog_opt::value(&Jx), "J_x (for xyz hamiltonian) [default 1.0]")
         ("jy", prog_opt::value(&Jy), "J_y (for xyz hamiltonian) [default 1.0]")
         ("jz", prog_opt::value(&Jz), "J_z (for xxz,xyz hamiltonians) [default 1.0]")
         ("V", prog_opt::value(&V), "nearest-neighbor coulomb (for spinless fermions) [default 0]")
         ("lambda", prog_opt::value(&Lambda), "transverse field strength"
          " (for itf hamiltonian) [default 1.0])")
          ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("timestep", prog_opt::value(&Timesteps), "timestep")
         ;

      prog_opt::positional_options_description p;
      p.add("timestep", -1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0 || vm.count("Hamiltonian") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-evolve-bonds [options] timestep1 timestep2 .... \n";
         std::cerr << desc << '\n';
         return 1;
      }

      if (Timesteps.size() == 0 || Timesteps.size() % 2 != 0)
      {
         std::cerr << "mp-evolve-bonds: error: must have an even number of timesteps.\n";
         TRACE(Timesteps.size());
         return 1;
      }

      SimpleOperator Ham;
      if (Hamiltonian == "xyz")
      {
         std::cout << "Hamiltonian is XYZ model with spin S=" << Spin
                   << "\nJx=" << Jx << ", Jy= " << Jy << ", Jz=" << Jz << '\n';
         SiteBlock Site = CreateSpinSite(Spin);
         Ham = Jx * tensor_prod(Site["Sx"], Site["Sx"])
            + Jy * tensor_prod(Site["Sy"], Site["Sy"])
            + Jz * tensor_prod(Site["Sz"], Site["Sz"]);
      }
      else if (Hamiltonian == "xxz-u1")
      {
         std::cout << "Hamiltonian is XXZ model with spin S=" << Spin
                   << ", Jz=" << Jz << '\n';
         SiteBlock Site = CreateU1SpinSite(Spin);
         Ham = 0.5 * (tensor_prod(Site["Sp"], Site["Sm"])
                      + tensor_prod(Site["Sm"], Site["Sp"]))
            + Jz * tensor_prod(Site["Sz"], Site["Sz"]);
      }
      else if (Hamiltonian == "xxx-su2")
      {
         std::cout << "Hamiltonian is XXX model with spin S=" << Spin << '\n';
         SiteBlock Site = CreateSU2SpinSite(Spin);
         Ham = -sqrt(3.0) * tensor_prod(Site["S"], Site["S"], QuantumNumber(Site.GetSymmetryList()));
      }
      else if (Hamiltonian == "itf")
      {
         std::cout << "Hamiltonian is transverse-field ising.\n";
         SiteBlock Site = CreateSpinSite(0.5);
         Ham = 4.0 * tensor_prod(Site["Sz"], Site["Sz"])
            + Lambda * (tensor_prod(Site["Sx"], Site["I"]) + tensor_prod(Site["I"], Site["Sx"]));
      }
      else if (Hamiltonian == "sf")
      {
         std::cout << "Hamiltonian is spinless fermions.\n";
         SiteBlock Site = CreateU1SpinlessFermion();
         QuantumNumbers::QNConstructor<QuantumNumbers::U1> QN(Site.GetSymmetryList());
         Ham = -1.0 * (tensor_prod(prod(Site["CH"], Site["P"], QN(1)), Site["C"])
                       - tensor_prod(prod(Site["C"], Site["P"], QN(-1)), Site["CH"]))
            + V * tensor_prod(Site["N"], Site["N"]);
      }
      else
      {
         std::cerr << "mp-evolve-bonds: error: Hamiltonian parameter must be one of xyz, xxz-u1, xxx-su2, itf, sf\n";
         exit(1);
      }

      std::cout << "Starting evolution...\n";
      std::string InputWavefunction = vm["wavefunction"].as<std::string>();
      std::cout << "Input wavefunction: " << InputWavefunction << std::endl;

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = MaxStates;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;
      SInfo.TruncateRelative = Relative;

      // Open the wavefunction
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<MPWavefunction> P = pheap::OpenPersistent(InputWavefunction, CacheSize);
      LinearWavefunction Psi = *P;
      P = pvalue_ptr<MPWavefunction>();

      unsigned Size = Psi.size();

      std::complex<double> complex_i(0.0, -1.0);
      for (unsigned i = 0; i < Timesteps.size(); ++i)
      {
         std::list<SimpleOperator> EvolList;
         for (unsigned s = 0; s < Size-1; ++s)
         {
            if ((s % 2 == i % 2 || EvolveAll) && Timesteps[i] != 0.0)
               EvolList.push_back(Exponentiate(complex_i*Timesteps[i]*Ham));
            else
               EvolList.push_back(SimpleOperator::make_identity(Ham.Basis1()));
         }

         if (i % 2 == 0)
            SweepRightEvolve(Psi, EvolList, SInfo, true);
         else
            SweepLeftEvolve(Psi, EvolList, SInfo, true);
      }

      std::cout << "Finished." << std::endl;
      P = pvalue_ptr<MPWavefunction>(new MPWavefunction(Psi));
      pheap::ShutdownPersistent(P);

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
