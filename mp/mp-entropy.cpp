// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-entropy.cpp
//
// Copyright (C) 2025 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "common/environment.h"
#include "common/formatting.h"
#include "common/prog_opt_accum.h"
#include "common/prog_options.h"
#include "common/terminal.h"
#include "interface/inittemp.h"
#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string Filename;
      int MaxEigenvalues =-1;
      bool Base2 = false;
      bool ShowDensity = false;
      bool ShowDegen = false;
      bool Quiet = false;
      std::vector<int> Sites;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("density-matrix,d", prog_opt::bool_switch(&ShowDensity), "Show the eigenspectrum of the density matrix")
         ("degen", prog_opt::bool_switch(&ShowDegen), "Show degeneracies in the density matrix as repeated eigenvalues (implies -d)")
         ("limit,l", prog_opt::value<int>(&MaxEigenvalues), "Limit the density matrix display to N eigenvalues (implies -d)")
         ("base2,2", prog_opt::bool_switch(&Base2), "Show the entropy using base 2 instead of base e")
         ("quiet", prog_opt::bool_switch(&Quiet), "Do not show column headings")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value(&Filename), "psi")
         ("sites", prog_opt::value<std::vector<int>>(&Sites)->multitoken(), "sites")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);
      p.add("sites", -1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("sites") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi> <sites>" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (vm.count("limit"))
         ShowDensity = true;

      if (ShowDegen)
         ShowDensity = true;

      pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(Filename, mp_pheap::CacheSize(), true);

      FiniteWavefunctionLeft Psi = PsiPtr->get<FiniteWavefunctionLeft>();

      SimpleOperator Rho = make_vacuum_simple(Psi.GetSymmetryList());
      StateComponent C;
      int nPrev;
      bool FirstIter = true;
      for (auto const& n : Sites)
      {
         if (FirstIter)
         {
            C = Psi[n];
            FirstIter = false;
         }
         else
         {
            if (n != nPrev + 1)
            {
               C = prod(C, Psi.lambda(nPrev+1));
               Rho = tensor_prod(Rho, trace_prod(C, herm(C)));

               C = Psi[n];
            }
            else
            {
               C = local_tensor_prod(C, Psi[n]);
            }
         }
         nPrev = n;
      }

      C = prod(C, Psi.lambda(nPrev+1));
      Rho = tensor_prod(Rho, trace_prod(C, herm(C)));

      DensityMatrix<SimpleOperator> DM(Rho);
      if (!ShowDensity)
      {
         // Just print the entropy
         std::cout << DensityEntropy(DM.begin(), DM.end(), DM.EigenSum(), Base2) << std::endl;
      }
      else
      {
         // Print the full density matrix eigenspectrum
         DM.DensityMatrixReport(std::cout, MaxEigenvalues, Base2, ShowDegen, Quiet);
         std::cout << std::endl;
      }

      pheap::Shutdown();
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << std::endl;
      pheap::Cleanup();
      return 1;
   }
   catch (pheap::PHeapCannotCreateFile& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      if (e.Why == "File exists")
         std::cerr << "Note: Use --force (-f) option to overwrite." << std::endl;
      return 1;
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      pheap::Cleanup();
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!" << std::endl;
      pheap::Cleanup();
      return 1;
   }
}
