// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ea-change-k.cpp
//
// Copyright (C) 2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "wavefunction/mpwavefunction.h"
#include "mp/copyright.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "common/terminal.h"
#include "common/environment.h"
#include "common/prog_options.h"
#include "interface/inittemp.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      double K = 0.0;
      int LatticeUCSize = 1;
      bool Force = false;
      std::string InputFilename;
      std::string OutputFilename;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("momentum,k", prog_opt::value(&K), "Change the momentum to this value (in units of pi)")
         ("latticeucsize", prog_opt::value(&LatticeUCSize), "Lattice unit cell size [default wavefunction attribute \"LatticeUnitCellSize\" or 1]")
         ("force,f", prog_opt::bool_switch(&Force), "Force overwriting output file")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi-in", prog_opt::value(&InputFilename), "Input wavefunction")
         ("psi-out", prog_opt::value(&OutputFilename), "Output wavefunction")
         ;

      prog_opt::positional_options_description p;
      p.add("psi-in", 1);
      p.add("psi-out", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("psi-out") == 0 || vm.count("momentum") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi-in> <psi-out> -k <momentum>" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      // Load the wavefunction.
      pvalue_ptr<MPWavefunction> PsiPtr;
      if (InputFilename == OutputFilename)
         PsiPtr = pheap::OpenPersistent(InputFilename.c_str(), mp_pheap::CacheSize());
      else
      {
         pheap::Initialize(OutputFilename, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);
         PsiPtr = pheap::ImportHeap(InputFilename);
      }

      // Get the lattice unit cell size if unspecified.
      if (!vm.count("latticeucsize"))
         LatticeUCSize = PsiPtr->Attributes()["LatticeUnitCellSize"].get_or_default<int>(1);

      EAWavefunction Psi = PsiPtr->get<EAWavefunction>();

      int LatticeUCsPerPsiUC = Psi.left().size() / LatticeUCSize;
      std::complex<double> ExpIK = exp(std::complex<double>(0.0, math_const::pi) * (K * LatticeUCsPerPsiUC));

      EAWavefunction PsiNew(Psi.left(), Psi.window_vec(), Psi.right(), Psi.left_qshift(), Psi.right_qshift(),
                            Psi.left_index(), Psi.right_index(), ExpIK, Psi.gs_overlap());

      // Stream the boundaries if the input file does.
      PsiNew.set_left_filename(Psi.get_left_filename());
      PsiNew.set_right_filename(Psi.get_right_filename());

      PsiPtr.mutate()->Wavefunction() = PsiNew;
      PsiPtr.mutate()->AppendHistoryCommand(EscapeCommandline(argc, argv));
      PsiPtr.mutate()->SetDefaultAttributes();

      pheap::ShutdownPersistent(PsiPtr);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << std::endl;
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
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!" << std::endl;
      return 1;
   }
}
