// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-finite-create.cpp
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

#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"
#include "lattice/latticesite.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "lattice/infinitelattice.h"
#include "lattice/infinite-parser.h"
#include "mp-algorithms/triangular_mpo_solver.h"
#include "linearalgebra/arpack_wrapper.h"
#include "mps/packunpack.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string OutputFilename;
      std::string InputFilename;
      int Repeat = 1;
      bool Force = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("repeat,r", prog_opt::value(&Repeat), "Repeat the infinite unit cell this number of times")
         ("force,f", prog_opt::bool_switch(&Force), "Force overwriting output file")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value(&InputFilename), "psi")
         ("psi-out", prog_opt::value(&OutputFilename), "psi-out")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);
      p.add("psi-out", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("psi-out") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi> <psi-out>" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      pheap::Initialize(OutputFilename, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);

      pvalue_ptr<MPWavefunction> InPsi = pheap::ImportHeap(InputFilename);
      InfiniteWavefunctionLeft PsiInf = InPsi->get<InfiniteWavefunctionLeft>();

      QuantumNumber QShift = PsiInf.qshift();

      LinearWavefunction Psi;
      RealDiagonalOperator Lambda;

      std::tie(Psi, Lambda) = get_left_canonical(PsiInf);

      QuantumNumber Q = QShift;

      if (Repeat > 1)
      {
         LinearWavefunction PsiUC = Psi;

         for (int i = 1; i < Repeat; ++i)
         {
            auto I = PsiUC.end();
            while (I != PsiUC.begin())
            {
               --I;
               Psi.push_front(delta_shift(*I, Q));
            }
            Q = delta_shift(Q, QShift);
         }
      }

      // Choose random left/right vectors to make the ends of the MPS one-dimensional.
      VectorBasis Vacuum(make_vacuum_basis(Psi.Basis2().GetSymmetryList()));

      MatrixOperator VR = MakeRandomMatrixOperator(Psi.Basis2(), Vacuum); 
      MatrixOperator VL = MakeRandomMatrixOperator(delta_shift(Vacuum, Q), Psi.Basis1()); 

      Psi.set_front(VL * delta_shift(Lambda, Q) * Psi.get_front());
      Psi.set_back(Psi.get_back() * VR);

      FiniteWavefunctionLeft PsiOut = FiniteWavefunctionLeft::Construct(Psi);
      PsiOut.check_structure();

      normalize(PsiOut);

      MPWavefunction Result(PsiOut);

      // Attributes
      Result.SetDefaultAttributes();

      // History log
      Result.AppendHistoryCommand(EscapeCommandline(argc, argv));

      pvalue_ptr<MPWavefunction> P(new MPWavefunction(Result));
      pheap::ShutdownPersistent(P);
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
