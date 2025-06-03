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

#include "common/environment.h"
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

      // Extract unit cell and lambda matrix of infinite wave function
      LinearWavefunction PsiUC;
      RealDiagonalOperator Lambda;

      std::tie(PsiUC, Lambda) = get_left_canonical(PsiInf);

      // Find degree 1 quantum number in right basis with largest size
      QuantumNumber QR;
      int Max = 0;
      for (int i = 0; i < PsiUC.Basis2().size(); ++i)
      {
         if (PsiUC.Basis2()[i].degree() == 1)
         {
            if (PsiUC.Basis2().dim(i) > Max)
            {
               QR = PsiUC.Basis2()[i];
               Max = PsiUC.Basis2().dim(i);
            }
         }
      }

      if (Max == 0)
      {
         std::cerr << "fatal: No quantum number of the unit cell basis has degree 1: cannot create finite wave function.\n";
         return 1;
      }

      if (Verbose > 0)
         std::cout << "Setting right basis quantum number " << QR << '\n';

      //QR = adjoint(QR);

      // Output finite wave function
      LinearWavefunction Psi;

      // Left basis quantum number
      QuantumNumber QL = QR;

      for (int i = 0; i < Repeat; ++i)
      {
         auto I = PsiUC.end();
         while (I != PsiUC.begin())
         {
            --I;
            Psi.push_front(delta_shift(*I, QL));
         }
         QL = delta_shift(QL, QShift);
      }

      // Choose random left/right vectors to make the ends of the MPS one-dimensional.
      VectorBasis Vacuum(make_vacuum_basis(Psi.Basis2().GetSymmetryList()));

      MatrixOperator VR = MakeRandomMatrixOperator(Psi.Basis2(), Vacuum);
      MatrixOperator VL = MakeRandomMatrixOperator(delta_shift(delta_shift(Vacuum, QL), adjoint(QR)), Psi.Basis1());

      Psi.set_front(VL * delta_shift(Lambda, QL) * Psi.get_front());
      Psi.set_back(Psi.get_back() * VR);

      // Canonicalize output wave function
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
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      return 1;
   }
   catch (pheap::PHeapCannotCreateFile& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      if (e.Why == "File exists")
         std::cerr << "Note: Use --force (-f) option to overwrite.\n";
      return 1;
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
