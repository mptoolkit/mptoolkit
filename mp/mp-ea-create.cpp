// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ea-create.cpp
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

#include "common/environment.h"
#include "common/formatting.h"
#include "common/prog_opt_accum.h"
#include "common/prog_options.h"
#include "common/terminal.h"
#include "common/unique.h"
#include "interface/inittemp.h"
#include "lattice/unitcell-parser.h"
#include "lattice/infinitelattice.h"
#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string LeftFilename;
      std::string RightFilename;
      std::string OpStr;
      double K = 0.0;
      std::string OutputFilename;
      bool Streaming = false;
      bool NoStreaming = false;
      bool Force = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("momentum,k", prog_opt::value(&K), FormatDefault("The momentum (in units of pi)", K).c_str())
         ("output,o", prog_opt::value(&OutputFilename), "Output filename")
         ("streaming", prog_opt::bool_switch(&Streaming), "Store the left and right strips by reference to the input files")
         ("no-streaming", prog_opt::bool_switch(&NoStreaming), "Store the left and right strips into the output file [default]")
         ("force,f", prog_opt::bool_switch(&Force), "Force overwriting output file")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("op", prog_opt::value(&OpStr), "op")
         ("psi", prog_opt::value(&LeftFilename), "psi")
         ("psi2", prog_opt::value(&RightFilename), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("op", 1);
      p.add("psi", 1);
      p.add("psi2", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("psi") == 0 || vm.count("output") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <operator> <psi> [psi-right] -o <psi-out>" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      if (Streaming && NoStreaming)
      {
         std::cerr << "fatal: Cannot use --streaming and --no-streaming simultaneously!" << std::endl;
         return 1;
      }
      else if (!Streaming && !NoStreaming)
         NoStreaming = true; // This is the current default behavior.

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      pheap::Initialize(OutputFilename, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);

      pvalue_ptr<MPWavefunction> InPsiLeft = pheap::ImportHeap(LeftFilename);
      InfiniteWavefunctionLeft PsiLeft = InPsiLeft->get<InfiniteWavefunctionLeft>();

      InfiniteWavefunctionRight PsiRight;
      if (vm.count("psi2"))
      {
         pvalue_ptr<MPWavefunction> InPsiRight = pheap::ImportHeap(RightFilename);
         if (InPsiRight->is<InfiniteWavefunctionRight>())
            PsiRight = InPsiRight->get<InfiniteWavefunctionRight>();
         else if (InPsiRight->is<InfiniteWavefunctionLeft>())
         {
            if (Streaming)
            {
               std::cerr << "fatal: psi-right must be an InfiniteWavefunctionRight if streaming is enabled." << std::endl;
               return 1;
            }

            PsiRight = InfiniteWavefunctionRight(InPsiRight->get<InfiniteWavefunctionLeft>());
         }
         else
         {
            std::cerr << "fatal: psi-right must be an InfiniteWavefunctionLeft or InfiniteWavefunctionRight." << std::endl;
            return 1;
         }
      }
      else
      {
         // This condition could be relaxed, in that case we would only save
         // the left boundary by reference, but we assume that the user want
         // both boundaries saved by reference.
         if (Streaming)
         {
            std::cerr << "fatal: psi-right must be specified if streaming is enabled." << std::endl;
            return 1;
         }

         PsiRight = InfiniteWavefunctionRight(PsiLeft);
      }

      if (PsiLeft.size() > 1 || PsiRight.size() > 1)
      {
         // TODO
         std::cerr << "fatal: multisite unit cells are not implemented yet." << std::endl;
         return 1;
      }

      InfiniteLattice Lattice;
      UnitCellMPO Op;
      std::tie(Op, Lattice) = ParseUnitCellOperatorAndLattice(OpStr);

      if (Op.size() > 1)
      {
         // TODO
         std::cerr << "fatal: multisite operators are not implemented yet." << std::endl;
         return 1;
      }

      WavefunctionSectionLeft PsiWindow(PsiLeft);

      LinearWavefunction PsiWindowLinear;
      MatrixOperator Lambda;

      std::tie(PsiWindowLinear, Lambda) = get_left_canonical(PsiWindow);
      PsiWindowLinear.set_back(prod(PsiWindowLinear.get_back(), Lambda));

      auto J = Op.MPO().begin();
      for (auto I = PsiWindowLinear.begin(); I != PsiWindowLinear.end(); ++I, ++J)
         (*I) = aux_tensor_prod(*J, *I);

      MatrixOperator Identity = MatrixOperator::make_identity(PsiWindowLinear.Basis2());
      PsiWindow = WavefunctionSectionLeft::ConstructFromLeftOrthogonal(std::move(PsiWindowLinear), Identity, Verbose);

      EAWavefunction PsiEA(PsiLeft, std::vector<WavefunctionSectionLeft>(1, PsiWindow), PsiRight,
                           Op.qn1(), Op.qn2(), 0, 0, std::exp(std::complex<double>(0.0, math_const::pi) * K));

      if (Streaming)
      {
         PsiEA.set_left_filename(LeftFilename);
         PsiEA.set_right_filename(RightFilename);
      }

      MPWavefunction Wavefunction;
      Wavefunction.Wavefunction() = std::move(PsiEA);
      Wavefunction.AppendHistoryCommand(EscapeCommandline(argc, argv));
      Wavefunction.SetDefaultAttributes();

      pvalue_ptr<MPWavefunction> PsiPtr(new MPWavefunction(Wavefunction));
      pheap::ShutdownPersistent(PsiPtr);
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
