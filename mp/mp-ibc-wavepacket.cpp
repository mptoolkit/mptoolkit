// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ibc-wavepacket.cpp
//
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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
#include "lattice/infinite-parser.h"
#include "lattice/infinitelattice.h"
#include "linearalgebra/arpack_wrapper.h"
#include "mp-algorithms/excitation-ansatz.h"
#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      double KMax = 0;
      double KMin = 0;
      int KNum = 1;
      std::string InputPrefix;
      int InputDigits = 0;
      std::string OutputFilename;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("kmax,k", prog_opt::value(&KMax), FormatDefault("Maximum momentum (divided by pi)", KMax).c_str())
         ("kmin", prog_opt::value(&KMin), FormatDefault("Minimum momentum (divided by pi)", KMin).c_str())
         ("knum", prog_opt::value(&KNum), "Number of momentum steps to use")
         ("wavefunction,w", prog_opt::value(&InputPrefix), "Prefix for input filenames (of the form [prefix].k[k])")
         ("output,o", prog_opt::value(&OutputFilename), "Output filename")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help"))
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options]" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      mp_pheap::InitializeTempPHeap();

      double KStep = (KMax-KMin)/(KNum-1);
      InputDigits = std::max(formatting::digits(KMax), formatting::digits(KStep));

      // Read input wavefunctions.
      std::vector<EAWavefunction> PsiVec;
      for (int n = 0; n < KNum; ++n)
      {
         std::string InputFilename = InputPrefix + ".k" + formatting::format_digits(KMin + KStep*n, InputDigits);
         pvalue_ptr<MPWavefunction> InPsi = pheap::ImportHeap(InputFilename);
         PsiVec.push_back(InPsi->get<EAWavefunction>());
         // We only handle single-site EAWavefunctions at the moment.
         CHECK(PsiVec.back().window_size() == 1);
      }
      
      std::vector<std::vector<StateComponent>> BSymVec;
      for (EAWavefunction Psi : PsiVec)
      {
         BSymVec.push_back(std::vector<StateComponent>());

         // Get the right null space matrices corresponding to each A-matrix in PsiRight.
         LinearWavefunction PsiLinearLeft;
         std::tie(PsiLinearLeft, std::ignore) = get_left_canonical(Psi.Left);

         LinearWavefunction PsiLinearRight;
         std::tie(std::ignore, PsiLinearRight) = get_right_canonical(Psi.Right);

         std::vector<StateComponent> NullRightVec;
         for (StateComponent C : PsiLinearRight)
            NullRightVec.push_back(NullSpace1(C));

         auto NR = NullRightVec.begin();
         auto AL = PsiLinearLeft.begin();
         auto AR = PsiLinearRight.begin();
         for (WavefunctionSectionLeft Window : Psi.WindowVec)
         {
            LinearWavefunction PsiLinear; 
            MatrixOperator U;
            std::tie(PsiLinear, U) = get_left_canonical(Window);
            // Note that we assume that the window is single-site.
            //BSymVec.back().push_back(PsiLinear.get_front()*U);
            StateComponent BL = PsiLinear.get_front()*U;

            // Find the B-matrix satisfying the right-gauge fixing condition.
            MatrixOperator XR = scalar_prod(BL, herm(*NR));
            StateComponent BR = prod(XR, *NR);
            // Scale norm to match BL
            //BR *= norm_frob(BL) / norm_frob(BR);

            TRACE(inner_prod(BR, *AR))(inner_prod(BR, *AL));
            TRACE(inner_prod(BL, *AL))(inner_prod(BL, *AR));
            TRACE(inner_prod(BR, BL));
            TRACE(norm_frob(BL))(norm_frob(BR))(norm_frob(XR));

            ++NR, ++AR, ++AL;
         }
      }
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
