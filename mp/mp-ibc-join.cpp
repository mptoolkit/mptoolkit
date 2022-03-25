// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ibc-join.cpp
//
// Copyright (C) 2020 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "lattice/latticesite.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp-algorithms/lanczos.h"
#include "lattice/infinitelattice.h"
#include "lattice/infinite-parser.h"
#include "mp-algorithms/triangular_mpo_solver.h"

namespace prog_opt = boost::program_options;

struct CMultiply
{
   CMultiply(StateComponent const& E_, StateComponent const& F_)
      : E(E_), F(F_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& C) const
   {
      return operator_prod(E, C, herm(F));
   }

   StateComponent const& E;
   StateComponent const& F;
};

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string OutputFile;
      std::string InputFileLeft;
      std::string InputFileRight;
      std::string HamStr;
      bool Random = false;
      bool Force = false;
      double GMRESTol = 1E-13;    // tolerance for GMRES for the initial H matrix elements.
      int Iter = 50;
      int MinIter = 4;
      double Tol = 1E-16;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("left,l", prog_opt::value(&InputFileLeft),
          "input iMPS wavefunction for the left semi-infinite strip [required]")
         ("right,r", prog_opt::value(&InputFileRight),
          "input iMPS wavefunction for the right semi-infinite strip [required]")
         ("output,o", prog_opt::value(&OutputFile), "output IBC wavefuction [required]")
         ("force,f", prog_opt::bool_switch(&Force), "overwrite output files")
         ("random", prog_opt::bool_switch(&Random), "Don't minimize the Hamiltonian, make a random state instead")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "model Hamiltonian, of the form lattice:operator")
         ("maxiter", prog_opt::value<int>(&Iter),
          FormatDefault("Maximum number of Lanczos iterations", Iter).c_str())
         ("miniter", prog_opt::value<int>(&MinIter),
          FormatDefault("Minimum number of Lanczos iterations", MinIter).c_str())
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Error tolerance for the Lanczos eigensolver", Tol).c_str())
         ("gmrestol", prog_opt::value(&GMRESTol),
          FormatDefault("Error tolerance for the GMRES algorithm", GMRESTol).c_str())
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
         ;

      prog_opt::options_description hidden("Hidden options");

      prog_opt::positional_options_description p;

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || InputFileLeft.empty() || InputFileRight.empty() || OutputFile.empty())
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] -l <left_psi> -r <right_psi> -o <output_psi>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      pheap::Initialize(OutputFile, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);

      pvalue_ptr<MPWavefunction> InPsiLeft = pheap::ImportHeap(InputFileLeft);
      pvalue_ptr<MPWavefunction> InPsiRight = pheap::ImportHeap(InputFileRight);

      InfiniteWavefunctionLeft PsiLeft = InPsiLeft->get<InfiniteWavefunctionLeft>();

      // There are some situations where the first method does not work
      // properly, so temporarily use the second method as a workaround.
#if 0
      InfiniteWavefunctionRight PsiRight = InPsiRight->get<InfiniteWavefunctionLeft>();
#else
      InfiniteWavefunctionLeft PsiRightLeft = InPsiRight->get<InfiniteWavefunctionLeft>();

      MatrixOperator U;
      RealDiagonalOperator D;
      LinearWavefunction PsiRightLinear;
      std::tie(U, D, PsiRightLinear) = get_right_canonical(PsiRightLeft);

      InfiniteWavefunctionRight PsiRight(U*D, PsiRightLinear, PsiRightLeft.qshift());
#endif

      // If the bases of the two boundary unit cells have only one quantum
      // number sector, manually ensure that they match.
      // FIXME: This workaround probably will not work for non-Abelian symmetries.
      if (PsiLeft.Basis2().size() == 1 && PsiRight.Basis1().size() == 1)
         inplace_qshift(PsiRight, delta_shift(PsiLeft.Basis2()[0], adjoint(PsiRight.Basis1()[0])));

      WavefunctionSectionLeft Window;

      // To join the wavefunctions we need a matrix to represent the wavefunction in the
      // left/right bipartite basis. We can do that by minimizing the energy, or we
      // can make a random state. We need a random state anyway, as input for the Lanczos
      MatrixOperator C = MakeRandomMatrixOperator(PsiLeft.Basis2(), PsiRight.Basis1());
      C *= 1.0 / norm_frob(C);
      if (Random)
      {
         // we don't need to do anything here
      }
      else if (!HamStr.empty())
      {
         InfiniteLattice Lattice;
         BasicTriangularMPO HamMPO, HamMPOLeft, HamMPORight;
         std::tie(HamMPO, Lattice) = ParseTriangularOperatorAndLattice(HamStr);

         HamMPOLeft = HamMPO;
         if (HamMPOLeft.size() < PsiLeft.size())
            HamMPOLeft = repeat(HamMPOLeft, PsiLeft.size() / HamMPOLeft.size());

         std::cout << "Solving fixed-point Hamiltonian..." << std::endl;
         StateComponent BlockHamL = Initial_E(HamMPOLeft, PsiLeft.Basis1());
         std::complex<double> LeftEnergy = SolveSimpleMPO_Left(BlockHamL, PsiLeft, HamMPOLeft,
                                                               GMRESTol, Verbose);
         std::cout << "Starting energy (left eigenvalue) = " << LeftEnergy << std::endl;
         BlockHamL = delta_shift(BlockHamL, adjoint(PsiLeft.qshift()));

         HamMPORight = HamMPO;
         if (HamMPORight.size() < PsiRight.size())
            HamMPORight = repeat(HamMPORight, PsiRight.size() / HamMPORight.size());

         StateComponent BlockHamR = Initial_F(HamMPORight, PsiRight.Basis2());
         std::complex<double> RightEnergy = SolveSimpleMPO_Right(BlockHamR, PsiRight, HamMPORight,
                                                                 GMRESTol, Verbose);
         std::cout << "Starting energy (right eigenvalue) = " << RightEnergy << std::endl;
         BlockHamR = delta_shift(BlockHamR, PsiRight.qshift());

         double Energy = Lanczos(C, CMultiply(BlockHamL, BlockHamR),
                                 Iter, Tol, MinIter, Verbose);

         C *= 1.0 / norm_frob(C);

         std::cout << "Energy is " << Energy << '\n';
      }
      else
      {
         std::cerr << basename(argv[0]) << ": fatal: no Hamiltonian specified\n";
         return 1;
      }

      // Make the window from our centre matrix
      Window = WavefunctionSectionLeft(C);

      IBCWavefunction ResultPsi(PsiLeft, Window, PsiRight, 0);

      MPWavefunction Result(ResultPsi);

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
