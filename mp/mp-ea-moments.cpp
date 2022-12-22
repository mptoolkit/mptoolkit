// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ea-moments.cpp
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

#include "wavefunction/mpwavefunction.h"
#include "mp/copyright.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "common/terminal.h"
#include "common/environment.h"
#include "common/prog_options.h"
#include "lattice/infinite-parser.h"
#include "mp-algorithms/triangular_mpo_solver.h"
#include "mp-algorithms/transfer.h"

namespace prog_opt = boost::program_options;

// Tolerance for the overlap of PsiLeft/PsiRight when calculating the leading
// eigenvectors of the mixed transfer matrix: we treat the states as orthogonal
// if the overlap - 1 is greater than this value.
double const OverlapTol = 1e-8;

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      double GMRESTol = 1e-10;
      double K = 0.0;
      std::string InputFilename;
      std::string OpStr;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("gmrestol", prog_opt::value(&GMRESTol), FormatDefault("Error tolerance for the GMRES algorithm", GMRESTol).c_str())
         ("momentum,k", prog_opt::value(&K), "Use this momentum for the EA wavefunction instead of the one in the file (in units of pi)")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value<std::string>(&InputFilename), "wavefunction")
         ("operator", prog_opt::value<std::string>(&OpStr), "operator")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);
      p.add("operator", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("psi") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi> [operator]" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      // Load the wavefunction.
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(InputFilename, CacheSize, true);
      EAWavefunction Psi = PsiPtr->get<EAWavefunction>();

      CHECK(Psi.window_size() == 1);

      // Extract the left and right semi-infinite boundaries.
      InfiniteWavefunctionLeft PsiLeft = Psi.left();
      inplace_qshift(PsiLeft, Psi.left_qshift());
      PsiLeft.rotate_left(Psi.left_index());

      QuantumNumber QShift = PsiLeft.qshift();

      InfiniteWavefunctionRight PsiRight = Psi.right();
      inplace_qshift(PsiRight, Psi.right_qshift());
      PsiRight.rotate_left(Psi.right_index());

      CHECK(PsiRight.qshift() == PsiLeft.qshift());

      LinearWavefunction PsiLinearLeft, PsiLinearRight;
      RealDiagonalOperator LambdaLeft, LambdaRight;
      std::tie(PsiLinearLeft, LambdaLeft) = get_left_canonical(PsiLeft);
      std::tie(LambdaRight, PsiLinearRight) = get_right_canonical(PsiRight);

      MatrixOperator RhoLeft = delta_shift(LambdaLeft*LambdaLeft, QShift);
      MatrixOperator RhoRight = LambdaRight*LambdaRight;

      // Extract the (single-site) windows.
      std::vector<StateComponent> BVec;
      for (WavefunctionSectionLeft Window : Psi.window_vec())
      {
         LinearWavefunction PsiLinear;
         MatrixOperator U;
         std::tie(PsiLinear, U) = get_left_canonical(Window);
         // Note that we assume that the window is single-site.
         BVec.push_back(PsiLinear.get_front()*U);
      }

      // The "triangular" wavefunction containing all of the windows.
      LinearWavefunction PsiTri;

      if (PsiLeft.size() == 1)
         PsiTri.push_back(BVec.back());
      else
      {
         auto CL = PsiLinearLeft.begin();
         auto CR = PsiLinearRight.begin();
         auto B = BVec.begin();
         SumBasis<VectorBasis> NewBasis0((*CL).Basis2(), (*B).Basis2());
         PsiTri.push_back(tensor_row_sum(*CL, *B, NewBasis0));
         ++CL, ++CR, ++B;
         for (int i = 1; i < PsiLeft.size()-1; ++i)
         {
            StateComponent Z = StateComponent((*CL).LocalBasis(), (*CR).Basis1(), (*CL).Basis2());
            SumBasis<VectorBasis> NewBasis1((*CL).Basis2(), (*B).Basis2());
            SumBasis<VectorBasis> NewBasis2((*CL).Basis1(), (*CR).Basis1());
            PsiTri.push_back(tensor_col_sum(tensor_row_sum(*CL, *B, NewBasis1), tensor_row_sum(Z, *CR, NewBasis1), NewBasis2));
            ++CL, ++CR, ++B;
         }
         SumBasis<VectorBasis> NewBasis3((*B).Basis1(), (*CR).Basis1());
         PsiTri.push_back(tensor_col_sum(*B, *CR, NewBasis3));
      }

      // Calculate the eigenvectors of the mixed left/right and right/left transfer operators.
      std::complex<double> OverlapLR, OverlapRL;
      MatrixOperator TLRLeft, TLRRight;
      MatrixOperator TRLLeft, TRLRight;

      std::tie(OverlapLR, TLRLeft, TLRRight) = get_transfer_eigenpair(PsiLinearLeft, PsiLinearRight, PsiLeft.qshift());
      TLRRight = delta_shift(TLRRight, QShift);
      if (Verbose > 1)
         std::cout << "LR overlap = " << OverlapLR << std::endl;

      std::tie(OverlapRL, TRLLeft, TRLRight) = get_transfer_eigenpair(PsiLinearRight, PsiLinearLeft, PsiLeft.qshift());
      TRLRight = delta_shift(TRLRight, QShift);
      if (Verbose > 1)
         std::cout << "RL overlap = " << OverlapRL << std::endl;

      if (std::abs(OverlapLR - std::conj(OverlapRL)) > OverlapTol)
         PANIC("OverlapLR and conj(OverlapRL) are different!");

      // If the left and right boundaries are different states, then we do not need these eigenvectors.
      if (std::abs(std::abs(OverlapLR) - 1.0) > OverlapTol)
      {
         TLRLeft = MatrixOperator();
         TLRRight = MatrixOperator();
         TRLLeft = MatrixOperator();
         TRLRight = MatrixOperator();
      }

      // If the operator is not specified, use the one specified in the file attributes.
      if (!vm.count("operator"))
      {
         OpStr = PsiPtr->Attributes()["Hamiltonian"].get_or_default(std::string());
         if (OpStr.empty())
         {
            std::cerr <<  basename(argv[0]) << ": fatal: no operator specified, and wavefunction "
               "attribute Hamiltonian does not exist or is empty." << std::endl;
            return 1;
         }
      }

      // Load the operator.
      BasicTriangularMPO Op;
      InfiniteLattice Lattice;
      std::tie(Op, Lattice) = ParseTriangularOperatorAndLattice(OpStr);

      // Ensure Op is the correct size.
      if (Op.size() < PsiLeft.size())
         Op = repeat(Op, PsiLeft.size() / Op.size());
      CHECK_EQUAL(Op.size(), PsiLeft.size());

      int LatticeUCsPerPsiUC = PsiLeft.size() / Lattice.GetUnitCell().size();

      // Use the phase factor from the input file by default, otherwise, use the specified value.
      std::complex<double> ExpIK;
      if (!vm.count("momentum"))
         ExpIK = Psi.exp_ik();
      else
         ExpIK = exp(std::complex<double>(0.0, math_const::pi) * (K * LatticeUCsPerPsiUC));

      // Calculate the full E-matrix.
      double UnityEpsilon = DefaultEigenUnityEpsilon;

      std::vector<MatrixPolyType> EMat0;
      MatrixOperator IdentLeft = MatrixOperator::make_identity(PsiLeft.Basis1());
      SolveSimpleMPO_Left(EMat0, PsiLinearLeft, QShift, Op, IdentLeft, RhoLeft,
                          true, 0, GMRESTol, UnityEpsilon, Verbose);

      std::vector<KMatrixPolyType> EMatKTop;
      SolveSimpleMPO_EA_Left(EMatKTop, EMat0, std::vector<KMatrixPolyType>(), std::vector<KMatrixPolyType>(),
                             PsiLinearLeft, PsiLinearRight, PsiTri,
                             QShift, Op, TRLLeft, TRLRight, ExpIK,
                             true, 0, GMRESTol, UnityEpsilon, Verbose, "top");

      std::vector<KMatrixPolyType> EMatKBot;
      SolveSimpleMPO_EA_Left(EMatKBot, EMat0, std::vector<KMatrixPolyType>(), std::vector<KMatrixPolyType>(),
                             PsiLinearLeft, PsiLinearRight, PsiTri,
                             QShift, Op, TLRLeft, TLRRight, ExpIK,
                             true, 0, GMRESTol, UnityEpsilon, Verbose, "bottom");

      std::vector<KMatrixPolyType> EMatK1;
      MatrixOperator IdentRight = MatrixOperator::make_identity(PsiLinearRight.Basis1());
      SolveSimpleMPO_EA_Left(EMatK1, EMat0, EMatKTop, EMatKBot,
                             PsiLinearLeft, PsiLinearRight, PsiTri,
                             QShift, Op, RhoRight, IdentRight, ExpIK,
                             false, 0, GMRESTol, UnityEpsilon, Verbose, "final");

      // Print the expectation values.
      std::cout << "#degree "
                << "#real                   "
                << "#imag                   "
                << std::endl << std::left;

      for (auto I = EMatK1.back()[1.0].begin(); I != EMatK1.back()[1.0].end(); ++I)
      {
         std::complex<double> x = inner_prod(I->second, IdentRight);
         std::cout << std::setw(7) << I->first << " "
                   << std::setw(20) << std::real(x) << "    "
                   << std::setw(20) << std::imag(x) << std::endl;
      }

      pheap::Shutdown();
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
