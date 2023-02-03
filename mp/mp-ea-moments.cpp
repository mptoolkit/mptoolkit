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
#include "mp-algorithms/triangular_mpo_solver_helpers.h"
#include "mp-algorithms/transfer.h"

namespace prog_opt = boost::program_options;

void PrintFormat(std::complex<double> const& Value, bool ShowRealPart, bool ShowImagPart,
                 bool ShowMagnitude, bool ShowArgument, bool ShowRadians)
{
   if (ShowRealPart)
   {
      std::cout << std::setw(20) << Value.real() << "    ";
   }
   if (ShowImagPart)
   {
      std::cout << std::setw(20) << Value.imag() << "    ";
   }
   if (ShowMagnitude)
   {
      std::cout << std::setw(20) << LinearAlgebra::norm_frob(Value) << "    ";
   }
   if (ShowArgument)
   {
      double Arg = std::arg(Value);
      if (!ShowRadians)
         Arg *= 180.0 / math_const::pi;
      std::cout << std::setw(20) << Arg << "    ";
   }
   std::cout << std::endl;
}

void ShowHeading(bool ShowRealPart, bool ShowImagPart,
                 bool ShowMagnitude, bool ShowArgument, bool ShowRadians)
{
   if (ShowRealPart)
      std::cout << "#real                   ";
   if (ShowImagPart)
      std::cout << "#imag                   ";
   if (ShowMagnitude)
      std::cout << "#magnitude              ";
   if (ShowArgument)
      std::cout << "#argument" << (ShowRadians ? "(rad)" : "(deg)") << "          ";
   std::cout << std::endl;
}

// Get the off-diagonal elements of the triangular unit cell for the excitation
// specified by the vector of (single-site) windows.
LinearWavefunction
ConstructPsiTri(LinearWavefunction PsiLeft, LinearWavefunction PsiRight,
                std::vector<WavefunctionSectionLeft> WindowVec)
{
   // Extract the (single-site) windows.
   std::vector<StateComponent> BVec;
   for (WavefunctionSectionLeft Window : WindowVec)
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
      auto CL = PsiLeft.begin();
      auto CR = PsiRight.begin();
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

   return PsiTri;
}

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      double K = 0.0;
      int Power = 1;
      bool CalculateMoments = false;
      bool CalculateMomentsFull = false;
      bool CalculateCumulants = false;
      bool ShowRealPart = false;
      bool ShowImagPart = false;
      bool ShowMagnitude = false;
      bool ShowArgument = false;
      bool ShowRadians = false;
      bool ShowPolar = false;
      bool ShowCartesian = false;
      int UnitCellSize = 0;
      int Degree = 0;
      bool Quiet = false;
      double Tol = 1e-15;
      double UnityEpsilon = DefaultEigenUnityEpsilon;
      bool Right = false;
      bool ShowAll = false;
      bool Product = false;
      std::string InputFilename;
      std::string InputFilename2;
      std::string OpStr;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("momentum,k", prog_opt::value(&K),
          "Use this momentum for the EA wavefunction instead of the one in the file (in units of pi)")
         ("power", prog_opt::value(&Power),
          FormatDefault("Calculate expectation value of operator to this power", Power).c_str())
         ("moments", prog_opt::bool_switch(&CalculateMoments),
          "Calculate the moments [default, unless --cumulants is specified]")
         ("moments-full", prog_opt::bool_switch(&CalculateMomentsFull),
          "Print the full moment polynomials containing the ground state and excitation components")
         ("cumulants,t", prog_opt::bool_switch(&CalculateCumulants),
          "Calculate the cumulants")
         ("cart,c", prog_opt::bool_switch(&ShowCartesian),
          "Show the result in cartesian coordinates [equivalent to --real --imag]")
         ("polar,p", prog_opt::bool_switch(&ShowPolar),
          "Show the result in polar coodinates [equivalent to --mag --arg]")
         ("real,r", prog_opt::bool_switch(&ShowRealPart),
          "Display the real part of the result")
         ("imag,i", prog_opt::bool_switch(&ShowImagPart),
          "Display the imaginary part of the result")
         ("mag,m", prog_opt::bool_switch(&ShowMagnitude),
          "Display the magnitude of the result")
         ("arg,a", prog_opt::bool_switch(&ShowArgument),
          "Display the argument (angle) of the result")
         ("radians", prog_opt::bool_switch(&ShowRadians),
          "Display the argument in radians instead of degrees")
         ("unitcell,u", prog_opt::value(&UnitCellSize),
          "Scale the results to use this unit cell size [default wavefunction unit cell]")
         ("degree,d", prog_opt::value(&Degree),
          "Force setting the degree of the MPO")
         ("quiet,q", prog_opt::bool_switch(&Quiet), "Don't show column headings")
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Linear solver convergence tolerance", Tol).c_str())
         ("unityepsilon", prog_opt::value(&UnityEpsilon),
          FormatDefault("Epsilon value for testing eigenvalues near unity", UnityEpsilon).c_str())
         ("right", prog_opt::bool_switch(&Right), "Calculate the moments in the opposite direction")
         ("showall", prog_opt::bool_switch(&ShowAll), "Show all columns of the fixed-point solutions (mostly for debugging)")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value<std::string>(&InputFilename), "wavefunction")
         ("operator", prog_opt::value<std::string>(&OpStr), "operator")
         ("psi2", prog_opt::value<std::string>(&InputFilename2), "wavefunction2")
         ("product", prog_opt::bool_switch(&Product), "Read a product operator instead of a triangular operator")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);
      p.add("operator", 1);
      p.add("psi2", 1);

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

      // If no output switches are used, default to showing everything
      if (!ShowRealPart && !ShowImagPart && !ShowMagnitude
          && !ShowCartesian && !ShowPolar && !ShowArgument)
      {
         ShowCartesian = true;
         ShowPolar = true;
      }

      if (ShowCartesian)
      {
         ShowRealPart = true;
         ShowImagPart = true;
      }
      if (ShowPolar)
      {
         ShowMagnitude = true;
         ShowArgument = true;
      }
      if (ShowRadians)
         ShowArgument = true;

      // If neither of --cumulants or --moments are specified, then the default is to calculate moments
      if (!CalculateMoments && !CalculateCumulants)
         CalculateMoments = true;

      if (CalculateMomentsFull)
         CalculateMoments = true;

      if (ShowAll)
      {
         CalculateMoments = true;
         CalculateMomentsFull = true;
      }

      // Load the wavefunction.
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(InputFilename, CacheSize, true);
      EAWavefunction Psi = PsiPtr->get<EAWavefunction>();

      EAWavefunction Psi2;
      if (vm.count("psi2"))
      {
         pvalue_ptr<MPWavefunction> Psi2Ptr = pheap::ImportHeap(InputFilename2);
         Psi2 = Psi2Ptr->get<EAWavefunction>();
      }
      else
         Psi2 = Psi;

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

      // Get the transfer matrix eigenvectors for the left and right boundaries.
      MatrixOperator RhoLeft = delta_shift(LambdaLeft*LambdaLeft, QShift);
      MatrixOperator RhoRight = LambdaRight*LambdaRight;

      MatrixOperator IdentLeft = MatrixOperator::make_identity(PsiLinearLeft.Basis1());
      MatrixOperator IdentRight = MatrixOperator::make_identity(PsiLinearRight.Basis1());

      if (Right)
      {
         RhoLeft.delta_shift(adjoint(QShift));
         RhoRight.delta_shift(adjoint(QShift));
         IdentLeft.delta_shift(adjoint(QShift));
         IdentRight.delta_shift(adjoint(QShift));
      }

      // Get the mixed transfer matrix eigenvectors.
      MatrixOperator TTopLeft, TTopRight, TBotLeft, TBotRight;

      if (!Right)
      {
         std::tie(std::ignore, TTopLeft, TTopRight) = get_transfer_unit_eigenpair(PsiLinearRight, PsiLinearLeft, QShift, Tol, UnityEpsilon, Verbose);
         TTopRight = delta_shift(TTopRight, QShift);
         std::tie(std::ignore, TBotLeft, TBotRight) = get_transfer_unit_eigenpair(PsiLinearLeft, PsiLinearRight, QShift, Tol, UnityEpsilon, Verbose);
         TBotRight = delta_shift(TBotRight, QShift);
      }
      else
      {
         std::tie(std::ignore, TTopLeft, TTopRight) = get_transfer_unit_eigenpair(PsiLinearLeft, PsiLinearRight, QShift, Tol, UnityEpsilon, Verbose);
         TTopLeft = delta_shift(TTopLeft, adjoint(QShift));
         std::tie(std::ignore, TBotLeft, TBotRight) = get_transfer_unit_eigenpair(PsiLinearRight, PsiLinearLeft, QShift, Tol, UnityEpsilon, Verbose);
         TBotLeft = delta_shift(TBotLeft, adjoint(QShift));
      }

      LinearWavefunction PsiTri = ConstructPsiTri(PsiLinearLeft, PsiLinearRight, Psi.window_vec());
      LinearWavefunction PsiTri2 = ConstructPsiTri(PsiLinearLeft, PsiLinearRight, Psi2.window_vec());

      // If the operator is not specified, use the one specified in the file attributes.
      if (!vm.count("operator"))
      {
         OpStr = PsiPtr->Attributes()["Hamiltonian"].get_or_default(std::string());
         if (OpStr.empty())
         {
            std::cerr << basename(argv[0]) << ": fatal: no operator specified, and wavefunction "
               "attribute Hamiltonian does not exist or is empty." << std::endl;
            return 1;
         }
      }

      // Load the operator.
      BasicTriangularMPO Op;
      InfiniteLattice Lattice;
      if (!Product)
         std::tie(Op, Lattice) = ParseTriangularOperatorAndLattice(OpStr);
      else
      {
         ProductMPO ProductOp;
         std::tie(ProductOp, Lattice) = ParseProductOperatorAndLattice(OpStr);
         Op = BasicTriangularMPO(std::vector<OperatorComponent>(ProductOp.begin(), ProductOp.end()));
      }

      // Ensure Op is the correct size.
      if (Op.size() < PsiLeft.size())
         Op = repeat(Op, PsiLeft.size() / Op.size());
      CHECK_EQUAL(Op.size(), PsiLeft.size());

      BasicTriangularMPO OriginalOp = Op;

      int LatticeUCsPerPsiUC = PsiLeft.size() / Lattice.GetUnitCell().size();

      // The default UnitCellSize for output is the wavefunction size
      if (UnitCellSize == 0)
         UnitCellSize = PsiLeft.size();
      double ScaleFactor = double(UnitCellSize) / double(PsiLeft.size());

      // Use the phase factor from the input file by default, otherwise, use the specified value.
      std::complex<double> ExpIK;
      if (!vm.count("momentum"))
         ExpIK = Psi.exp_ik();
      else
         ExpIK = exp(std::complex<double>(0.0, math_const::pi) * (K * LatticeUCsPerPsiUC));

      if (!Quiet)
      {
         print_preamble(std::cout, argc, argv);
         if (!vm.count("operator"))
         {
            std::cout << "#operator " << EscapeArgument(OpStr) << '\n';
         }
         std::cout << "#quantities are calculated per unit cell size of " << UnitCellSize
                   << (UnitCellSize == 1 ? " site\n" : " sites\n");

         // Print column headings.
         if (CalculateMoments)
            std::cout << "#moment ";
         if (ShowAll)
         {
            if (!Right)
               std::cout << "#matrix #column #momentum(rad/pi)    ";
            else
               std::cout << "#matrix #row    #momentum(rad/pi)    ";
         }
         if (CalculateMomentsFull)
            std::cout << "#degree ";
         if (ShowAll)
            std::cout << "#norm                   ";
         if (CalculateCumulants & !CalculateMoments)
            std::cout << "#cumulant ";
         ShowHeading(ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
         std::cout << std::left;
      }

      std::vector<std::complex<double>> Moments, Cumulants;

      // Loop over the powers of the operator.
      for (int p = 1; p <= Power; ++p)
      {
         Polynomial<std::complex<double>> FullMoment;
         std::vector<std::tuple<std::string, std::vector<KMatrixPolyType>, MatrixOperator>> FullData;

         if (!Right)
         {
            // Calculate the full E-matrix.
            std::vector<KMatrixPolyType> EMatK0;
            EMatK0.push_back(KMatrixPolyType());
            EMatK0[0][1.0] = MatrixPolyType(IdentLeft);
            SolveMPO_EA_Left(EMatK0, std::vector<KMatrixPolyType>(), PsiLinearLeft, PsiLinearLeft,
                             QShift, Op, IdentLeft, RhoLeft, 1.0,
                             Degree*p, Tol, UnityEpsilon, true, false, Verbose);

            std::vector<KMatrixPolyType> EMatKTop;
            std::vector<KMatrixPolyType> CTriKTop = CalculateCTriK_Left(std::vector<KMatrixPolyType>(), EMatK0, std::vector<KMatrixPolyType>(),
                                                                        PsiLinearRight, PsiLinearLeft, PsiTri, PsiTri2, QShift, Op, 1.0, 1.0);
            SolveMPO_EA_Left(EMatKTop, CTriKTop, PsiLinearRight, PsiLinearLeft,
                             QShift, Op, TTopLeft, TTopRight, ExpIK,
                             Degree*p, Tol, UnityEpsilon, true, false, Verbose);

            std::vector<KMatrixPolyType> EMatKBot;
            std::vector<KMatrixPolyType> CTriKBot = CalculateCTriK_Left(EMatK0, std::vector<KMatrixPolyType>(), std::vector<KMatrixPolyType>(),
                                                                        PsiLinearLeft, PsiLinearRight, PsiTri, PsiTri2, QShift, Op, 1.0, 1.0);
            SolveMPO_EA_Left(EMatKBot, CTriKBot, PsiLinearLeft, PsiLinearRight,
                             QShift, Op, TBotLeft, TBotRight, std::conj(ExpIK),
                             Degree*p, Tol, UnityEpsilon, true, false, Verbose);

            std::vector<KMatrixPolyType> EMatK1;
            std::vector<KMatrixPolyType> CTriK1 = CalculateCTriK_Left(EMatKTop, EMatKBot, EMatK0, PsiLinearRight, PsiLinearRight,
                                                                      PsiTri, PsiTri2, QShift, Op, ExpIK, ExpIK);
            SolveMPO_EA_Left(EMatK1, CTriK1, PsiLinearRight, PsiLinearRight,
                             QShift, Op, RhoRight, IdentRight, 1.0,
                             Degree*p, Tol, UnityEpsilon, false, false, Verbose);

            // Get the moment for the excitation.
            Moments.push_back(inner_prod(EMatK1.back()[1.0].coefficient(1), IdentRight));

            // Get the full moment if we need it.
            if (CalculateMomentsFull && !ShowAll)
               FullMoment = ExtractOverlap(EMatK1.back()[1.0], IdentRight);

            if (ShowAll)
            {
               FullData.push_back(std::tie("Initial ", EMatK0, RhoLeft));
               FullData.push_back(std::tie("Top     ", EMatKTop, TTopRight));
               FullData.push_back(std::tie("Bottom  ", EMatKBot, TBotRight));
               FullData.push_back(std::tie("Final   ", EMatK1, IdentRight));
            }
         }
         else
         {
            // Calculate the full F-matrix.
            std::vector<KMatrixPolyType> FMatK0;
            FMatK0.push_back(KMatrixPolyType());
            FMatK0[0][1.0] = MatrixPolyType(IdentRight);
            SolveMPO_EA_Right(FMatK0, std::vector<KMatrixPolyType>(), PsiLinearRight, PsiLinearRight,
                              QShift, Op, RhoRight, IdentRight, 1.0,
                              Degree*p, Tol, UnityEpsilon, true, false, Verbose);

            std::vector<KMatrixPolyType> FMatKTop;
            std::vector<KMatrixPolyType> CTriKTop = CalculateCTriK_Right(std::vector<KMatrixPolyType>(), FMatK0, std::vector<KMatrixPolyType>(),
                                                                         PsiLinearLeft, PsiLinearRight, PsiTri, PsiTri2, QShift, Op, 1.0, 1.0);
            SolveMPO_EA_Right(FMatKTop, CTriKTop, PsiLinearLeft, PsiLinearRight,
                              QShift, Op, TTopLeft, TTopRight, ExpIK,
                              Degree*p, Tol, UnityEpsilon, true, false, Verbose);

            std::vector<KMatrixPolyType> FMatKBot;
            std::vector<KMatrixPolyType> CTriKBot = CalculateCTriK_Right(FMatK0, std::vector<KMatrixPolyType>(), std::vector<KMatrixPolyType>(),
                                                                         PsiLinearRight, PsiLinearLeft, PsiTri, PsiTri2, QShift, Op, 1.0, 1.0);
            SolveMPO_EA_Right(FMatKBot, CTriKBot, PsiLinearRight, PsiLinearLeft,
                              QShift, Op, TBotLeft, TBotRight, std::conj(ExpIK),
                              Degree*p, Tol, UnityEpsilon, true, false, Verbose);

            std::vector<KMatrixPolyType> FMatK1;
            std::vector<KMatrixPolyType> CTriK1 = CalculateCTriK_Right(FMatKTop, FMatKBot, FMatK0, PsiLinearLeft, PsiLinearLeft,
                                                                       PsiTri, PsiTri2, QShift, Op, ExpIK, ExpIK);
            SolveMPO_EA_Right(FMatK1, CTriK1, PsiLinearLeft, PsiLinearLeft,
                              QShift, Op, IdentLeft, RhoLeft, 1.0,
                              Degree*p, Tol, UnityEpsilon, false, false, Verbose);

            // Get the moment for the excitation.
            Moments.push_back(inner_prod(IdentLeft, FMatK1.front()[1.0].coefficient(1)));

            if (ShowAll)
            {
               FullData.push_back(std::tie("Initial ", FMatK0, RhoRight));
               FullData.push_back(std::tie("Top     ", FMatKTop, TTopLeft));
               FullData.push_back(std::tie("Bottom  ", FMatKBot, TBotLeft));
               FullData.push_back(std::tie("Final   ", FMatK1, IdentLeft));
            }
            else if (CalculateMomentsFull)
               FullMoment = ExtractOverlap(FMatK1.front()[1.0], IdentLeft);
         }

         // Get the cumulant if desired.
         if (CalculateCumulants)
         {
            Cumulants.push_back(Moments.back());
            int n = Moments.size()-1; // Number of moments we already have.
            for (int i = 0; i < n; ++i)
               Cumulants.back() -= double(Binomial(n, i)) * Cumulants[i] * Moments[n-i-1];

            // Print the cumulant if we are not printing moments as well.
            if (!CalculateMoments)
            {
               std::cout << std::setw(9) << Cumulants.size() << " ";
               PrintFormat(Cumulants.back(), ShowRealPart, ShowImagPart, ShowMagnitude,
                           ShowArgument, ShowRadians);
            }
         }

         if (ShowAll)
         {
            // Print all of the columns of each E-matrix for each degree and momentum.
            for (auto const& I : FullData)
            {
               int Dim = std::get<1>(I).size();
               for (int i = 0; i < Dim; ++i)
               {
                  int Index = Right ? Dim-i-1 : i;
                  for (auto const& J : std::get<1>(I)[Index])
                  {
                     for (auto const& E : J.second)
                     {
                        if (E.second.TransformsAs() != std::get<2>(I).TransformsAs())
                           break;

                        std::cout << std::setw(7) << p << " "
                                  << std::get<0>(I)
                                  << std::setw(7) << Index << " "
                                  << std::setw(20) << std::arg(J.first)/math_const::pi << " "
                                  << std::setw(7) << E.first << " "
                                  << std::setw(20) << norm_frob(E.second) << "    ";
                        std::complex<double> x = inner_prod(std::get<2>(I), E.second)
                                               * std::pow(ScaleFactor, double(E.first-1));
                        if (Right)
                           x = std::conj(x);
                        PrintFormat(x, ShowRealPart, ShowImagPart, ShowMagnitude,
                                    ShowArgument, ShowRadians);
                     }
                  }
               }
            }
         }
         else if (CalculateMomentsFull)
         {
            // Print the full moment polynomials.
            for (auto const& I : FullMoment)
            {
               std::cout << std::setw(7) << p << " "
                         << std::setw(7) << I.first-1 << " ";
               std::complex<double> x = I.second * std::pow(ScaleFactor, double(I.first-1));
               if (Right)
                  x = std::conj(x);
               PrintFormat(x, ShowRealPart, ShowImagPart, ShowMagnitude,
                           ShowArgument, ShowRadians);
            }
         }
         else if (CalculateMoments)
         {
            // Print the moment.
            std::cout << std::setw(7) << Moments.size() << " ";
            PrintFormat(Moments.back(), ShowRealPart, ShowImagPart, ShowMagnitude,
                        ShowArgument, ShowRadians);
         }

         // Get the next power of the MPO.
         if (p < Power)
            Op = Op * OriginalOp;
      }

      // If we wanted both moments and cumulants, print the cumulants now.
      if (CalculateMoments && CalculateCumulants)
      {
         std::cout << std::endl;
         if (!Quiet)
         {
            std::cout << "#cumulant ";
            ShowHeading(ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
         }
         for (int i = 0; i < Cumulants.size(); ++i)
         {
            std::cout << std::setw(9) << i << " ";
            PrintFormat(Cumulants[i], ShowRealPart, ShowImagPart, ShowMagnitude,
                        ShowArgument, ShowRadians);
         }
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
