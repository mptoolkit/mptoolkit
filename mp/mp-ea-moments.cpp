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
#include "common/terminal.h"
#include "common/prog_options.h"
#include "interface/inittemp.h"
#include "lattice/infinite-parser.h"
#include "mp-algorithms/ef-matrix.h"
#include "mp-algorithms/triangular_mpo_solver_helpers.h"

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

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      double K = 0.0;
      int LatticeUCSize = 1;
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
      bool Right = false;
      bool ShowAll = false;
      bool String = false;
      bool Overlap = false;
      std::string InputFilename;
      std::string InputFilename2;
      std::string OpStr;

      EFMatrixSettings Settings;
      Settings.Tol = 1e-15;
      Settings.NeedFinalMatrix = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("operator,o", prog_opt::value(&OpStr), "Calculate the expectation value of this triangular operator")
         ("string", prog_opt::value(&OpStr), "Calculate the expectation value of this string operator")
         ("overlap", prog_opt::bool_switch(&Overlap), "Calculate the overlap")
         ("momentum,k", prog_opt::value(&K),
          "Use this momentum for the EA wavefunction instead of the one in the file (in units of pi)")
         ("latticeucsize", prog_opt::value(&LatticeUCSize), "Lattice unit cell size [default wavefunction attribute \"LatticeUnitCellSize\" or 1]")
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
         ("tol", prog_opt::value(&Settings.Tol),
          FormatDefault("Linear solver convergence tolerance", Settings.Tol).c_str())
         ("unityepsilon", prog_opt::value(&Settings.UnityEpsilon),
          FormatDefault("Epsilon value for testing eigenvalues near unity", Settings.UnityEpsilon).c_str())
         ("right", prog_opt::bool_switch(&Right), "Calculate the moments in the opposite direction")
         ("showall", prog_opt::bool_switch(&ShowAll), "Show all columns of the fixed-point solutions (mostly for debugging)")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value(&InputFilename), "wavefunction")
         ("psi2", prog_opt::value(&InputFilename2), "wavefunction2")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);
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
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi> [psi2]" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      Settings.Verbose = Verbose;
      Settings.Degree = Degree;

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

      if (vm.count("operator") && vm.count("string"))
      {
         std::cerr << "fatal: cannot use --operator and --string simultaneously." << std::endl;
         return 1;
      }

      if (vm.count("string") || Overlap)
      {
         String = true;

         if (Power > 1)
         {
            std::cerr << "fatal: cannot use --power with string operators." << std::endl;
            return 1;
         }
      }

      // Load the wavefunctions.
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(InputFilename, CacheSize, true);
      EAWavefunction Psi = PsiPtr->get<EAWavefunction>();

      // Get the lattice unit cell size if unspecified.
      if (!vm.count("latticeucsize"))
         LatticeUCSize = PsiPtr->Attributes()["LatticeUnitCellSize"].get_or_default<int>(1);

      EAWavefunction Psi2;
      if (vm.count("psi2"))
      {
         pvalue_ptr<MPWavefunction> Psi2Ptr = pheap::ImportHeap(InputFilename2);
         Psi2 = Psi2Ptr->get<EAWavefunction>();
      }
      else
         Psi2 = Psi;

      // If the operator is not specified, use the one specified in the file attributes.
      if (!vm.count("operator") && !String)
      {
         OpStr = PsiPtr->Attributes()["Hamiltonian"].get_or_default(std::string());
         if (OpStr.empty())
         {
            std::cerr << "fatal: no operator specified, and wavefunction attribute Hamiltonian does not exist." << std::endl;
            return 1;
         }
      }

      // Load the operator.
      InfiniteLattice Lattice;
      InfiniteMPO Op;

      if (!String)
      {
         BasicTriangularMPO TriOp;
         std::tie(TriOp, Lattice) = ParseTriangularOperatorAndLattice(OpStr);

         // Ensure TriOp is the correct size.
         if (TriOp.size() < Psi.left().size())
            TriOp = repeat(TriOp, Psi.left().size() / TriOp.size());
         CHECK_EQUAL(TriOp.size(), Psi.left().size());

         Op = TriOp;
      }
      else
      {
         ProductMPO StringOp = ProductMPO::make_identity(ExtractLocalBasis(Psi.left()));

         if (!Overlap)
            std::tie(StringOp, Lattice) = ParseProductOperatorAndLattice(OpStr);

         // Ensure StringOp is the correct size.
         if (StringOp.size() < Psi.left().size())
            StringOp = repeat(StringOp, Psi.left().size() / StringOp.size());
         CHECK_EQUAL(StringOp.size(), Psi.left().size());

         Op = StringOp;
      }

      // Extract the left and right semi-infinite boundaries.
      InfiniteWavefunctionLeft PsiLeft = Psi.left();
      inplace_qshift(PsiLeft, Psi.left_qshift());
      PsiLeft.rotate_left(Psi.left_index());

      InfiniteWavefunctionRight PsiRight = Psi.right();
      inplace_qshift(PsiRight, Psi.right_qshift());
      PsiRight.rotate_left(Psi.right_index());

      // Extract the windows.
      int WindowSize = Psi.window_size();
      int WindowSize2 = Psi2.window_size();
      std::vector<LinearWavefunction> WindowVec, WindowVec2;

      for (WavefunctionSectionLeft Window : Psi.window_vec())
      {
         LinearWavefunction PsiLinear;
         MatrixOperator U;
         std::tie(PsiLinear, U) = get_left_canonical(Window);
         PsiLinear.set_back(PsiLinear.get_back()*U);
         WindowVec.push_back(PsiLinear);
      }

      for (WavefunctionSectionLeft Window : Psi2.window_vec())
      {
         LinearWavefunction PsiLinear;
         MatrixOperator U;
         std::tie(PsiLinear, U) = get_left_canonical(Window);
         PsiLinear.set_back(PsiLinear.get_back()*U);
         WindowVec2.push_back(PsiLinear);
      }

      InfiniteMPO OriginalOp = Op;

      // The default UnitCellSize for output is the wavefunction size
      if (UnitCellSize == 0)
         UnitCellSize = PsiLeft.size();
      double ScaleFactor = double(UnitCellSize) / double(PsiLeft.size());

      // If we are calculating the overlap, then we do not have the lattice,
      // but ExpIK shouldn't matter in that case anyway, so we can set this to zero.
      int LatticeUCsPerPsiUC = PsiLeft.size() / LatticeUCSize;

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
            std::cout << "#operator " << EscapeArgument(OpStr) << std::endl;
         }
         std::cout << "#quantities are calculated per unit cell size of " << UnitCellSize
                   << (UnitCellSize == 1 ? " site" : " sites") << std::endl;

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

      EFMatrix EF(Op, Settings);
      EF.SetPsi({0}, PsiLeft);
      EF.SetPsi({Infinity}, PsiRight, ExpIK);
      EF.SetWindowUpper(WindowVec);
      EF.SetWindowLower(WindowVec2);

      // Loop over the powers of the operator.
      for (int p = 1; p <= Power; ++p)
      {
         if (p > 1)
         {
            // Set Op to the next power.
            Op = Op * OriginalOp;
            EF.SetOp(Op, Degree * p);
         }

         // The index for the final element
         ZeroInf FI(!Right);

         // Get the moment for the excitation.
         MatrixPolyType FinalElement = !Right ? EF.GetElement({FI}, {FI}, Right).back()[1.0] : EF.GetElement({FI}, {FI}, Right).front()[1.0];

         Moments.push_back(inner_prod(EF.GetRho({FI}, {FI}, Right), FinalElement.coefficient(1)));
         if (Right) // Conjugate for the right.
            Moments.back() = std::conj(Moments.back());

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
            for (int i = 0; i <= 1; ++i)
            {
               for (int j = 0; j <= 1; ++j)
               {
                  ZeroInf I(!Right != (i == 0));
                  ZeroInf J(!Right != (j == 0));

                  std::vector<KMatrixPolyType> Element = EF.GetElement({I}, {J}, Right);
                  int Dim = Element.size();
                  for (int n = 0; n < Dim; ++n)
                  {
                     int Index = Right ? Dim-n-1 : n;
                     for (auto const& K : Element[Index])
                     {
                        for (auto const& E : K.second)
                        {
                           if (EF.GetRho({I}, {J}, Right).is_null())
                              break;

                           if (E.second.TransformsAs() != EF.GetRho({I}, {J}, Right).TransformsAs())
                              break;

                           std::cout << std::setw(7) << p << " "
                                     << std::setw(3) << int(I.is_inf()) << " "
                                     << std::setw(3) << int(J.is_inf()) << " "
                                     << std::setw(7) << Index << " "
                                     << std::setw(20) << std::arg(K.first)/math_const::pi << " "
                                     << std::setw(7) << E.first << " "
                                     << std::setw(20) << norm_frob(E.second) << "    ";
                           std::complex<double> x = inner_prod(EF.GetRho({I}, {J}, Right), E.second)
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
         }
         else if (CalculateMomentsFull)
         {
            // Print the full moment polynomials.
            for (auto const& I : ExtractOverlap(FinalElement, EF.GetRho({FI}, {FI}, Right)))
            {
               std::cout << std::setw(7) << p << " "
                         << std::setw(7) << I.first << " ";
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
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << std::endl;
      pheap::Cleanup();
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
