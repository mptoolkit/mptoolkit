// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-ispectrum.cpp
//
// Copyright (C) 2004-2022 Ian McCulloch <ian@qusim.net>
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

#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"
#include "mp-algorithms/transfer.h"
#include "wavefunction/operator_actions.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "lattice/siteoperator-parser.h"
#include "linearalgebra/arpack_wrapper.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcell-parser.h"
#include "lattice/infinite-parser.h"

namespace prog_opt = boost::program_options;


void PrintFormat(QuantumNumber const& q, std::complex<double> x, int n, bool ShowRealPart, bool ShowImagPart,
                 bool ShowCorrLength, bool ShowMagnitude, bool ShowArgument,
                 bool ShowRadians, double ScaleFactor)
{
   std::string SectorStr = boost::lexical_cast<std::string>(q);
   std::complex<double> Value = std::pow(x, ScaleFactor);
   std::cout << std::setw(11) << SectorStr << ' ';
   std::cout << std::setw(4) << n << ' ';
   if (ShowRealPart)
   {
      std::cout << std::setw(20) << Value.real() << "  ";
   }
   if (ShowImagPart)
   {
      std::cout << std::setw(20) << Value.imag() << "  ";
   }
   if (ShowCorrLength)
   {
      std::cout << std::setw(20) << (-1.0/std::log(LinearAlgebra::norm_frob(Value)))
                << "  ";
   }
   if (ShowMagnitude)
   {
      std::cout << std::setw(20) << LinearAlgebra::norm_frob(Value) << "  ";
   }
   if (ShowArgument)
   {
      double Arg =  std::atan2(Value.imag(), Value.real());
      if (!ShowRadians)
         Arg *= 180.0 / math_const::pi;
      std::cout << std::setw(20) << Arg << "  ";
   }
}

// inject_left for a BasicFiniteMPO.  This can have support on multiple wavefunction unit cells
int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      bool ShowRealPart = false, ShowImagPart = false, ShowMagnitude = false;
      bool ShowCartesian = false, ShowPolar = false, ShowArgument = false;
      bool ShowRadians = false, ShowCorrLength = false;
      int UnitCellSize = 0;
      std::string PsiStr;
      std::vector<std::string> Sector;
      double Tol = 1E-15;
      bool Sort = false;
      bool Quiet = false;
      bool Print = false;
      bool Symmetric = false;
      int KrylovLength = 0;
      std::string String;
      int MaxEigen = 10;
      std::vector<std::string> LeftOpStr;
      std::vector<std::string> RightOpStr;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("cart,c", prog_opt::bool_switch(&ShowCartesian),
          "show the result in cartesian coordinates [equivalent to --real --imag]")
         ("polar,p", prog_opt::bool_switch(&ShowPolar),
          "show the result in polar coodinates [equivalent to --mag --arg]")
         ("real,r", prog_opt::bool_switch(&ShowRealPart),
          "display the real part of the result")
         ("imag,i", prog_opt::bool_switch(&ShowImagPart),
          "display the imaginary part of the result")
         ("mag,m", prog_opt::bool_switch(&ShowMagnitude),
          "display the magnitude of the result")
         ("arg,a", prog_opt::bool_switch(&ShowArgument),
          "display the argument (angle) of the result")
         ("radians", prog_opt::bool_switch(&ShowRadians),
          "display the argument in radians instead of degrees")
         ("corr,x", prog_opt::bool_switch(&ShowCorrLength),
          "display the equivalent correlation length")
         ("unitcell,u", prog_opt::value(&UnitCellSize),
          "scale the results to use this unit cell size [default wavefunction unit cell]")
         ("numeigen,n", prog_opt::value(&MaxEigen),
          FormatDefault("Number of eigenvalues to calculate in each sector", MaxEigen).c_str())
         ("left", prog_opt::value(&LeftOpStr),
          "Calculate the expansion coefficients of this operator acting on the left")
         ("right", prog_opt::value(&RightOpStr),
          "Calculate the expansion coefficients of this operator acting on the right")
         ("string", prog_opt::value(&String),
          "use this product operator as a string operator")
         ("quantumnumber,q", prog_opt::value(&Sector),
          "calculate the overlap only in this quantum number sector, "
          "can be used multiple times [default is to calculate all sectors]")
         ("symmetric", prog_opt::bool_switch(&Symmetric),
          "transform into the symmetric canonical form")
         ("sort,s", prog_opt::bool_switch(&Sort),
          "sort the eigenvalues by magnitude (not yet implemented)")
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Tolerance of the Arnoldi eigensolver", Tol).c_str())
         //("iter", prog_opt::value(&Iter),
         //FormatDefault("Maximum subspace size in the Arnoldi basis", Iter).c_str())
         ("krylov,k", prog_opt::value(&KrylovLength),
          "Length of the Krylov sequence [default 2*num-eigenvalues]")
         ("quiet", prog_opt::bool_switch(&Quiet),
          "don't show the preamble and column headings")
         ("print", prog_opt::bool_switch(&Print), "with --string, Print the MPO to standard output")
         //         ("overlaps", prog_opt::bool_switch(&ShowEigenvectorOverlaps),
         //"Write the matrix of overlaps of the left/right eigenvectors to cerr (for debugging)")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value(&PsiStr), "psi1")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("psi") < 1)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (!Quiet)
         print_preamble(std::cout, argc, argv);

      // If no output switches are used, default to showing everything
      if (!ShowRealPart && !ShowImagPart && !ShowMagnitude
          && !ShowCartesian && !ShowPolar && !ShowArgument
          && !ShowCorrLength)
      {
         ShowCartesian = true;
         ShowPolar = true;
         ShowCorrLength = true;
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

      SimpleOperator MyOpL, MyOpR, MyOpL2;
      LatticeSite Site;

      if (Verbose)
         std::cout << "Loading wavefunction...\n";

      pvalue_ptr<MPWavefunction> Psi1
         = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);

      // firstly, get the LinearWavefunction
      InfiniteWavefunctionLeft InfPsi = Psi1->get<InfiniteWavefunctionLeft>();

      LinearWavefunction Psi;
      QuantumNumber QShift = InfPsi.qshift();

      RealDiagonalOperator D;
      std::tie(Psi, D) = get_left_canonical(InfPsi);

      if (Symmetric)
      {
         PANIC("not supported");
         // MatrixOperator LambdaSqrt = D;
         // LambdaSqrt = SqrtDiagonal(LambdaSqrt);
         // MatrixOperator LambdaInvSqrt = InvertDiagonal(LambdaSqrt, InverseTol);
         // Psi.set_back(prod(Psi.get_back(), LambdaSqrt));
         // Psi.set_front(prod(delta_shift(LambdaInvSqrt, QShift), Psi.get_front()));
      }


      // now get the principal eigenpair
      MatrixOperator LeftIdent, RightIdent;
      if (Symmetric)
      {
         if (Verbose)
            std::cout << "Solving principal eigenpair...\n";
         std::complex<double> EValue;
         std::tie(EValue, LeftIdent, RightIdent) = get_transfer_eigenpair(Psi, Psi, QShift, Tol, Verbose);
         RightIdent = delta_shift(RightIdent, QShift);
      }
      else
      {
         // if we're in the left canonical basis, then we know what the eigenpair is.
         LeftIdent = MatrixOperator::make_identity(Psi.Basis1());
         RightIdent = D*D;
      }

      std::complex<double> IdentNormalizationFactor = inner_prod(LeftIdent, delta_shift(RightIdent, QShift));

      // Get the string operator
      ProductMPO StringOp;
      if (vm.count("string"))
      {
         InfiniteLattice Lattice;
         std::tie(StringOp, Lattice) = ParseProductOperatorAndLattice(String);
         if (Print)
         {
            std::cout << "String MPO is:\n" << StringOp << '\n';
         }
         CHECK(Psi.size() % StringOp.size() == 0)
            ("Wavefunction size must be a multiple of the string operator size")
            (Psi.size())(StringOp.size());
         StringOp = repeat(StringOp, Psi.size() / StringOp.size());
      }
      else
      {
         StringOp = ProductMPO::make_identity(ExtractLocalBasis(Psi));
      }

      // The default UnitCellSize for output is the wavefunction size
      if (UnitCellSize == 0)
      {
         UnitCellSize = Psi.size();
      }
      double ScaleFactor = double(UnitCellSize) / double(Psi.size());

      // get the set of quantum numbers to show
      typedef std::set<QuantumNumbers::QuantumNumber> QSetType;
      QSetType QL;
      if (vm.count("quantumnumber") != 0)
      {
         while (!Sector.empty())
         {
            QL.insert(QuantumNumbers::QuantumNumber(Psi.GetSymmetryList(), Sector.back()));
            Sector.pop_back();
         }
      }
      else
      {
         // Assemble the list of all possible quantum numbers
         QSetType Q1 = QuantumNumbersInBasis(Psi.Basis2());
         for (QSetType::const_iterator I = Q1.begin(); I != Q1.end(); ++I)
         {
            for (QSetType::const_iterator J = Q1.begin(); J != Q1.end(); ++J)
            {
               QuantumNumberList NextQ = inverse_transform_targets(*I, *J);
               QL.insert(NextQ.begin(), NextQ.end());
            }
         }
      }

      // show the title
      if (!Quiet)
      {
         std::cout << "\n#quantities are calculated per unit cell size of " << UnitCellSize
                   << (UnitCellSize == 1 ? " site\n" : " sites\n");
         std::cout << "#sector     #n   ";
         if (ShowRealPart)
            std::cout << "#real                 ";
         if (ShowImagPart)
            std::cout << "#imag                 ";
         if (ShowCorrLength)
            std::cout << "#corr_length          ";
         if (ShowMagnitude)
            std::cout << "#magnitude            ";
         if (ShowArgument)
            std::cout << "#argument" << (ShowRadians ? "(rad)" : "(deg)") << "        ";

         // titles for the overlaps
         for (unsigned i = 0; i < LeftOpStr.size(); ++i)
         {
            for (unsigned j = 0; j < RightOpStr.size(); ++j)
            {
               std::string Title = "#overlap_" + boost::lexical_cast<std::string>(i) + '_'
                  + boost::lexical_cast<std::string>(j);
               if (ShowRealPart)
                  std::cout << std::setw(20) << std::left << (Title+"_real") << "  ";
               if (ShowImagPart)
                  std::cout << std::setw(20) << std::left << (Title+"_imag") << "  ";
               if (ShowMagnitude)
                  std::cout << std::setw(20) << std::left << (Title+"_mag") << "  ";
               if (ShowArgument)
                  std::cout << std::setw(20) << std::left << (Title+"_arg" + (ShowRadians ? "(rad)" : "(deg)")) << "  ";
            }
         }
         std::cout << '\n';
         std::cout << std::left;
      }

      // initialize the finite operators for the observables
      std::vector<MatrixOperator> LeftOp;
      std::vector<MatrixOperator> RightOp;

      //LeftIdent = delta_shift(LeftIdent, QShift);

      for (unsigned i = 0; i < LeftOpStr.size(); ++i)
      {
         // TODO: the MPO types here are a bit weird, converting to a finite MPO then back to a GenericMPO...
         UnitCellMPO Op = ParseUnitCellOperatorAndLattice(LeftOpStr[i]).first;
         Op.ExtendToCoverUnitCell(Psi.size());
         // Adjust the MPO so that the left basis is the identity
         BasicFiniteMPO Mpo = Op.MPO() * BasicFiniteMPO::make_identity(Op.MPO().LocalBasis2List(),
                                                              adjoint(Op.TransformsAs()));
         Mpo = project(Mpo, QuantumNumbers::QuantumNumber(Op.GetSymmetryList()));
         // Don't do a qshift here; this leaves LeftOp in the wavefunction Basis2()
         LeftOp.push_back(inject_left(LeftIdent, Psi, Mpo.data(), Psi));
      }

      for (unsigned i = 0; i < RightOpStr.size(); ++i)
      {
         UnitCellMPO Op = ParseUnitCellOperatorAndLattice(RightOpStr[i]).first;
         Op.ExtendToCoverUnitCell(Psi.size());
         // Don't do a qshift here; this leaves RightOp in the wavefunction Basis1()
         RightOp.push_back(inject_right(RightIdent, Psi, Op.MPO(), Psi));
      }

      // iterate over the relevant quantum number sectors
      for (QSetType::const_iterator qI = QL.begin(); qI != QL.end(); ++qI)
      {
         LinearAlgebra::Vector<MatrixOperator> LeftEigenvectors;
         LinearAlgebra::Vector<MatrixOperator> RightEigenvectors;
         LinearAlgebra::Vector<std::complex<double> > EValues;

         LinearAlgebra::Vector<MatrixOperator>* LeftEigenvectorsPtr = NULL;
         LinearAlgebra::Vector<MatrixOperator>* RightEigenvectorsPtr = NULL;

         if (!RightOp.empty())
         {
            RightEigenvectorsPtr = &RightEigenvectors;
            LeftEigenvectorsPtr = &LeftEigenvectors;
         }

         // determine the spectrum
         EValues = get_spectrum_string(Psi, QShift, StringOp * ProductMPO::make_identity(StringOp.LocalBasis2List(), *qI),
                                       MaxEigen, Tol, LeftEigenvectorsPtr,
                                       RightEigenvectorsPtr, KrylovLength, true, Verbose);

         for (int i = 0; i < int(size(EValues)); ++i)
         {
            PrintFormat(*qI, EValues[i], i, ShowRealPart, ShowImagPart, ShowCorrLength, ShowMagnitude,
                        ShowArgument, ShowRadians, ScaleFactor);

            // show eigenvector info?
            for (unsigned iL = 0; iL < LeftOp.size(); ++iL)
            {
               for (unsigned iR = 0; iR < RightOp.size(); ++iR)
               {
                  if (LeftOp[iL].TransformsAs() == LeftEigenvectors[i].TransformsAs()
                      && RightOp[iR].TransformsAs() == RightEigenvectors[i].TransformsAs())
                  {
                     std::complex<double> Overlap = inner_prod(LeftOp[iL], RightEigenvectors[i])
                        * inner_prod(LeftEigenvectors[i], RightOp[iR])
                        / (inner_prod(LeftEigenvectors[i], delta_shift(RightEigenvectors[i], QShift)) * IdentNormalizationFactor);

                     if (ShowRealPart)
                        std::cout << std::setw(20) << Overlap.real() << "  ";
                     if (ShowImagPart)
                        std::cout << std::setw(20) << Overlap.imag() << "  ";
                     double Magnitude = norm_frob(Overlap);
                     double Arg = std::atan2(Overlap.imag(), Overlap.real());
                     if (!ShowRadians)
                        Arg *= 180.0 / math_const::pi;
                     if (ShowMagnitude)
                        std::cout << std::setw(20) << Magnitude << "  ";
                     if (ShowArgument)
                        std::cout << std::setw(20) << Arg << "  ";
                  }
               }
            }
            std::cout << '\n';

         }// for i

      } // for qI

      pheap::Shutdown();
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
