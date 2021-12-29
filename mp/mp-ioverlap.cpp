// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ioverlap.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "lattice/unitcell.h"
#include "lattice/infinite-parser.h"
#include "common/statistics.h"

namespace prog_opt = boost::program_options;

struct TransEigenInfo
{
   TransEigenInfo() {}
   TransEigenInfo(QuantumNumber const& q_, std::complex<double> const& x_) : q(q_), n(1), x(x_) {}
   TransEigenInfo(QuantumNumber const& q_, int n_, std::complex<double> const& x_) : q(q_), n(n_), x(x_) {}

   QuantumNumber q;
   int n;
   std::complex<double> x;
};

// sort the magnitude of the eigenvalue in reverse order
bool operator<(TransEigenInfo const& x, TransEigenInfo const& y)
{
   return LinearAlgebra::norm_frob_sq(x.x) > LinearAlgebra::norm_frob_sq(y.x);
}

void PrintFormat(QuantumNumber const& q, int n, std::complex<double> x, int NumEigen,
                 bool ShowRealPart, bool ShowImagPart,
                 bool ShowCorrLength, bool ShowRate, bool ShowMagnitude, bool ShowArgument,
                 bool ShowRadians, double ScaleFactor)
{
   std::string SectorStr = boost::lexical_cast<std::string>(q);
   std::complex<double> Value = std::pow(x, ScaleFactor);
   std::cout << std::setw(11) << SectorStr << ' ';
   if (NumEigen > 1)
   {
      std::cout << std::setw(6) << n << ' ';
   }
   if (ShowRealPart)
   {
      std::cout << std::setw(20) << Value.real() << "    ";
   }
   if (ShowImagPart)
   {
      std::cout << std::setw(20) << Value.imag() << "    ";
   }
   if (ShowCorrLength)
   {
      std::cout << std::setw(20) << (-1.0/std::log(LinearAlgebra::norm_frob(Value)))
                << "    ";
   }
   if (ShowRate)
   {
      std::cout << std::setw(20) << (-std::log(LinearAlgebra::norm_frob(Value)))
                << "    ";
   }
   if (ShowMagnitude)
   {
      std::cout << std::setw(20) << LinearAlgebra::norm_frob(Value) << "    ";
   }
   if (ShowArgument)
   {
      double Arg =  std::atan2(Value.imag(), Value.real());
      if (!ShowRadians)
         Arg *= 180.0 / math_const::pi;
      std::cout << std::setw(20) << Arg << "    ";
   }
   std::cout << std::endl;
}

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      bool UseTempFile = false;
      bool ShowRealPart = false, ShowImagPart = false, ShowMagnitude = false;
      bool ShowCartesian = false, ShowPolar = false, ShowArgument = false;
      bool ShowRadians = false, ShowCorrLength = false, ShowRate = false;
      int Rotate = 0;
      int UnitCellSize = 0;
      std::string LhsStr, RhsStr;
      std::vector<std::string> Sector;
      double Tol = 1E-15;
      int Iter = 40;         // 2017-04-17: increased default from 30 to 40
      bool Sort = false;
      bool Quiet = false;
      bool Reflect = false;
      bool Conj = false;
      bool Print = false;
      int CoarseGrain1 = 1;
      int CoarseGrain2 = 1;
      int NumEigen = 1; // number of eigenvalues to calculate
      std::string String;

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
         ("rate", prog_opt::bool_switch(&ShowRate),
          "display the rate function -*log(real_part)")
         ("unitcell,u", prog_opt::value(&UnitCellSize),
          "scale the results to use this unit cell size [default wavefunction unit cell]")
         ("tempfile", prog_opt::bool_switch(&UseTempFile),
          "a temporary data file for workspace (path set by environment MP_BINPATH)")
         ("rotate", prog_opt::value(&Rotate),
          "rotate the unit cell of psi1 this many sites to the left before calculating the overlap [default 0]")
         ("reflect", prog_opt::bool_switch(&Reflect),
          "reflect psi1 (gives parity eigenvalue)")
         ("coarsegrain1", prog_opt::value(&CoarseGrain1),
          "coarse-grain wavefunction 1 by this amount")
         ("coarsegrain2", prog_opt::value(&CoarseGrain2),
          "coarse-grain wavefunction 2 by this amount")
         ("string", prog_opt::value(&String),
          "use this product operator as a string operator for the overlap")
         ("conj", prog_opt::bool_switch(&Conj),
          "complex conjugate psi1")
         ("q,quantumnumber", prog_opt::value(&Sector),
          "calculate the overlap only in this quantum number sector, "
          "can be used multiple times [default is to calculate all sectors]")
         ("numeigen,n", prog_opt::value(&NumEigen), "Calculate this many eigenvalues in each sector [default 1]")
         ("sort,s", prog_opt::bool_switch(&Sort),
          "sort the eigenvalues by magnitude")
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Tolerance of the Arnoldi eigensolver", Tol).c_str())
         ("iter", prog_opt::value(&Iter),
          FormatDefault("Maximum subspace size in the Arnoldi basis", Iter).c_str())
         ("quiet", prog_opt::bool_switch(&Quiet),
          "don't show the column headings")
         ("print", prog_opt::bool_switch(&Print), "with --string, Print the MPO to standard output")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("lhs", prog_opt::value<std::string>(&LhsStr), "psi1")
         ("rhs", prog_opt::value<std::string>(&RhsStr), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("lhs", 1);
      p.add("rhs", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("lhs") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi1> [<psi2>]\n";
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
          && !ShowRadians && !ShowCorrLength && !ShowRate)
      {
         ShowCartesian = true;
         ShowPolar = true;
         ShowCorrLength = true;
         ShowRate = true;
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

      if (Verbose)
         std::cout << "Loading LHS wavefunction...\n";

      pvalue_ptr<MPWavefunction> Psi1Ptr;
      if (UseTempFile)
      {
         mp_pheap::InitializeTempPHeap(Verbose);
         Psi1Ptr = pheap::ImportHeap(LhsStr);
      }
      else
      {
         Psi1Ptr = pheap::OpenPersistent(LhsStr, mp_pheap::CacheSize(), true);
      }
      InfiniteWavefunctionLeft Psi1 = Psi1Ptr->get<InfiniteWavefunctionLeft>();

      if (Verbose)
         std::cout << "Loading RHS wavefunction...\n";

      if (CoarseGrain1 != 1)
      {
         if (Psi1.size() % CoarseGrain1 != 0)
         {
            std::cerr << "Wavefunction 1 size is not a multiple of the coarsegran size, expanding...\n";
            Psi1 = repeat(Psi1, statistics::lcm(Psi1.size(), CoarseGrain1) / Psi1.size());
         }

         LinearWavefunction PsiL;
         RealDiagonalOperator Lambda;
         std::tie(PsiL, Lambda) = get_left_canonical(Psi1);

         Psi1 = InfiniteWavefunctionLeft::ConstructFromOrthogonal(coarse_grain(PsiL, CoarseGrain1), Lambda, Psi1.qshift());
      }

      int Size = Psi1.size();

      InfiniteWavefunctionLeft Psi2;
      if (RhsStr.empty())
         Psi2 = Psi1;
      else
      {
         pvalue_ptr<MPWavefunction> Psi2Ptr = pheap::ImportHeap(RhsStr);
         Psi2 = Psi2Ptr->get<InfiniteWavefunctionLeft>();

         if (CoarseGrain2 != 1)
         {
            if (Psi2.size() % CoarseGrain1 != 0)
            {
               std::cerr << "Wavefunction 2 size is not a multiple of the coarsegran size, expanding...\n";
               Psi2 = repeat(Psi1, statistics::lcm(Psi2.size(), CoarseGrain2) / Psi2.size());
            }

            LinearWavefunction PsiL;
            RealDiagonalOperator Lambda;
            std::tie(PsiL, Lambda) = get_left_canonical(Psi2);

            Psi2 = InfiniteWavefunctionLeft::ConstructFromOrthogonal(coarse_grain(PsiL, CoarseGrain2), Lambda, Psi2.qshift());
         }
      }

      if (Verbose)
         std::cout << "Calculating overlap...\n";

      // Rotate as necessary.  Do this BEFORE determining the quantum number sectors!
      if (Verbose)
      {
         std::cout << "Rotating Psi1 right by" << Rotate << " sites\n";
      }
      Psi1.rotate_right(Rotate);
      if (Reflect)
      {
         if (Verbose)
            std::cout << "Reflecting psi1..." << std::endl;
         inplace_reflect(Psi1);
      }
      if (Conj)
      {
         if (Verbose)
            std::cout << "Conjugating psi1..." << std::endl;
         inplace_conj(Psi1);
      }

      UnitCell Cell;
      LatticeSite Site;
      ProductMPO StringOp;
      if (vm.count("string"))
      {
         InfiniteLattice Lattice;
         std::tie(StringOp, Lattice) = ParseProductOperatorAndLattice(String);
         if (Print)
         {
            std::cout << "String MPO is:\n" << StringOp << '\n';
         }
      }
      else
      {
         StringOp = ProductMPO::make_identity(ExtractLocalBasis(Psi2));
      }

      Size = statistics::lcm(Psi1.size(), Psi2.size(), StringOp.size());
      StringOp = repeat(StringOp, Size / StringOp.size());
      if (Verbose > 0 && Psi1.size() != Psi2.size())
      {
         std::cerr << "Wavefunction unit cells differ, extending wavefunctions to size " << Size << '\n';
      }
      Psi1 = repeat(Psi1, Size/Psi1.size());
      Psi2 = repeat(Psi2, Size/Psi2.size());

      if (ExtractLocalBasis2(StringOp) != ExtractLocalBasis(Psi2))
      {
         std::cerr << "mp-ioverlap: fatal: local basis for RHS wavefunction does not match!\n";
         return 1;
      }
      if (ExtractLocalBasis1(StringOp) != ExtractLocalBasis(Psi1))
      {
         std::cerr << "mp-ioverlap: fatal: local basis for LHS wavefunction does not match!\n";
         return 1;
      }

      // The default UnitCellSize for output is the wavefunction size
      if (UnitCellSize == 0)
         UnitCellSize = Size;

      // get the list of quantum number sectors
      std::set<QuantumNumber> Sectors;

      if (!Sector.empty())
      {
         for (std::vector<std::string>::const_iterator I = Sector.begin(); I != Sector.end(); ++I)
         {
            Sectors.insert(QuantumNumber(Psi1.GetSymmetryList(),*I));
         }
      }
      else
      {
         // auto-detect the quantum number sectors
         std::set<QuantumNumber> B1 = QuantumNumbersInBasis(Psi1.Basis1().Basis());

         // Merge B2 with the operator basis
         std::set<QuantumNumber> B2 = QuantumNumbersInBasis(Psi2.Basis1().Basis());
         std::set<QuantumNumber> OpBasis = QuantumNumbersInBasis(StringOp.Basis1());
         std::set<QuantumNumber> OpB2;
         for (QuantumNumber const& q2 : B2)
         {
            for (QuantumNumber const& qOp : OpBasis)
            {
               transform_targets(q2, qOp, std::inserter(OpB2, OpB2.begin()));
            }
         }

         // finally determine the target quantum numbers as B1 * target = B2
         for (QuantumNumber const& qi : B1)
         {
            for (QuantumNumber const& qj : OpB2)
            {
               inverse_transform_targets(qj, qi, std::inserter(Sectors, Sectors.begin()));
            }
         }
      }

      if (!Quiet)
      {
         std::cout << "#quantities are calculated per unit cell size of " << UnitCellSize
                   << (UnitCellSize == 1 ? " site\n" : " sites\n");
         if (Size != UnitCellSize)
         {
            std::cout << "#actual size of wavefunction is " << Size << '\n';
         }
         std::cout << "#sector     ";
         if (NumEigen > 1)
            std::cout << "#n     ";
         if (ShowRealPart)
            std::cout << "#real                   ";
         if (ShowImagPart)
            std::cout << "#imag                   ";
         if (ShowCorrLength)
            std::cout << "#corr_length            ";
         if (ShowRate)
            std::cout << "#rate                   ";
         if (ShowMagnitude)
            std::cout << "#magnitude              ";
         if (ShowArgument)
            std::cout << "#argument" << (ShowRadians ? "(rad)" : "(deg)") << "          ";
         std::cout << '\n';
      }
      std::cout << std::left;

      double ScaleFactor = 0;  // sentinel value - we set this once we get the acual length from overlap_arpack

      // Calculate the actual overlaps
      std::vector<TransEigenInfo> EigenList;
      for (std::set<QuantumNumber>::const_iterator I = Sectors.begin(); I != Sectors.end(); ++I)
      {
         //BasicFiniteMPO StringOp = BasicFiniteMPO::make_identity(ExtractLocalBasis(Psi2.Psi));
         std::vector<std::complex<double>> Eigen;
         int Length;
         std::tie(Eigen, Length) = overlap_arpack(Psi1, StringOp, Psi2, NumEigen, *I, Iter, Tol, Verbose);
         ScaleFactor = double(UnitCellSize) / double(Length);

         if (Sort)
         {
            int n = 0;
            for (auto const& e : Eigen)
            {
               EigenList.emplace_back(*I, n++, e);
            }
         }
         else
         {
            int n = 0;
            for (auto const& e : Eigen)
            {
               PrintFormat(*I, n++, e, NumEigen, ShowRealPart, ShowImagPart,
                  ShowCorrLength, ShowRate, ShowMagnitude, ShowArgument, ShowRadians, ScaleFactor);
            }

         }
      }

      if (Sort)
      {
         std::sort(EigenList.begin(), EigenList.end());
         for (auto const& e : EigenList)
         {
            PrintFormat(e.q, e.n, e.x, NumEigen, ShowRealPart, ShowImagPart,
               ShowCorrLength, ShowRate, ShowMagnitude, ShowArgument, ShowRadians, ScaleFactor);
         }
      }


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
