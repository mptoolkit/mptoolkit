// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-ioverlap.cpp
//
// Copyright (C) 2004-2024 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2017 Tomohiro <tomohiro.hashizume@uq.net.au>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

void PrintFormat(QuantumNumber const& q, int n, std::complex<double> x,
                 bool ShowSector, bool ShowNum, bool ShowRealPart, bool ShowImagPart,
                 bool ShowCorrLength, bool ShowRate, bool ShowMagnitude, bool ShowArgument,
                 bool ShowRadians, double ScaleFactor)
{
   std::string SectorStr = boost::lexical_cast<std::string>(q);
   std::complex<double> Value = std::pow(x, ScaleFactor);
   if (ShowSector)
      std::cout << std::setw(11) << SectorStr << ' ';
   if (ShowNum)
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
      // NOTE: we don't negate the log term here
      std::cout << std::setw(20) << (std::log(LinearAlgebra::norm_frob(Value)))
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
      bool Sort = false;
      bool Quiet = false;
      bool Reflect = false;
      bool Conj = false;
      bool Print = false;
      int CoarseGrain1 = 1;
      int CoarseGrain2 = 1;
      int NumEigen = 1; // number of eigenvalues to calculate
      bool Scale = false; // scale the overlap by 1 / sqrt(<psi1|psi1> <psi2|psi2>)
      bool ShowSector = false;
      bool NoShowSector = false;
      bool ShowNum = false;
      bool NoShowNum = false;
      bool OneLine = false;
      bool PadZero = false;
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
          "display the rate function log(magnitude) (no minus sign!)")
         ("unitcell,u", prog_opt::value(&UnitCellSize),
          "scale the results to use this unit cell size [default wavefunction unit cell]")
         ("scale", prog_opt::bool_switch(&Scale), "scale the results by the wavefunction norms 1/sqrt(<psi1|psi1><psi2|psi2>)")
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
         // ("iter", prog_opt::value(&Iter),
         //  FormatDefault("Maximum subspace size in the Arnoldi basis", Iter).c_str())
         ("quiet", prog_opt::bool_switch(&Quiet), "don't show the column headings")
         ("showsector", prog_opt::bool_switch(&ShowSector), "Always show the 'sector' column")
         ("no-showsector", prog_opt::bool_switch(&NoShowSector), "Never show the 'sector' column")
         ("shownum", prog_opt::bool_switch(&ShowNum), "Always show the 'n' column (eigenvalue number)")
         ("no-shownum", prog_opt::bool_switch(&NoShowNum), "Never show the 'n' column (eigenvalue number)")
         ("oneline", prog_opt::bool_switch(&OneLine), "Show all output on one line")
         ("pad-zero", prog_opt::bool_switch(&PadZero), "If there are less than n eigenvalues, then pad with zeros")
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

      if (ShowSector && NoShowSector) // if both are specified, revert to default
      {
         ShowSector = NoShowSector = false;
      }
      if (ShowNum && NoShowNum) // if both are specified, revert to default
      {
         ShowNum = NoShowNum = false;
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
            std::cerr << "Wavefunction 1 size is not a multiple of the coarsegrain size, expanding...\n";
            Psi1 = repeat(Psi1, statistics::lcm(Psi1.size(), CoarseGrain1) / Psi1.size());
         }

         Psi1 = coarse_grain(Psi1, CoarseGrain1);
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
            if (Psi2.size() % CoarseGrain2 != 0)
            {
               std::cerr << "Wavefunction 2 size is not a multiple of the coarsegrain size, expanding...\n";
               Psi2 = repeat(Psi1, statistics::lcm(Psi2.size(), CoarseGrain2) / Psi2.size());
            }

            Psi2 = coarse_grain(Psi2, CoarseGrain2);
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
      if (Sectors.size() > 1 && !NoShowSector)
         ShowSector = true;

      if (NumEigen > 1 && !NoShowNum)
         ShowNum = true;

      if (!Quiet)
      {
         std::cout << "#quantities are calculated per unit cell size of " << UnitCellSize
                   << (UnitCellSize == 1 ? " site\n" : " sites\n");
         if (Size != UnitCellSize)
         {
            std::cout << "#actual size of wavefunction is " << Size << '\n';
         }
         std::cout << "#psi1 log amplitude per unit cell size is "
            << Psi1.log_amplitude() * (double(UnitCellSize)/Psi2.size()) << '\n';
         std::cout << "#psi2 log amplitude per unit cell size is "
            << Psi2.log_amplitude() * (double(UnitCellSize)/Psi2.size()) << '\n';
         if (Scale)
            std::cout << "#overlap will be scaled by 1/amplitude\n";
         if (!ShowSector && Sectors.size() == 1)
            std::cout << "#quantum number sector is " << *Sectors.begin() << '\n';
         if (!OneLine)
         {
            if (ShowSector)
               std::cout << "#sector     ";
            if (ShowNum)
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
      }
      std::cout << std::left;

      double ScaleFactor = 0;  // sentinel value - we set this once we get the acual length from overlap()

      // Calculate the actual overlaps
      std::vector<TransEigenInfo> EigenList;
      for (std::set<QuantumNumber>::const_iterator I = Sectors.begin(); I != Sectors.end(); ++I)
      {
         if (Verbose >= 2)
            std::cerr << "Evaluating sector " << (*I) << '\n';
         //BasicFiniteMPO StringOp = BasicFiniteMPO::make_identity(ExtractLocalBasis(Psi2.Psi));
         std::vector<std::complex<double>> Eigen;
         int Length;
         std::tie(Eigen, Length) = overlap(Psi1, StringOp, Psi2, NumEigen, *I, !Scale, Tol, Verbose);
         ScaleFactor = double(UnitCellSize) / double(Length);
         while (PadZero && Eigen.size() < NumEigen)  // pad the eigenvalues to NumEigen, if requested
         {
            Eigen.push_back({0.0, 0.0});
         }

         if (Sort || OneLine)
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
               PrintFormat(*I, n++, e, ShowSector, ShowNum, ShowRealPart, ShowImagPart,
                  ShowCorrLength, ShowRate, ShowMagnitude, ShowArgument, ShowRadians, ScaleFactor);
               std::cout << std::endl;
            }

         }
      }

      if (Sort || OneLine)
      {
         if (Sort)
            std::sort(EigenList.begin(), EigenList.end());

         if (OneLine && !Quiet)
         {
            int i = 0;
            for (auto const& e : EigenList)
            {
               std::string Suffix = "_" + std::to_string(i++);
               if (ShowSector)
                  std::cout << std::setw(12) << std::string("#sector" + Suffix);
               if (ShowNum)
                  std::cout << std::setw(7) << std::string("#n" + Suffix);
               if (ShowRealPart)
                  std::cout << std::setw(24) << std::string("#real" + Suffix);
               if (ShowImagPart)
                  std::cout << std::setw(24) << std::string("#imag" + Suffix);
               if (ShowCorrLength)
                  std::cout << std::setw(24) << std::string("#corr_length" + Suffix);
               if (ShowRate)
                  std::cout << std::setw(24) << std::string("#rate" + Suffix);
               if (ShowMagnitude)
                  std::cout << std::setw(24) << std::string("#magnitude" + Suffix);
               if (ShowArgument)
                  std::cout << std::setw(24) << std::string("#argument" + Suffix + (ShowRadians ? "(rad)" : "(deg)"));
            }
            std::cout << '\n';
         }
         for (auto const& e : EigenList)
         {
            PrintFormat(e.q, e.n, e.x, ShowSector, ShowNum, ShowRealPart, ShowImagPart, ShowCorrLength, ShowRate, ShowMagnitude, ShowArgument, ShowRadians, ScaleFactor);
            if (!OneLine)
               std::cout << '\n';
         }
         if (OneLine)
            std::cout << std::endl;
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
