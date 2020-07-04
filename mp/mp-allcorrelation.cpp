// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-iunitcorrelation.cpp
//
// Copyright (C) 2004-2020 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "lattice/latticesite.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "mp/copyright.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "common/prog_options.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "lattice/unitcell-parser.h"
#include "wavefunction/operator_actions.h"
#include "lattice/infinite-parser.h"



namespace prog_opt = boost::program_options;

void PrintFormat(std::complex<double> const& Value, bool ShowRealPart, bool ShowImagPart,
                 bool ShowMagnitude, bool ShowArgument,
                 bool ShowRadians)
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
      double Arg =  std::atan2(Value.imag(), Value.real());
      if (!ShowRadians)
         Arg *= 180.0 / math_const::pi;
      std::cout << std::setw(20) << Arg << "    ";
   }
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
}

int main(int argc, char** argv)
{
   try
   {
      bool ShowRealPart = false, ShowImagPart = false, ShowMagnitude = false;
      bool ShowCartesian = false, ShowPolar = false, ShowArgument = false;
      bool ShowRadians = false;
      std::string PsiStr;
      std::string Op1Str, Op2Str;
      int Verbose = 0;
      bool Quiet = false;
      bool UseTempFile = false;
      std::string LatticeFile;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("cart,c", prog_opt::bool_switch(&ShowCartesian),
          "show the result in cartesian coordinates [equivalent to --real --imag]")
         ("polar,p", prog_opt::bool_switch(&ShowPolar),
          "show the result in polar coodinates [equivalent to --mag --arg]")
         ("real,r", prog_opt::bool_switch(&ShowRealPart),
          "display only the real part of the result")
         ("imag,i", prog_opt::bool_switch(&ShowImagPart),
          "display only the imaginary part of the result")
         ("mag,m", prog_opt::bool_switch(&ShowMagnitude),
          "display the magnitude of the result")
         ("arg,a", prog_opt::bool_switch(&ShowArgument),
          "display the argument (angle) of the result")
         ("radians", prog_opt::bool_switch(&ShowRadians),
          "display the argument in radians instead of degrees")
         ("lattice,l", prog_opt::value(&LatticeFile),
          "Use this lattice file, instead of specifying the lattice for each operator")
         ("tempfile", prog_opt::bool_switch(&UseTempFile),
          "a temporary data file for workspace (path set by environment MP_BINPATH)")
         ("quiet", prog_opt::bool_switch(&Quiet),
          "don't show the preamble and column headings")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose),
          "Extra debug output [can be used multiple times])")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value(&PsiStr), "psi")
         ("op1", prog_opt::value(&Op1Str), "op1")
         ("op2", prog_opt::value(&Op2Str), "op2")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);
      p.add("op1", 1);
      p.add("op2", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("op2") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi> <operator1> <operator2>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (!Quiet)
         print_preamble(std::cout, argc, argv);

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

      // Load the wavefunction
      pvalue_ptr<MPWavefunction> PsiPtr;
      if (UseTempFile)
      {
          mp_pheap::InitializeTempPHeap(Verbose);
          PsiPtr = pheap::ImportHeap(PsiStr);
      }
      else
      {
         PsiPtr = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);
      }
      InfiniteWavefunctionLeft Psi = PsiPtr->get<InfiniteWavefunctionLeft>();

      // If a lattice was explicitly specified, then load it
      pvalue_ptr<InfiniteLattice> Lattice;
      if (!LatticeFile.empty())
         Lattice = pheap::ImportHeap(LatticeFile);

      // Load the operators
      UnitCellMPO Op1;
      int UnitCellSize;
      {
         if (Lattice)
         {
            Op1 = ParseUnitCellOperator(Lattice->GetUnitCell(), 0, Op1Str);
            UnitCellSize = Lattice->GetUnitCell().size();
         }
         else
         {
            InfiniteLattice Lattice1;
            std::tie(Op1, Lattice1) = ParseUnitCellOperatorAndLattice(Op1Str);
            UnitCellSize = Lattice1.GetUnitCell().size();
         }
      }

      UnitCellMPO Op2;
      int UnitCellSize2;
      {
         if (Lattice)
         {
            Op2 = ParseUnitCellOperator(Lattice->GetUnitCell(), 0, Op2Str);
            UnitCellSize2 = Lattice->GetUnitCell().size();
         }
         else
         {
            InfiniteLattice Lattice2;
            std::tie(Op2, Lattice2) = ParseUnitCellOperatorAndLattice(Op2Str);
            UnitCellSize2 = Lattice2.GetUnitCell().size();
         }
      }

      if (adjoint(Op1.TransformsAs()) != Op2.TransformsAs())
      {
         std::cerr << "mp-iunitcorrelation: fatal: product of Op1 and Op2 is not a scalar!\n";
         return 1;
      }

      int const PsiSize = Psi.size();
      int const NumUnitCells = PsiSize / UnitCellSize;

      // Some sanity checks
      if (NumUnitCells == 1)
      {
         std::cerr << "mp-iunitcorrelation: warning: wavefunction unit cell is the same size as lattice; no output.\n";
         return 1;
      }
      if (Op1.size() != UnitCellSize || Op1.offset() != 0)
      {
         std::cerr << "mp-iunitcorrelation: error: operator1 must have support over a single unit cell.\n";
      }
      if (Op2.size() != UnitCellSize || Op2.offset() != 0)
      {
         std::cerr << "mp-iunitcorrelation: error: operator1 must have support over a single unit cell.\n";
      }
      if (UnitCellSize2 != UnitCellSize)
      {
         std::cerr << "mp-iunitcorrelation: error: operator1 and operator2 have different unit cell sizes!\n";
      }



      // show the heading
      if (!Quiet)
      {
         std::cout << "#x    #y    ";
         ShowHeading(ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
         std::cout << std::endl;
      }

      int const LeftUnitCells = NumUnitCells / 2;
      int const RightUnitCells = NumUnitCells - LeftUnitCells;

      // TODO: OP1 will be Op1 * identity (in the Op2.TransformsAs() string)
      // and we will also have that identity as the 'string' term between the MPO's
      GenericMPO OP1 = Op1.MPO();
      GenericMPO OP2 = Op2.MPO();

      // Assemble the E-matrices
      std::vector<MatrixOperator> E;
      E.reserve(LeftUnitCells);

      auto it = Psi.begin();
      for (int i = 0; i < LeftUnitCells; ++i)
      {
         // LinearWavefunction for the next unit cell
         auto next = it+UnitCellSize;
         LinearWavefunction ThisCell = LinearWavefunction::FromContainer(it, next);
         it = next;
         // for existing terms, construct the right operator and the final expectation value
         MatrixOperator Rho = Psi.lambda((i+1)*UnitCellSize);
         Rho = Rho*Rho;
         for (int j = 0; j < E.size(); ++j)
         {
            MatrixOperator e = E[j];
            e = inject_left(e, ThisCell, OP2, ThisCell);
            // TODO: we could make this slightly more efficient by closing off the last contraction with rho
            std::complex<double> v = inner_prod(Rho, e);
            std::cout << std::left << std::setw(5) << j << ' '
                      << std::left << std::setw(5) << i << ' ';
            PrintFormat(v, ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
            std::cout << '\n';
            std::cout << std::flush;
         }

         // extend the operators
         for (int j = 0; j < E.size(); ++j)
         {
            E[j] = inject_left(E[j], ThisCell, ThisCell);
         }
         // next site
         E.push_back(inject_left(MatrixOperator::make_identity(ThisCell.Basis1()), ThisCell, OP1, ThisCell));

      }

      // Assemble the F-matrices
      std::vector<MatrixOperator> F;
      F.reserve(RightUnitCells);
      it = Psi.end();
      for (int i = NumUnitCells-1; i >= LeftUnitCells; --i)
      {
         auto next = it-UnitCellSize;
         LinearWavefunction ThisCell = LinearWavefunction::FromContainer(next, it);
         it = next;
         // for existing terms, construct the right operator and the final expectation value
         MatrixOperator Rho = Psi.lambda((i+1)*UnitCellSize);
         Rho = Rho*Rho;
         for (int j = 0; j < F.size(); ++j)
         {
            MatrixOperator f = F[j];
            // TODO: we could make this slightly more efficient by closing off the last contraction with a trace
            f = inject_right(f, ThisCell, OP1, ThisCell);
            using std::conj;
            std::complex<double> v = conj(trace(f));
            std::cout << std::left << std::setw(5) << i << ' '
                      << std::left << std::setw(5) << (NumUnitCells-j-1) << ' ';
            PrintFormat(v, ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
            std::cout << '\n';
            std::cout << std::flush;
         }

         // extend the operators
         for (int j = 0; j < F.size(); ++j)
         {
            F[j] = inject_right(F[j], ThisCell, ThisCell);
         }
         // next site
         F.push_back(inject_right(Rho, ThisCell, OP2, ThisCell));
      }

      // now the matrix of operators at the half-way point
      for (int i = 0; i < E.size(); ++i)
      {
         for (int j = 0; j < F.size(); ++j)
         {
            std::complex<double> v = inner_prod(F[j], E[i]);
            std::cout << std::left << std::setw(5) << i << ' '
                      << std::left << std::setw(5) << (NumUnitCells-j-1) << ' ';
            PrintFormat(v, ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
            std::cout << '\n';
            std::cout << std::flush;
         }
      }

      pheap::Shutdown();

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
