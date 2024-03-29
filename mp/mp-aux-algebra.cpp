// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-aux-algebra.cpp
//
// Copyright (C) 2014-2021 Ian McCulloch <ian@qusim.net>
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
#include "mp-algorithms/transfer.h"
#include "interface/inittemp.h"
#include "lattice/latticesite.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "common/unique.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "lattice/unitcell.h"
#include "lattice/infinite-parser.h"
#include <boost/algorithm/string/predicate.hpp>
#include "common/statistics.h"
#include <tuple>
#include "parser/matrix-parser.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string PsiStr;
      std::string LatticeFile;
      double Tol = 1e-14;
      int Verbose = 0;
      bool Quiet = false;
      std::vector<std::string> OperatorStr;
      std::vector<std::string> CommutatorStr;
      bool UseTempFile = false;
      std::string QSector;
      bool Square = false;
      std::vector<std::string> Expressions;
      bool NoRandomPhase = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("wavefunction,w", prog_opt::value(&PsiStr), "Wavefunction [required]")
         ("lattice,l", prog_opt::value(&LatticeFile), "use this lattice file for the operators")
         ("sector,s", prog_opt::value(&QSector), "select a different quantum number sector [don't use this unless you know what you are doing]")
         ("square", prog_opt::bool_switch(&Square), "calculate also the square of the commutator [for debugging/information]")
         //      ("commutator,c", prog_opt::value(&CommutatorStr),
         //       "calculate the commutator phase angle, U X X^\\dagger = exp(i*theta) X")
         ("expression", prog_opt::value(&Expressions),
          "Evaluate this expression of the matrices")
         ("tempfile", prog_opt::bool_switch(&UseTempFile),
          "a temporary data file for workspace (path set by environment MP_BINPATH)")
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Tolerance of the Arnoldi eigensolver", Tol).c_str())
         ("norandomize", prog_opt::bool_switch(&NoRandomPhase), "don't randomize the phase of the eigenvectors")
         ("quiet,q", prog_opt::bool_switch(&Quiet), "suppress informational preamble about each operator")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("operator", prog_opt::value(&OperatorStr), "operator")
         ;

      prog_opt::positional_options_description p;
      p.add("operator", -1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("wavefunction") < 1)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " -w <psi> [options] Operator1 [Operator2] ...\n";
         std::cerr << desc << '\n';
         std::cerr << "If -l [--lattice] is specified, then the operators must all come from the specified lattice file\n";
         std::cerr << "Otherwise all operators must be of the form lattice:operator\n";
         std::cerr << "The operators must be of the ProductMPO form.\n";
         std::cerr << "\nThis tool calculates the commutator phase of operator pairs <X Y X\u2020 Y\u2020>,\n";
         std::cerr << "and <X X*> phase.\n";
         std::cerr << "For complex conjugation, prefix the operator expression with c&\n";
         std::cerr << "For spatial reflection, prefix with r& (cr& or rc& for conjugate-reflection)\n";
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (!Quiet)
         print_preamble(std::cout, argc, argv);

      if (!NoRandomPhase)
         srand(ext::get_unique() % RAND_MAX);

      // if the number of eigenvalues is specified but
      // the cutoff is not, then set the cutoff to zero
      //      if (vm.count("num-eigenvalues") == 1 && vm.count("eigen-cutoff") == 0)
      //      {
      //         EigenCutoff = 0;
      //      }


      // Load the wavefunction
      pvalue_ptr<MPWavefunction> Psi;
      if (UseTempFile)
      {
          mp_pheap::InitializeTempPHeap(Verbose);
          Psi = pheap::ImportHeap(PsiStr);
      }
      else
      {
         Psi = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);
      }
      InfiniteWavefunctionLeft InfPsi = Psi->get<InfiniteWavefunctionLeft>();

      // Load the lattice, if it was specified
      pvalue_ptr<InfiniteLattice> Lattice = pheap::ImportHeap(LatticeFile);

      UnitCell Cell = Lattice->GetUnitCell();
      int UnitCellSize = Cell.size();
      //      int const NumUnitCells = InfPsi.size() / UnitCellSize;

      // orthogonalize the wavefunction
      LinearWavefunction Psi1;
      RealDiagonalOperator D;
      std::tie(Psi1, D) = get_left_canonical(InfPsi);
      MatrixOperator Rho = D;
      Rho = scalar_prod(Rho, herm(Rho));
      MatrixOperator Identity = MatrixOperator::make_identity(Psi1.Basis1());

      // reflected and conjugated versions of the wavefunction - we leave them as null
      // until they are actually needed
      LinearWavefunction PsiR, PsiC, PsiRC;

      // The list of U matrices, for each operator
      std::vector<MatrixOperator> U;

      for (unsigned i = 0; i < OperatorStr.size(); ++i)
      {
         std::string OpStr = OperatorStr[i];
         LinearWavefunction* Psi2 = &Psi1;
         // Do we have time reversal or reflection?
         if (boost::starts_with(OpStr, "r&"))
         {
            OpStr = std::string(OpStr.begin()+2, OpStr.end());
            if (PsiR.empty())
            {
               InfiniteWavefunctionLeft PR = InfPsi;
               inplace_reflect(PR);
               // The wavefunction must be in a reflection-symmetric basis, or this isn't valid
               if (PR.Basis1() != InfPsi.Basis1())
               {
                  std::cerr << "mp-aux-algebra: cannot reflect operator because the basis is not reflection symmetric.\n"
                            << "mp-aux-algebra: ignoring operator " << OpStr << '\n';
                  continue;
               }
               PsiR = get_left_canonical(PR).first;
            }
            Psi2 = &PsiR;
         }
         else if (boost::starts_with(OpStr,"c&"))
         {
            OpStr = std::string(OpStr.begin()+2, OpStr.end());
            if (PsiC.empty())
            {
               PsiC = conj(Psi1);
            }
            Psi2 = &PsiC;
         }
         else if (boost::starts_with(OpStr,"rc&") || boost::starts_with(OpStr,"cr&"))
         {
            OpStr = std::string(OpStr.begin()+3, OpStr.end());
            if (PsiR.empty())
            {
               InfiniteWavefunctionLeft PR = InfPsi;
               inplace_reflect(PR);

               // The wavefunction must be in a reflection-symmetric basis, or this isn't valid
               if (PR.Basis1() != InfPsi.Basis1())
               {
                  std::cerr << "mp-aux-algebra: cannot reflect operator because the basis is not reflection symmetric.\n"
                            << "mp-aux-algebra: ignoring operator " << OpStr << '\n';
                  continue;
               }
               PsiR = get_left_canonical(PR).first;
            }
            if (PsiRC.empty())
            {
               PsiRC = conj(PsiR);
            }
            Psi2 = &PsiRC;
         }

         ProductMPO StringOperator = ParseProductOperator(*Lattice, OpStr);

         // Hack for non-identity symmetry sector
         if (!QSector.empty())
         {
            QuantumNumber q(Psi1.GetSymmetryList(), QSector);
            StringOperator = StringOperator * ProductMPO::make_identity(StringOperator.LocalBasis2List(), q);
         }

         int n = Psi1.size(); // unit cell size

         std::complex<double> e;
         MatrixOperator v;
         std::tie(e, v) = get_left_transfer_eigenvector(*Psi2, Psi1, InfPsi.qshift(), StringOperator,
                                                            Tol, Verbose);


         // Normalization
         // it might not be unitary, eg anti-unitary.  So we need to take the 4th power
         MatrixOperator v2 = scalar_prod(herm(v), v);
         v2 = scalar_prod(herm(v2), v2);
         std::complex<double> x = inner_prod(v2, Rho);

         v *= std::sqrt(std::sqrt(1.0 / x));
         // randomize phase

         if (!NoRandomPhase)
         {
            v *= std::polar(1.0, LinearAlgebra::random<double>() * 2 * math_const::pi);
         }

         //v *= std::sqrt(Dim); // make it properly unitary
         U.push_back(v);

         if (!Quiet)
         {
            if (i > 0)
               std::cout << '\n';
            std::cout << "#Operator " << i << " = " << OperatorStr[i] << '\n'
                      << "#unit cell size " << n << '\n'
                      << "#eigenvalue = " << e << '\n'
                      << "#magnitude = " << std::abs(e) << '\n';
            if (n != int(Psi1.size()))
            {
               std::cout << "#magnitude per wavefunction unit cell size (" << Psi1.size() << " sites) = "
                         << std::pow(std::abs(e), double(Psi1.size())/double(n)) << '\n';
            }
            if (n != 1)
            {
               std::cout << "#magnitude per site = "
                         << std::pow(std::abs(e), 1.0/double(n)) << '\n';
            }

            std::cout << "#UU\u2020 = " << inner_prod(Rho, scalar_prod(v,herm(v))) << "\n";
            //std::cout << "#<U> = " << inner_prod(Rho, v) << "\n";

            if (is_scalar(v.TransformsAs()) && v.Basis2() == v.Basis1())
            {
               std::cout << "#UU* = " << inner_prod(Rho, v*conj(v)) << '\n';
               std::complex<double> U2 =  inner_prod(Rho, v*v);
               std::cout << "#U^2 = " << U2 << '\n';
               std::cout << "#U^2 magnitude = " << std::abs(U2) << '\n';
            }
         }

         // The eigenvalue should be nearly 1, or it isn't a unitary operator
         if (LinearAlgebra::norm_frob(std::pow(LinearAlgebra::norm_frob(e),1.0/n)-1.0) > 1E-3)
         {
            std::cerr << basename(argv[0]) << ": warning: operator eigenvalue modulus per site "
                      << std::pow(LinearAlgebra::norm_frob(e),1.0/n)
                      << " differs from 1.0.  Results may be unreliable!\n";
         }
      }
      if (!Quiet)
         std::cout << '\n';

      // Now go through each operator pair
      if (!Quiet)
      {
         if (Square)
         {
            std::cout << "#Op1  #Op2  #Commutator-Real  #Commutator-Imag      #Square-Real     #Square-Imag\n";
         }
         else
         {
            std::cout << "#Op1  #Op2  #Commutator-Real  #Commutator-Imag\n";
         }
      }
      for (unsigned i = 0; i < U.size(); ++i)
      {
         for (unsigned j = i+1; j < U.size(); ++j)
         {
            MatrixOperator X = herm(U[j]) * (herm(U[i]) * U[j] * U[i]);
//            scalar_prod(herm(U[j]), operator_prod(herm(U[i]), U[j], U[i]));
            //            std::complex<double> x = inner_prod(U[i]*U[j], U[j]*U[i]) / Dim;
            std::complex<double> x = inner_prod(X, Rho);
            // scalar_prod(herm(U[i]*U[j]), U[j]*U[i]), Rho);
            //            std::complex<double> tr = trace(U[i]*U[j]*Rho);
            //      TRACE(scalar_prod(herm(U[i]*U[j]), U[i]*U[j]));
            std::cout << std::setw(5) << std::left << i << " "
                      << std::setw(5) << std::left << j << " "
                      << std::setw(17) << std::right << std::fixed << x.real() << " "
                      << std::setw(17) << std::right << std::fixed << x.imag() << " ";

            if (Square)
            {
               x = inner_prod(X*X, Rho);
               std::cout << std::setw(17) << std::right << std::fixed << x.real() << " "
                         << std::setw(17) << std::right << std::fixed << x.imag() << " "
                         << std::endl;
            }
            else
            {
               std::cout << std::endl;
            }
         }
      }

      // Now evaluate expressions, if any
      if (!Expressions.empty())
      {
         std::map<std::string, MatrixOperator> Matrices;
         for (unsigned n = 0; n < U.size(); ++n)
         {
            Matrices[std::string(1, char('A'+n))] = U[n];
         }

         Function::ArgumentList Args;

         for (unsigned i = 0; i < Expressions.size(); ++i)
         {
            MatrixOperator X = ParseMatrixOperator(Expressions[i], Args, Matrices);
            std::cout << "Expression " << Expressions[i] << " = " << inner_prod(X, Rho) << '\n';
         }
      }


#if 0
      // Now calculate the commutator phase angles
      for (unsigned j = 0; j < CommutatorStr.size(); ++j)
      {
         ProductMPO Op = ParseProductOperator(*Lattice, CommutatorStr[j]);

         Op = repeat(Op, Psi.size() / UnitCellSize);
         MatrixOperator O = MatrixOperator::make_identity(Psi.Basis2());
         O = inject_left_qshift(O, Op, Psi, Psi1->shift());

         for (unsigned i = 0; i < U.size(); ++i)
         {

            MatrixOperator Obar = triple_prod(U[i], O, herm(U[i]));
            std::complex<double> Phase = inner_prod(O, Obar) / inner_prod(O, O);
            std::cout << "Phase " << CommutatorStr[j] << " with " << OperatorStr[i] << " = " << Phase << '\n';
         }
      }
#endif

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
