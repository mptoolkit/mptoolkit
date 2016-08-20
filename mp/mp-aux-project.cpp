// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-aux-project.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "mps/packunpack.h"
#include "interface/inittemp.h"
#include "lattice/latticesite.h"
#include "wavefunction/operator_actions.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "mp-algorithms/arnoldi.h"
#include "lattice/unitcell.h"
#include "lattice/infinite-parser.h"
#include "linearalgebra/arpack_wrapper.h"
#include <boost/algorithm/string/predicate.hpp>
#include "common/statistics.h"
#include <tuple>

namespace prog_opt = boost::program_options;

MatrixOperator
ConstructProjectorOntoPositiveDiagonal(MatrixOperator const& M)
{
   CHECK_EQUAL(M.Basis1(), M.Basis2());

   // Determine which states we should keep
   std::vector<std::vector<int>> KeptStates;
   for (unsigned i = 0; i < M.Basis1().size(); ++i)
   {
      KeptStates.push_back(std::vector<int>());
      MatrixOperator::const_inner_iterator I = iterate_at(M.data(), i,i);
      if (I)
      {
         for (unsigned j = 0; j < size1(*I); ++j)
         {
            if ((*I)(j,j).real() > 0)
               KeptStates[i].push_back(j);
         }
      }
   }

   // Assemble the new basis
   VectorBasis B(M.GetSymmetryList());
   for (unsigned i = 0; i < KeptStates.size(); ++i)
   {
      // Skip over any sectors where we have no kept states
      if (!KeptStates[i].empty())
         B.push_back(M.Basis1()[i], KeptStates[i].size());
   }

   // Assemble the projector
   MatrixOperator U(B, M.Basis1(), QuantumNumber(M.GetSymmetryList()));
   int n = 0; // index into the new basis
   for (unsigned i = 0; i < KeptStates.size(); ++i)
   {
      if (!KeptStates[i].empty())
      {
         LinearAlgebra::Matrix<std::complex<double>> m(KeptStates[i].size(), M.Basis1().dim(i), 0.0);
         for (unsigned k = 0; k < KeptStates[i].size(); ++k)
         {
            m(k,KeptStates[i][k]) = 1.0;
         }
         U(n,i) = m;
         ++n;
      }
   }

   U.check_structure();
   return U;
}


template <typename Func>
struct PackApplyFunc
{
   PackApplyFunc(PackStateComponent const& Pack_, Func f_) : Pack(Pack_), f(f_) {}

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      StateComponent x = Pack.unpack(In);
      x = f(x);
      Pack.pack(x, Out);
   }
   PackStateComponent const& Pack;
   Func f;
};

template <typename Func>
PackApplyFunc<Func>
MakePackApplyFunc(PackStateComponent const& Pack_, Func f_)
{
   return PackApplyFunc<Func>(Pack_, f_);
}

std::tuple<std::complex<double>, int, StateComponent>
get_left_eigenvector(LinearWavefunction const& Psi1, QuantumNumber const& QShift1,
                     LinearWavefunction const& Psi2, QuantumNumber const& QShift2,
                     ProductMPO const& StringOp,
                     double tol = 1E-14, int Verbose = 0)
{
   int ncv = 0;
   int Length = statistics::lcm(Psi1.size(), Psi2.size(), StringOp.size());
   PackStateComponent Pack(StringOp.Basis1(), Psi1.Basis2(), Psi2.Basis2());
   int n = Pack.size();
   //   double tolsave = tol;
   //   int ncvsave = ncv;
   int NumEigen = 1;

   std::vector<std::complex<double> > OutVec;
      LinearAlgebra::Vector<std::complex<double> > LeftEigen =
         LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(Pack,
                                                            LeftMultiplyOperator(Psi1, QShift1,
                                                                                 StringOp,
                                                                                 Psi2, QShift2, Length, Verbose-2)),
                                          n, NumEigen, tol, &OutVec, ncv, false, Verbose);

   StateComponent LeftVector = Pack.unpack(&(OutVec[0]));

   return std::make_tuple(LeftEigen[0], Length, LeftVector);
}

int main(int argc, char** argv)
{
   try
   {
      double Tol = 1e-14;
      int Verbose = 0;
      bool Quiet = false;
      std::string Operator1;
      std::string Operator2;
      bool UseTempFile = false;
      double Sign = 0;
      std::string InputFile;
      std::string OutputFile;
      bool Force = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("tempfile", prog_opt::bool_switch(&UseTempFile),
          "a temporary data file for workspace (path set by environment MP_BINPATH)")
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Tolerance of the Arnoldi eigensolver", Tol).c_str())
         ("force,f", prog_opt::bool_switch(&Force),
          "allow overwriting the output file, if it already exists")
         ("quiet,q", prog_opt::bool_switch(&Quiet), "suppress informational preamble about each operator")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("operator1", prog_opt::value(&Operator1), "operator1")
         ("operator2", prog_opt::value(&Operator2), "operator2")
         ("sign", prog_opt::value(&Sign), "sign")
         ("psi1", prog_opt::value(&InputFile), "psi1")
         ("psi2", prog_opt::value(&OutputFile), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("operator1", 1);
      p.add("operator2", 1);
      p.add("sign", 1);
      p.add("psi1", 1);
      p.add("psi2", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("psi2") < 1)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <operator1> <operator2> <sign> <input-psi> <output-psi>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      // if the number of eigenvalues is specified but
      // the cutoff is not, then set the cutoff to zero
      //      if (vm.count("num-eigenvalues") == 1 && vm.count("eigen-cutoff") == 0)
      //      {
      //         EigenCutoff = 0;
      //      }


      if (Verbose > 0)
         std::cout << "Loading wavefunction..." << std::endl;

      pvalue_ptr<MPWavefunction> PsiPtr;
      if (InputFile == OutputFile)
         PsiPtr = pheap::OpenPersistent(InputFile.c_str(), mp_pheap::CacheSize());
      else
      {
         pheap::Initialize(OutputFile, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);
         PsiPtr = pheap::ImportHeap(InputFile);
      }

      InfiniteWavefunctionLeft Psi = PsiPtr->get<InfiniteWavefunctionLeft>();


      // orthogonalize the wavefunction
      LinearWavefunction Psi1;
      RealDiagonalOperator D;
      std::tie(Psi1, D) = get_left_canonical(Psi);
      MatrixOperator Rho = D;
      Rho = scalar_prod(Rho, herm(Rho));
      MatrixOperator Identity = MatrixOperator::make_identity(Psi1.Basis1());
      double Dim = Psi1.Basis1().total_degree();

      ProductMPO Op1 = ParseProductOperatorAndLattice(Operator1).first;
      ProductMPO Op2 = ParseProductOperatorAndLattice(Operator2).first;

      if (Verbose > 0)
         std::cout << "Calculating eigenmatrix of operator 1..." << std::endl;

      std::complex<double> e1;
      StateComponent v1;
      int n1;
      std::tie(e1, n1, v1) = get_left_eigenvector(Psi1, Psi.qshift(), Psi1, Psi.qshift(), Op1, Tol, Verbose);

      // Normalization
      // it might not be unitary, eg anti-unitary.  So we need to take the 4th power
#if 0
      std::complex<double> x1 = inner_prod(Rho, scalar_prod(v1,herm(v1)));

#else
      std::complex<double> x1 = inner_prod(scalar_prod(herm(v1), operator_prod(herm(v1), v1, v1)), Rho);
#endif
      v1 *= std::sqrt(std::sqrt(1.0 / x1));

      if (!Quiet)
      {
         std::cout << "#Operator " << Operator1 << '\n'
                   << "#eigenvalue = " << e1 << '\n';
      }

      if (Verbose > 0)
         std::cout << "Calculating eigenmatrix of operator 2..." << std::endl;

      std::complex<double> e2;
      StateComponent v2;
      int n2;
      std::tie(e2, n2, v2) = get_left_eigenvector(Psi1, Psi.qshift(), Psi1, Psi.qshift(), Op2, Tol, Verbose);

      // Normalization
#if 0
      std::complex<double> x2 = inner_prod(Rho, scalar_prod(v2,herm(v2)));

#else
      std::complex<double> x2 = inner_prod(scalar_prod(herm(v2), operator_prod(herm(v2), v2, v2)), Rho);
#endif
      v2 *= std::sqrt(std::sqrt(1.0 / x2));

      if (!Quiet)
      {
         std::cout << "#Operator " << Operator2 << '\n'
                   << "#eigenvalue = " << e2 << '\n';
      }

      if (Verbose > 0)
         std::cout << "Calculating commutator..." << std::endl;

      // A^\dagger B^\dagger A B
      MatrixOperator U = scalar_prod(herm(v1), operator_prod(herm(v2), v1, v2));

#if 1
      // construct the projector onto the subspace
      MatrixOperator P = 0.5 * (U + Sign * MatrixOperator::make_identity(U.Basis1()));



#else
      // project onto the parts of U that are positive
      if (Verbose > 0)
         std::cout << "Calculating projector..." << std::endl;
      MatrixOperator P = ConstructProjectorOntoPositiveDiagonal(U);
#endif

      if (!Quiet)
         std::cout << "Number of states in projected basis: " << P.Basis1().total_dimension() << std::endl;

      // apply the projector
      Psi1.set_front(prod(P, Psi1.get_front()));
      Psi1.set_back(prod(Psi1.get_back(), herm(P)));

      // save the wavefunction
      if (Verbose > 0)
         std::cout << "Orthogonalizing wavefunction..." << std::endl;
      PsiPtr.mutate()->Wavefunction() = InfiniteWavefunctionLeft::Construct(Psi1, Psi.qshift(), Verbose);
      PsiPtr.mutate()->AppendHistory(EscapeCommandline(argc, argv));

      if (Verbose > 0)
         std::cout << "Finished." << std::endl;

      pheap::ShutdownPersistent(PsiPtr);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      pheap::Cleanup();
      return 1;
   }
}
