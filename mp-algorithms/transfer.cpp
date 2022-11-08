// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/transfer.h
//
// Copyright (C) 2015-2021 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "transfer.h"
#include "mps/packunpack.h"
#include "mp-algorithms/arnoldi.h"
#include "common/statistics.h"
#include "wavefunction/operator_actions.h"
#include "linearalgebra/arpack_wrapper.h"

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

std::tuple<std::complex<double>, MatrixOperator>
get_right_transfer_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift,
                      ProductMPO const& StringOp,
                      double tol, int Verbose)
{
   int ncv = 0;
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   CHECK_EQUAL(Psi1.size() % StringOp.size(), 0);
   PackStateComponent Pack(StringOp.Basis1(), Psi1.Basis2(), Psi2.Basis2());
   int n = Pack.size();
   //   double tolsave = tol;
   //   int ncvsave = ncv;
   int NumEigen = 1;

   std::vector<std::complex<double>> OutVec;
      LinearAlgebra::Vector<std::complex<double>> LeftEigen =
         LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(Pack,
                                                            RightMultiplyOperator(Psi1, QShift,
                                                                                 StringOp,
                                                                                 Psi2, QShift, Psi1.size(), Verbose-1)),
                                          n, NumEigen, tol, &OutVec, ncv, false, Verbose);

   StateComponent LeftVector = Pack.unpack(&(OutVec[0]));

   return std::make_tuple(LeftEigen[0], LeftVector[0]);
}

std::tuple<std::complex<double>, MatrixOperator>
get_left_transfer_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift,
                     ProductMPO const& StringOp,
                     double tol, int Verbose)
{
   int ncv = 0;
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   CHECK_EQUAL(Psi1.size() % StringOp.size(), 0);
   PackStateComponent Pack(StringOp.Basis1(), Psi1.Basis2(), Psi2.Basis2());
   int n = Pack.size();
   int NumEigen = 1;

   std::vector<std::complex<double>> OutVec;
      LinearAlgebra::Vector<std::complex<double>> LeftEigen =
         LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(Pack, LeftMultiplyOperator(Psi1, QShift, StringOp, Psi2, QShift,
            Psi1.size(), Verbose-1)), n, NumEigen, tol, &OutVec, ncv, false, Verbose);

   StateComponent LeftVector = Pack.unpack(&(OutVec[0]));

   return std::make_tuple(LeftEigen[0], LeftVector[0]);
}

std::tuple<std::vector<std::complex<double>>, std::vector<MatrixOperator>>
get_left_transfer_eigenvectors(int N, LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift,
                     ProductMPO const& StringOp,
                     double tol, int Verbose)
{
   int ncv = 0;
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   CHECK_EQUAL(Psi1.size() % StringOp.size(), 0);
   PackStateComponent Pack(StringOp.Basis1(), Psi1.Basis2(), Psi2.Basis2());
   int n = Pack.size();
   int NumEigen = N;

   std::vector<std::complex<double>> OutVec;
      LinearAlgebra::Vector<std::complex<double>> LeftEigen =
         LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(Pack, LeftMultiplyOperator(Psi1, QShift, StringOp, Psi2, QShift,
            Psi1.size(), Verbose-1)), n, NumEigen, tol, &OutVec, ncv, false, Verbose);

   std::vector<MatrixOperator> LeftVectors(N);
   for (int i = 0; i < N; ++i)
   {
      StateComponent C = Pack.unpack(&(OutVec[n*i]));
      LeftVectors[i] = C[0];
   }
   std::vector<std::complex<double>> Eigen(LeftEigen.begin(), LeftEigen.end());
   return std::make_tuple(Eigen, LeftVectors);

}

std::tuple<std::complex<double>, MatrixOperator, MatrixOperator>
get_transfer_eigenpair(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift,
                       ProductMPO const& StringOp, double tol, int Verbose)
{
   std::complex<double> eLeft, eRight;
   MatrixOperator Left, Right;
   std::tie(eLeft, Left) = get_left_transfer_eigenvector(Psi1, Psi2, QShift, StringOp, tol, Verbose);
   std::tie(eRight, Right) = get_right_transfer_eigenvector(Psi1, Psi2, QShift, StringOp, tol, Verbose);
   CHECK(norm_frob(eLeft-conj(eRight)) < 1E-10)(eLeft)(eRight);
   // normalize
   Right *= 1.0 / norm_frob(Right);
   Left *= 1.0 / inner_prod(delta_shift(Right, QShift), Left);
   return std::make_tuple(eLeft, Left, Right);
}

std::tuple<std::complex<double>, MatrixOperator, MatrixOperator>
get_transfer_eigenpair(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, QuantumNumber const& q, double tol, int Verbose)
{
   return get_transfer_eigenpair(Psi1, Psi2, QShift, ProductMPO::make_identity(ExtractLocalBasis(Psi1), q), tol, Verbose);
}

std::tuple<std::complex<double>, MatrixOperator, MatrixOperator>
get_transfer_eigenpair(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, double tol, int Verbose)
{
   return get_transfer_eigenpair(Psi1, Psi2, QShift, QuantumNumbers::QuantumNumber(Psi1.GetSymmetryList()), tol, Verbose);
}

std::tuple<std::complex<double>, MatrixOperator, MatrixOperator>
get_transfer_eigenpair(InfiniteWavefunctionLeft const& Psi1, InfiniteWavefunctionLeft const& Psi2,
                       QuantumNumber const& q, double tol, int Verbose)
{
   return get_transfer_eigenpair(get_left_canonical(Psi1).first, get_left_canonical(Psi2).first, Psi1.qshift(), ProductMPO::make_identity(ExtractLocalBasis(Psi1), q), tol, Verbose);
}


LinearAlgebra::Vector<std::complex<double>>
get_spectrum_string(LinearWavefunction const& Psi, QuantumNumber const& QShift,
                    ProductMPO const& StringOp,
                    int NumEigen, double tol,
                    LinearAlgebra::Vector<MatrixOperator>* LeftVectors,
                    LinearAlgebra::Vector<MatrixOperator>* RightVectors,
                    int ncv, bool Sort, int Verbose)
{
   PackStateComponent PackL(StringOp.Basis1(), Psi.Basis1(), Psi.Basis1());
   PackStateComponent PackR(StringOp.Basis1(), Psi.Basis2(), Psi.Basis2());
   int n = PackL.size();

   if (Verbose >= 1)
   {
      std::cerr << "Calculating left eigenvalues\n";
   }

   std::vector<std::complex<double>>* OutVecL
      = LeftVectors ? new std::vector<std::complex<double>>() : NULL;
   LinearAlgebra::Vector<std::complex<double>>  LeftEigen =
   LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(PackL, LeftMultiplyOperator(Psi, QShift, StringOp, Psi, QShift,
      Psi.size(), Verbose-1)), n, NumEigen, tol, OutVecL, ncv, Sort, Verbose);

   if (RightVectors)
   {
      if (Verbose >= 1)
      {
         std::cerr << "Calculating right eigenvalues\n";
      }
      std::vector<std::complex<double>>* OutVecR = new std::vector<std::complex<double>>();
      LinearAlgebra::Vector<std::complex<double>> RightEigen =
      LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(PackR, RightMultiplyOperator(Psi, QShift, StringOp, Psi, QShift,
      Psi.size(), Verbose-1)), n, NumEigen, tol, OutVecR, ncv, Sort, Verbose);

      // The right vectors are the hermitian conjugate, not the transpose.  So conjugate eigenvalues
      RightEigen = conj(RightEigen);
      *RightVectors = LinearAlgebra::Vector<MatrixOperator>(RightEigen.size());

      // If we have both left & right eigenvectors, then match the corresponding eigenvalues
      if (LeftVectors)
         LinearAlgebra::MatchEigenvectors(n, LeftEigen, *OutVecL, RightEigen, *OutVecR, tol*10);

      // Unpack the eigenvectors into the output array.
      // note that we do not conjugate the eigenvector, since we want this to be a 'column vector',
      // that we will use on the left-hand side of an inner product (eg inner_prod(left, right)).
      for (unsigned i = 0; i < RightEigen.size(); ++i)
      {
         (*RightVectors)[i] = PackR.unpack(&((*OutVecR)[n*i]))[0];
      }
   }

   // eigenvectors
   if (LeftVectors)
   {
      *LeftVectors = LinearAlgebra::Vector<MatrixOperator>(LeftEigen.size());
      for (unsigned i = 0; i < LeftEigen.size(); ++i)
      {
         (*LeftVectors)[i] = PackL.unpack(&((*OutVecL)[n*i]))[0];
      }
   }

   return LeftEigen;
}
