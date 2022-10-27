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
   PackStateComponent Pack(StringOp.Basis2(), Psi1.Basis2(), Psi2.Basis2());
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
   PackStateComponent Pack(StringOp.Basis1(), Psi1.Basis1(), Psi2.Basis1());
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
get_transfer_eigenpair(InfiniteWavefunctionLeft const& Psi1, InfiniteWavefunctionLeft const& Psi2,
                       QuantumNumber const& q, double tol, int Verbose)
{
   return get_transfer_eigenpair(get_left_canonical(Psi1).first, get_left_canonical(Psi2).first, Psi1.qshift(), ProductMPO::make_identity(ExtractLocalBasis(Psi1), q), tol, Verbose);
}
