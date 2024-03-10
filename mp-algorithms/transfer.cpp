// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/transfer.cpp
//
// Copyright (C) 2015-2024 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022-2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "transfer.h"
#include "mps/packunpack.h"
#include "mp-algorithms/arnoldi.h"
#include "common/statistics.h"
#include "wavefunction/operator_actions.h"
#include "linearalgebra/arpack_wrapper.h"

// utility functions

typedef std::complex<double> complex;

void swap_vectors(int n, std::complex<double>* x, std::complex<double>* y)
{
   std::vector<std::complex<double>> Temp(x, x+n);
   LinearAlgebra::fast_copy(y, y+n, x);
   LinearAlgebra::fast_copy(&Temp[0], &Temp[0]+n, y);
}

void
MatchEigenvectors(int n,
                  LinearAlgebra::Vector<std::complex<double> >& LeftValues,
                  std::vector<std::complex<double> >& LeftVectors,
                  LinearAlgebra::Vector<std::complex<double> >& RightValues,
                  std::vector<std::complex<double> >& RightVectors, double tol, bool IgnoreFinalMismatch)
{
   CHECK_EQUAL(LeftValues.size(), RightValues.size());

   for (unsigned i = 0; i < LeftValues.size(); ++i)
   {
      // find the right eigenvalue closest to LeftValues[i]
      unsigned Pos = i;
      double Dist = norm_frob(LeftValues[i] - RightValues[Pos]);
      for (unsigned j = i+1; j < RightValues.size(); ++j)
      {
         double jDist = norm_frob(LeftValues[i] - RightValues[j]);
         if (jDist < Dist)
         {
            Pos = j;
            Dist = jDist;
         }
      }
      if (Dist > tol && (i < LeftValues.size()-1 || !IgnoreFinalMismatch))
         std::cerr << "MatchEigenvalues: warning: left & right eigenvalues differ by > tol.  Left="
                   << LeftValues[i] << ", Right=" << RightValues[Pos] << ", delta=" << Dist << '\n';
      std::swap(RightValues[i], RightValues[Pos]);
      swap_vectors(n, &RightVectors[n*i], &RightVectors[n*Pos]);
   }
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

template <typename Func>
struct PackApplyFuncMatrix
{
   PackApplyFuncMatrix(PackMatrixOperator const& Pack_, Func f_) : Pack(Pack_), f(f_) {}

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      MatrixOperator x = Pack.unpack(In);
      x = f(x);
      Pack.pack(x, Out);
   }
   PackMatrixOperator const& Pack;
   Func f;
};

template <typename Func>
PackApplyFuncMatrix<Func>
MakePackApplyFunc(PackMatrixOperator const& Pack_, Func f_)
{
   return PackApplyFuncMatrix<Func>(Pack_, f_);
}

std::tuple<std::complex<double>, MatrixOperator>
get_left_transfer_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, ProductMPO const& StringOp, double tol, int Verbose)
{
   int ncv = 0;
   TRACE(StringOp.Basis1());  // FIXME: this is the error here, the StringOp is invalid
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   CHECK_EQUAL(Psi1.size() % StringOp.size(), 0);
   PackStateComponent Pack(StringOp.Basis1(), Psi1.Basis1(), Psi2.Basis1());
   int n = Pack.size();
   int NumEigen = 1;

   //TRACE(StringOp.Basis1())(Psi1.Basis1())(Psi2.Basis1())(n);
   std::vector<std::complex<double>> OutVec;
      LinearAlgebra::Vector<std::complex<double>> LeftEigen =
         LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(Pack, LeftMultiplyOperator(Psi1, QShift, StringOp, Psi2, QShift,
            Psi1.size(), Verbose-1)), n, NumEigen, nullptr, tol, &OutVec, ncv, false, Verbose);

   StateComponent LeftVector = Pack.unpack(&(OutVec[0]));

   return std::make_tuple(LeftEigen[0], LeftVector[0]);
}

std::tuple<std::complex<double>, MatrixOperator>
get_left_transfer_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, ProductMPO const& StringOp, MatrixOperator InitialGuess, double tol, int Verbose)
{
   int ncv = 0;
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   CHECK_EQUAL(Psi1.size() % StringOp.size(), 0);
   PackStateComponent Pack(StringOp.Basis1(), Psi1.Basis1(), Psi2.Basis1());
   int n = Pack.size();
   int NumEigen = 1;
   std::vector<std::complex<double>> Initial(n);
   PackMatrixOperator(InitialGuess).pack(InitialGuess, Initial.data());

   std::vector<std::complex<double>> OutVec;
      LinearAlgebra::Vector<std::complex<double>> LeftEigen =
         LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(Pack, LeftMultiplyOperator(Psi1, QShift, StringOp, Psi2, QShift,
            Psi1.size(), Verbose-1)), n, NumEigen, Initial.data(), tol, &OutVec, ncv, false, Verbose);

   StateComponent LeftVector = Pack.unpack(&(OutVec[0]));

   return std::make_tuple(LeftEigen[0], LeftVector[0]);
}

std::tuple<std::complex<double>, MatrixOperator>
get_left_transfer_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, double tol, int Verbose)
{
   QuantumNumbers::QuantumNumber q(Psi1.GetSymmetryList());
   return get_left_transfer_eigenvector(Psi1, Psi2, QShift, ProductMPO::make_identity(ExtractLocalBasis(Psi1), q), tol, Verbose);
}

std::tuple<std::complex<double>, MatrixOperator>
get_left_transfer_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, MatrixOperator InitialGuess, double tol, int Verbose)
{
   return get_left_transfer_eigenvector(Psi1, Psi2, QShift, ProductMPO::make_identity(ExtractLocalBasis(Psi1), InitialGuess.TransformsAs()), InitialGuess, tol, Verbose);
}

std::tuple<std::vector<std::complex<double>>, std::vector<MatrixOperator>>
get_left_transfer_eigenvectors(int N, LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, ProductMPO const& StringOp, double tol, int Verbose)
{
   int ncv = 0;
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   CHECK_EQUAL(Psi1.size() % StringOp.size(), 0);
   PackStateComponent Pack(StringOp.Basis1(), Psi1.Basis1(), Psi2.Basis1());
   int n = Pack.size();
   int NumEigen = N;

   std::vector<std::complex<double>> OutVec;
      LinearAlgebra::Vector<std::complex<double>> LeftEigen =
         LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(Pack, LeftMultiplyOperator(Psi1, QShift, StringOp, Psi2, QShift,
            Psi1.size(), Verbose-1)), n, NumEigen, nullptr, tol, &OutVec, ncv, false, Verbose);

   std::vector<MatrixOperator> LeftVectors(N);
   for (int i = 0; i < N; ++i)
   {
      StateComponent C = Pack.unpack(&(OutVec[n*i]));
      LeftVectors[i] = C[0];
   }
   std::vector<std::complex<double>> Eigen(LeftEigen.begin(), LeftEigen.end());
   return std::make_tuple(Eigen, LeftVectors);

}

std::tuple<std::vector<std::complex<double>>, std::vector<MatrixOperator>>
get_left_transfer_eigenvectors(int N, LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, ProductMPO const& StringOp, MatrixOperator InitialGuess, double tol, int Verbose)
{
   int ncv = 0;
   CHECK_EQUAL(InitialGuess.Basis1(), Psi1.Basis1());
   CHECK_EQUAL(InitialGuess.Basis2(), Psi2.Basis1());
   CHECK_EQUAL(InitialGuess.TransformsAs(), StringOp.Basis1()[0]);
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   CHECK_EQUAL(Psi1.size() % StringOp.size(), 0);
   PackStateComponent Pack(StringOp.Basis1(), Psi1.Basis1(), Psi2.Basis1());
   int n = Pack.size();
   int NumEigen = N;
   std::vector<std::complex<double>> Initial(n);
   PackMatrixOperator(InitialGuess).pack(InitialGuess, Initial.data());

   std::vector<std::complex<double>> OutVec;
      LinearAlgebra::Vector<std::complex<double>> LeftEigen =
         LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(Pack, LeftMultiplyOperator(Psi1, QShift, StringOp, Psi2, QShift,
            Psi1.size(), Verbose-1)), n, NumEigen, Initial.data(), tol, &OutVec, ncv, false, Verbose);

   std::vector<MatrixOperator> LeftVectors(N);
   for (int i = 0; i < N; ++i)
   {
      StateComponent C = Pack.unpack(&(OutVec[n*i]));
      LeftVectors[i] = C[0];
   }
   std::vector<std::complex<double>> Eigen(LeftEigen.begin(), LeftEigen.end());
   return std::make_tuple(Eigen, LeftVectors);

}

std::tuple<std::vector<std::complex<double>>, std::vector<MatrixOperator>>
get_left_transfer_eigenvectors(int N, LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, double tol, int Verbose)
{
   return get_left_transfer_eigenvectors(N, Psi1, Psi2, QShift, ProductMPO::make_identity(ExtractLocalBasis(Psi1)), tol, Verbose);
}

std::tuple<std::vector<std::complex<double>>, std::vector<MatrixOperator>>
get_left_transfer_eigenvectors(int N, LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, MatrixOperator InitialGuess, double tol, int Verbose)
{
   return get_left_transfer_eigenvectors(N, Psi1, Psi2, QShift, ProductMPO::make_identity(ExtractLocalBasis(Psi1), InitialGuess.TransformsAs()), InitialGuess, tol, Verbose);
}


//
// right
//

std::tuple<std::complex<double>, MatrixOperator>
get_right_transfer_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, ProductMPO const& StringOp, MatrixOperator InitialGuess, double tol, int Verbose)
{
   int ncv = 0;
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   CHECK_EQUAL(Psi1.size() % StringOp.size(), 0);
   PackStateComponent Pack(StringOp.Basis1(), Psi1.Basis2(), Psi2.Basis2());
   int n = Pack.size();
   int NumEigen = 1;
   std::vector<std::complex<double>> Initial(n);
   PackMatrixOperator(InitialGuess).pack(InitialGuess, Initial.data());

   std::vector<std::complex<double>> OutVec;
   LinearAlgebra::Vector<std::complex<double>> RightEigen =
      LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(Pack,
                                                         RightMultiplyOperator(Psi1, QShift,
                                                                              StringOp,
                                                                              Psi2, QShift, Psi1.size(), Verbose-1)),
                                          n, NumEigen, Initial.data(), tol, &OutVec, ncv, false, Verbose);

   StateComponent RightVector = Pack.unpack(&(OutVec[0]));

   return std::make_tuple(RightEigen[0], RightVector[0]);
}


std::tuple<std::complex<double>, MatrixOperator>
get_right_transfer_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, ProductMPO const& StringOp, double tol, int Verbose)
{
   int ncv = 0;
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   CHECK_EQUAL(Psi1.size() % StringOp.size(), 0);
   PackStateComponent Pack(StringOp.Basis1(), Psi1.Basis2(), Psi2.Basis2());
   int n = Pack.size();
   int NumEigen = 1;

   std::vector<std::complex<double>> OutVec;
   LinearAlgebra::Vector<std::complex<double>> RightEigen =
      LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(Pack,
                                                         RightMultiplyOperator(Psi1, QShift,
                                                                              StringOp,
                                                                              Psi2, QShift, Psi1.size(), Verbose-1)),
                                          n, NumEigen, nullptr, tol, &OutVec, ncv, false, Verbose);

   StateComponent RightVector = Pack.unpack(&(OutVec[0]));

   return std::make_tuple(RightEigen[0], RightVector[0]);
}

std::tuple<std::complex<double>, MatrixOperator>
get_right_transfer_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, double tol, int Verbose)
{
   return get_right_transfer_eigenvector(Psi1, Psi2, QShift, ProductMPO::make_identity(ExtractLocalBasis(Psi1)), tol, Verbose);
}

std::tuple<std::complex<double>, MatrixOperator>
get_right_transfer_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, MatrixOperator InitialGuess, double tol, int Verbose)
{
   return get_right_transfer_eigenvector(Psi1, Psi2, QShift, ProductMPO::make_identity(ExtractLocalBasis(Psi1), InitialGuess.TransformsAs()), InitialGuess, tol, Verbose);
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

std::tuple<std::complex<double>, MatrixOperator, MatrixOperator>
get_transfer_unit_eigenpair(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift,
                            ProductMPO const& StringOp, double tol, double UnityEpsilon, int Verbose)
{
   std::complex<double> eLeft, eRight;
   MatrixOperator Left, Right;

   std::tie(eLeft, Left) = get_left_transfer_eigenvector(Psi1, Psi2, QShift, StringOp, tol, Verbose);

   // If the norm of the left eigenvalue is less than one, return early.
   if (std::abs(std::abs(eLeft) - 1.0) > UnityEpsilon)
      return std::make_tuple(eLeft, MatrixOperator(), MatrixOperator());

   std::tie(eRight, Right) = get_right_transfer_eigenvector(Psi1, Psi2, QShift, StringOp, tol, Verbose);

   CHECK(norm_frob(eLeft-conj(eRight)) < 1E-10)(eLeft)(eRight);
   // normalize
   Right *= 1.0 / norm_frob(Right);
   Left *= 1.0 / inner_prod(delta_shift(Right, QShift), Left);
   return std::make_tuple(eLeft, Left, Right);
}

std::tuple<std::complex<double>, MatrixOperator, MatrixOperator>
get_transfer_unit_eigenpair(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift,
                            QuantumNumber const& q, double tol, double UnityEpsilon, int Verbose)
{
   return get_transfer_unit_eigenpair(Psi1, Psi2, QShift, ProductMPO::make_identity(ExtractLocalBasis(Psi1), q), tol, UnityEpsilon, Verbose);
}

std::tuple<std::complex<double>, MatrixOperator, MatrixOperator>
get_transfer_unit_eigenpair(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift,
                            double tol, double UnityEpsilon, int Verbose)
{
   return get_transfer_unit_eigenpair(Psi1, Psi2, QShift, QuantumNumbers::QuantumNumber(Psi1.GetSymmetryList()), tol, UnityEpsilon, Verbose);
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
   std::size_t n = PackL.size();

   if (Verbose >= 1)
   {
      std::cerr << "Calculating left eigenvalues\n";
   }

   // There is a problem here matching eigenvalues between the left and the right eigenvalues in the case
   // where the last eigenvalue is one of a complex conjugate pair. In that case, it is not uncommon that the final
   // left eigenvalue is the conjugate pair of the final right eigenvalue, and we get a warning that the eigenvalues
   // don't match.  The fix for this is to calculate n+1 eigenvalues, and throw the last one away at the end.
   std::vector<std::complex<double>> LeftVec;
   std::vector<std::complex<double>>* OutVecL
      = LeftVectors ? &LeftVec : NULL;
   LinearAlgebra::Vector<std::complex<double>>  LeftEigen =
   LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(PackL, LeftMultiplyOperator(Psi, QShift, StringOp, Psi, QShift,
      Psi.size(), Verbose-1)), n+1, NumEigen, nullptr, tol, OutVecL, ncv, Sort, Verbose);

   if (RightVectors)
   {
      if (Verbose >= 1)
      {
         std::cerr << "Calculating right eigenvalues\n";
      }
      std::vector<std::complex<double>> RightVec;
      LinearAlgebra::Vector<std::complex<double>> RightEigen =
      LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(PackR, RightMultiplyOperator(Psi, QShift, StringOp, Psi, QShift,
      Psi.size(), Verbose-1)), n+1, NumEigen, nullptr, tol, &RightVec, ncv, Sort, Verbose);

      // The right vectors are the hermitian conjugate, not the transpose.  So conjugate eigenvalues
      RightEigen = conj(RightEigen);
      *RightVectors = LinearAlgebra::Vector<MatrixOperator>(RightEigen.size());

      // If we have both left & right eigenvectors, then match the corresponding eigenvalues
      if (LeftVectors)
         MatchEigenvectors(n, LeftEigen, *OutVecL, RightEigen, RightVec, tol*10, true);

      // Unpack the eigenvectors into the output array.
      // note that we do not conjugate the eigenvector, since we want this to be a 'column vector',
      // that we will use on the left-hand side of an inner product (eg inner_prod(left, right)).
      for (unsigned i = 0; i < std::min(n,RightEigen.size()); ++i)
      {
         (*RightVectors)[i] = PackR.unpack(&(RightVec[n*i]))[0];
      }
   }

   // eigenvectors
   if (LeftVectors)
   {
      *LeftVectors = LinearAlgebra::Vector<MatrixOperator>(LeftEigen.size());
      for (unsigned i = 0; i < std::min(n,LeftEigen.size()); ++i)
      {
         (*LeftVectors)[i] = PackL.unpack(&((*OutVecL)[n*i]))[0];
      }
   }

   delete LeftVectors;
   delete RightVectors;

   return LeftEigen;
}
