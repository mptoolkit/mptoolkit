// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/infinitewavefunctionright.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "infinitewavefunctionright.h"
#include "infinitewavefunctionleft.h"
#include "tensor/tensor_eigen.h"
#include "mp-algorithms/arnoldi.h"
#include "common/environment.h"

#include "mps/packunpack.h"
#include <fstream>

#include "wavefunction/operator_actions.h"
#include "pheap/pheapstream.h"

// Streaming versions:
// Note: the base class CanonicalWavefunctionBase has a separate version number.
//
// Version 1:
//      CanonicalWavefunctionBase (base class)
//      QuantumNumber              QShift

extern double const ArnoldiTol;
extern double const InverseTol;
extern double const OrthoTol;

PStream::VersionTag
InfiniteWavefunctionRight::VersionT(1);

namespace
{

struct LeftMultiply
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   LeftMultiply(LinearWavefunction const& L_, QuantumNumber const& QShift_)
      : L(L_), QShift(QShift_) {}

   result_type operator()(argument_type const& x) const
   {
      result_type r = delta_shift(x, QShift);
      for (LinearWavefunction::const_iterator I = L.begin(); I != L.end(); ++I)
      {
         r = operator_prod(herm(*I), r, *I);
      }
      return r;
   }

   LinearWavefunction const& L;
   QuantumNumber QShift;
};

struct RightMultiply
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   RightMultiply(LinearWavefunction const& R_, QuantumNumber const& QShift_)
      : R(R_), QShift(QShift_) {}

   result_type operator()(argument_type const& x) const
   {
      result_type r = x;
      LinearWavefunction::const_iterator I = R.end();
      while (I != R.begin())
      {
         --I;
         r = operator_prod(*I, r, herm(*I));
      }
      return delta_shift(r, adjoint(QShift));
   }

   LinearWavefunction const& R;
   QuantumNumber QShift;
};

} // namespace

InfiniteWavefunctionRight::InfiniteWavefunctionRight(MatrixOperator const& Lambda,
                                                     LinearWavefunction const& Psi,
                                                     QuantumNumbers::QuantumNumber const& QShift_)
   : QShift(QShift_)
{
   this->Initialize(Lambda, Psi);
}

void
InfiniteWavefunctionRight::Initialize(MatrixOperator const& Lambda,
                                      LinearWavefunction const& Psi_)
{
   LinearWavefunction Psi = Psi_;
   MatrixOperator M = left_orthogonalize(Lambda, Psi);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   // we need to assemble the MPS in reverse order, so make a temporary container
   std::list<StateComponent> AMat;
   std::list<RealDiagonalOperator> Lam;

   LinearWavefunction::const_iterator I = Psi.end();
   while (I != Psi.begin())
   {
      --I;
      StateComponent A = prod(*I, M);
      M = ExpandBasis1(A);
      SingularValueDecomposition(M, U, D, Vh);
      AMat.push_front(prod(Vh, A));
      Lam.push_front(D);
      M = U*D;
   }
   U = delta_shift(U, adjoint(QShift));
   AMat.back() = prod(AMat.back(), U);
   Lam.push_back(delta_shift(D, adjoint(QShift)));

   for (auto const& A : AMat)
   {
      this->push_back(A);
   }

   for (auto const& L : Lam)
   {
      this->push_back_lambda(L);
   }

   this->setBasis1(D.Basis1());
   this->setBasis2(U.Basis2());


#if !defined(NDEBUG)
   // check
   LinearWavefunction PsiCheck(this->base_begin(), this->base_end());
   MatrixOperator II = MatrixOperator::make_identity(this->Basis2());
   MatrixOperator ICheck = delta_shift(inject_right(II, PsiCheck), adjoint(QShift));
   TRACE(norm_frob(II-ICheck))("Identity - should be epsilon");

   MatrixOperator Rho = this->lambda(0)*this->lambda(0);
   MatrixOperator RhoCheck = inject_left(Rho, PsiCheck);
   RhoCheck = delta_shift(RhoCheck, QShift);
   TRACE(norm_frob(Rho-RhoCheck))("Rho - should be epsilon");
#endif

   this->debug_check_structure();
}

InfiniteWavefunctionRight::InfiniteWavefunctionRight(LinearWavefunction const& Psi,
                                                     QuantumNumbers::QuantumNumber const& QShift_)
   : QShift(QShift_)
{
   LinearWavefunction PsiR = Psi;

   MatrixOperator Guess = MatrixOperator::make_identity(PsiR.Basis2());

   MatrixOperator RightEigen = Guess;

   // get the eigenmatrix.  Do some dodgy explict restarts.
   int Iterations = 20;
   double Tol = ArnoldiTol;
   RightEigen = 0.5 * (RightEigen + adjoint(RightEigen)); // make the eigenvector symmetric
   std::complex<double> EtaR = LinearSolvers::Arnoldi(RightEigen, RightMultiply(PsiR, QShift),
                                                      Iterations, Tol,
                                                      LinearSolvers::LargestAlgebraicReal, false);
   while (Tol < 0)
   {
      std::cerr << "RightEigen: Arnoldi not converged, restarting.  EValue="
                << EtaR << ", Tol=" << Tol << "\n";
      Iterations = 20; Tol = ArnoldiTol;
      RightEigen = 0.5 * (RightEigen + adjoint(RightEigen)); // make the eigenvector symmetric
      EtaR = LinearSolvers::Arnoldi(RightEigen, RightMultiply(PsiR, QShift),
                                    Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   }

   CHECK(EtaR.real() > 0)("Eigenvalue must be positive");

   DEBUG_TRACE(EtaR);

   if (trace(RightEigen).real() < 0)
      RightEigen = -RightEigen;

   DEBUG_CHECK(norm_frob(RightEigen - adjoint(RightEigen)) < 1e-12)
      (norm_frob(RightEigen - adjoint(RightEigen)));

   MatrixOperator D = RightEigen;
   MatrixOperator U = DiagonalizeHermitian(D);
   D = SqrtDiagonal(D, OrthoTol);
   MatrixOperator DInv = InvertDiagonal(D, OrthoTol);

   // RightEigen = triple_prod(U, D*D, herm(U))
   DEBUG_CHECK(norm_frob(RightEigen - triple_prod(herm(U), D*D, U)) < 1e-10)
      (norm_frob(RightEigen - triple_prod(herm(U), D*D, U)));

   MatrixOperator R = adjoint(U)*D;
   MatrixOperator RInv = delta_shift(DInv * U, QShift);

   // Incorporate into the MPS: PsiR -> R^{-1} * PsiR * R
   // We don't need to actually right-orthogonalize everything, that is done in Initialize() anyway
   PsiR.set_back(prod(PsiR.get_back(), R));
   PsiR.set_front(prod(RInv, PsiR.get_front()));

   // Get the left eigenvector, which is the density matrix
   MatrixOperator LeftEigen = Guess;

   // get the eigenmatrix
   Iterations = 20; Tol = ArnoldiTol;
   LeftEigen = 0.5 * (LeftEigen + adjoint(LeftEigen));
   std::complex<double> EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMultiply(PsiR, QShift),
                                                      Iterations, Tol,
                                                      LinearSolvers::LargestAlgebraicReal, false);
   //   DEBUG_TRACE(norm_frob(LeftEigen - adjoint(LeftEigen)));
   while (Tol < 0)
   {
      std::cerr << "LeftEigen: Arnoldi not converged, restarting.  EValue="
                << EtaL << ", Tol=" << Tol << "\n";
      Iterations = 20; Tol = ArnoldiTol;
      LeftEigen = 0.5 * (LeftEigen + adjoint(LeftEigen));
      EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMultiply(PsiR, QShift),
                                    Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   }

   D = LeftEigen;
   U = DiagonalizeHermitian(D);
   D = SqrtDiagonal(D, OrthoTol);

   // normalize
   D *= 1.0 / norm_frob(D);

   DEBUG_CHECK(norm_frob(LeftEigen - triple_prod(herm(U), D*D, U)) < 1e-10)
      (norm_frob(LeftEigen - triple_prod(herm(U), D*D, U)));

   // incorporate U into the MPS

   PsiR.set_front(prod(U, PsiR.get_front()));
   PsiR.set_back(prod(PsiR.get_back(), adjoint(U)));

   // And now we have the right-orthogonalized form and the left-most lambda matrix
   this->Initialize(D, PsiR);
}

InfiniteWavefunctionRight::InfiniteWavefunctionRight(InfiniteWavefunctionLeft const& Psi)
   : InfiniteWavefunctionRight(LinearWavefunction(Psi.base_begin(), Psi.base_end()), Psi.qshift())
{
}

void read_version(PStream::ipstream& in, InfiniteWavefunctionRight& Psi, int Version)
{
   if (Version == 1)
   {
      int BaseVersion = Psi.CanonicalWavefunctionBase::ReadStream(in);
      in >> Psi.QShift;
      CHECK(BaseVersion >= 3);
   }
   else
   {
      PANIC("This program is too old to read this wavefunction, expected Version <= 2")(Version);
   }

   Psi.debug_check_structure();
}

PStream::ipstream&
operator>>(PStream::ipstream& in, InfiniteWavefunctionRight& Psi)
{
   int Version = in.read<int>();
   PStream::VersionSentry Sentry(in, InfiniteWavefunctionRight::VersionT, Version);
   read_version(in, Psi, Version);
   return in;
}

PStream::opstream& operator<<(PStream::opstream& out, InfiniteWavefunctionRight const& Psi)
{
   out << InfiniteWavefunctionRight::VersionT.default_version();

   Psi.CanonicalWavefunctionBase::WriteStream(out);

   out << Psi.QShift;

   return out;
}

std::pair<RealDiagonalOperator, LinearWavefunction>
get_right_canonical(InfiniteWavefunctionRight const& Psi)
{
   return std::make_pair(Psi.lambda(0), LinearWavefunction(Psi.base_begin(), Psi.base_end()));
}

std::tuple<LinearWavefunction, RealDiagonalOperator, MatrixOperator>
get_left_canonical(InfiniteWavefunctionRight const& Psi)
{
   LinearWavefunction Result;
   RealDiagonalOperator D = Psi.lambda_l();
   MatrixOperator Vh = MatrixOperator::make_identity(D.Basis2());
   for (auto const& I : Psi)
   {
      StateComponent A = prod(D*Vh, I);
      std::tie(D,Vh) = OrthogonalizeBasis2(A);
      Result.push_back(A);
   }

   return std::make_tuple(Result, D, Vh);
}

void
InfiniteWavefunctionRight::rotate_left(int Count)
{
   // Rotation is fairly straightforward, we just rotate the vectors around
   if (Count < 0)
   {
      this->rotate_right(-Count);
      return;
   }

   Count = Count % this->size();
   if (Count == 0)
      return;

   // the first Count elements are going to get shifted to the right hand side, so we need to
   // delta_shift them
   for (mps_iterator I = this->begin_(); I != this->begin_()+Count; ++I)
   {
      I->delta_shift(adjoint(this->qshift()));
   }
   // now do the actual rotation
   std::rotate(this->base_begin_(), this->base_begin_()+Count, this->base_end_());

   // for the Lambda matrices, start by removing the double-counted boundary lambda
   this->pop_back_lambda();
   // and delta-shift
   for (lambda_iterator I = this->lambda_begin_(); I != this->lambda_begin_()+Count; ++I)
   {
      I->delta_shift(adjoint(this->qshift()));
   }
   // and rotate
   std::rotate(this->lambda_base_begin_(), this->lambda_base_begin_()+Count,
               this->lambda_base_end_());
   // and put back the boundary lambda
   this->push_back_lambda(delta_shift(this->lambda_l(), adjoint(this->qshift())));

   // set the right and right basis
   this->setBasis1(lambda_l().Basis1());
   this->setBasis2(lambda_r().Basis2());

   this->debug_check_structure();
}

void
InfiniteWavefunctionRight::rotate_right(int Count)
{
   if (Count < 0)
   {
      this->rotate_left(-Count);
      return;
   }

   Count = Count % this->size();
   if (Count == 0)
      return;

   this->rotate_left(this->size() - Count);
}

void
InfiniteWavefunctionRight::check_structure() const
{
   this->CanonicalWavefunctionBase::check_structure();

   CHECK_EQUAL(this->Basis1(), delta_shift(this->Basis2(), this->qshift()));
}

void inplace_conj(InfiniteWavefunctionRight& Psi)
{
   for (InfiniteWavefunctionRight::mps_iterator I = Psi.begin_(); I != Psi.end_(); ++I)
   {
      *I = conj(*I);
   }
}

void inplace_qshift(InfiniteWavefunctionRight& Psi, QuantumNumbers::QuantumNumber const& Shift)
{
   Psi.setBasis1(delta_shift(Psi.Basis1(), Shift));
   Psi.setBasis2(delta_shift(Psi.Basis2(), Shift));

   for (InfiniteWavefunctionRight::mps_iterator I = Psi.begin_(); I != Psi.end_(); ++I)
   {
      *I = delta_shift(*I, Shift);
   }

   for (InfiniteWavefunctionRight::lambda_iterator I = Psi.lambda_begin_(); I != Psi.lambda_end_(); ++I)
   {
      *I = delta_shift(*I, Shift);
   }

   Psi.check_structure();
}


InfiniteWavefunctionRight
reflect(InfiniteWavefunctionLeft const& Psi)
{
   PANIC("not implemented");
}

std::tuple<std::complex<double>, int, StateComponent>
overlap(InfiniteWavefunctionRight const& x, ProductMPO const& StringOp,
        InfiniteWavefunctionRight const& y,
        QuantumNumbers::QuantumNumber const& Sector, int Iter, double Tol, int Verbose)
{
   int Length = statistics::lcm(x.size(), y.size(), StringOp.size());

   LinearWavefunction xPsi = get_right_canonical(x).second;
   LinearWavefunction yPsi = get_right_canonical(y).second;

   ProductMPO Str = StringOp * ProductMPO::make_identity(StringOp.LocalBasis2List(), Sector);

   StateComponent Init = MakeRandomStateComponent(Str.Basis1(), x.Basis1(), y.Basis1());

   int Iterations = Iter;
   int TotalIterations = 0;
   double MyTol = Tol;
   if (Verbose > 1)
   {
      std::cerr << "Starting Arnoldi, Tol=" << MyTol << ", Iterations=" << Iter << '\n';
   }
   std::complex<double> Eta = LinearSolvers::Arnoldi(Init,
                                                     RightMultiplyOperator(xPsi, x.qshift(), Str,
                                                                           yPsi, y.qshift(), Length),
                                                     Iterations,
                                                     MyTol,
                                                     LinearSolvers::LargestMagnitude, false, Verbose);
   TotalIterations += Iterations;
   DEBUG_TRACE(Eta)(Iterations);

   while (MyTol < 0)
   {
      if (Verbose > 0)
         std::cerr << "Restarting Arnoldi, eta=" << Eta << ", Tol=" << -MyTol << '\n';
      Iterations = Iter;
      MyTol = Tol;
      Eta = LinearSolvers::Arnoldi(Init, RightMultiplyOperator(xPsi, x.qshift(), Str,
                                                               yPsi, y.qshift(), Length),
                                   Iterations, MyTol, LinearSolvers::LargestMagnitude, false, Verbose);
      TotalIterations += Iterations;
      DEBUG_TRACE(Eta)(Iterations);
   }
   if (Verbose > 0)
      std::cerr << "Converged.  TotalIterations=" << TotalIterations
                << ", Tol=" << MyTol << '\n';

   return std::make_tuple(Eta, Length, Init);
}

std::pair<std::complex<double>, StateComponent>
overlap(InfiniteWavefunctionRight const& x,  InfiniteWavefunctionRight const& y,
        QuantumNumbers::QuantumNumber const& Sector, int Iter, double Tol, int Verbose)
{
   CHECK_EQUAL(x.size(), y.size());
   std::tuple<std::complex<double>, int, StateComponent> Result =
      overlap(x, ProductMPO::make_identity(ExtractLocalBasis(y)), y, Sector, Iter, Tol, Verbose);
   return std::make_pair(std::get<0>(Result), std::get<2>(Result));
}
