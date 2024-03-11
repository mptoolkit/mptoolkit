// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// wavefunction/infinitewavefunctionright.cpp
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
//
// Version 2 adds:
//      double                     LogAmplitude

extern double const ArnoldiTol;
extern double const InverseTol;
extern double const OrthoTol;

std::string InfiniteWavefunctionRight::Type = "InfiniteWavefunctionRight";

PStream::VersionTag
InfiniteWavefunctionRight::VersionT(2);

InfiniteWavefunctionRight::InfiniteWavefunctionRight(InfiniteWavefunctionLeft const& Psi)
{
   this->setBasis1(Psi.Basis1());
   this->setBasis2(Psi.Basis2());
   QShift = Psi.qshift();
   LogAmplitude = Psi.log_amplitude();

   auto PsiI = Psi.begin();
   auto LambdaI = Psi.lambda_begin();
   ++LambdaI; // skip to the first bond within the unit cell
   while (PsiI != Psi.end())
   {
      StateComponent A = (*PsiI) * (*LambdaI);
      MatrixOperator U;
      RealDiagonalOperator Lambda;
      std::tie(U, Lambda) = OrthogonalizeBasis1(A);
      A = U*A;
      MatrixOperator LambdaP = U*Lambda*herm(U);
      Lambda = ExtractRealDiagonal(LambdaP);
      // if (norm_frob(LambdaP-Lambda) > 1E-14)
      // {
      //    std::cerr << "get_right_canonical: warning: Lambda matrix isn't diagonal, difference = " << norm_frob(LambdaP-Lambda) << '\n';
      // }
      this->push_back(A);
      this->push_back_lambda(Lambda);
      ++PsiI;
      ++LambdaI;
   }
   this->push_back_lambda(delta_shift(this->lambda(0), adjoint(QShift)));
   this->check_structure();
}

void read_version(PStream::ipstream& in, InfiniteWavefunctionRight& Psi, int Version)
{
   int BaseVersion = Psi.CanonicalWavefunctionBase::ReadStream(in);
   CHECK(BaseVersion >= 3);
   in >> Psi.QShift;
   if (Version == 1)
   {
      Psi.LogAmplitude = 0;
   }
   else if (Version == 2)
   {
      in >> Psi.LogAmplitude;
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
   out << Psi.LogAmplitude;

   return out;
}

std::pair<RealDiagonalOperator, LinearWavefunction>
get_right_canonical(InfiniteWavefunctionRight const& Psi)
{
   return std::make_pair(Psi.lambda(0), LinearWavefunction(Psi.base_begin(), Psi.base_end()));
}

std::tuple<LinearWavefunction, RealDiagonalOperator>
get_left_canonical(InfiniteWavefunctionRight const& Psi)
{
   LinearWavefunction Result;
   RealDiagonalOperator Lambda;
   auto PsiI = Psi.begin();
   auto LambdaI = Psi.lambda_begin();
   while (PsiI != Psi.end())
   {
      StateComponent A = (*LambdaI) * (*PsiI);
      MatrixOperator Vh;
      std::tie(Lambda, Vh) = OrthogonalizeBasis2(A);
      A = A*Vh;
      MatrixOperator LambdaP = herm(Vh)*(Lambda*Vh);
      Lambda = ExtractRealDiagonal(LambdaP);
      // if (norm_frob(LambdaP-Lambda) > 1E-14)
      // {
      //    std::cerr << "get_right_canonical: warning: Lambda matrix isn't diagonal, difference = " << norm_frob(LambdaP-Lambda) << '\n';
      // }
      Result.push_back(A);
      ++PsiI;
      ++LambdaI;
   }
   return std::make_tuple(Result, Lambda);
}

void InfiniteWavefunctionRight::scale(std::complex<double> x)
{
   double Amplitude = std::abs(x);
   x = x / Amplitude;  // x is now the phase
   LogAmplitude += std::log(Amplitude);
   // Because the components have value semantics we can't directly load() the pvalue_ptr
   // and then mutate it, because pvalue_handle's are immutable - we would end up referring
   // do a different object.
   pvalue_ptr<StateComponent> s = this->base_begin_()->load();
   *s.mutate() *= x;
   *this->base_begin_() = s;
}

void InfiniteWavefunctionRight::scale_log(std::complex<double> x)
{
   LogAmplitude += x.real();
   // Because the components have value semantics we can't directly load() the pvalue_ptr
   // and then mutate it, because pvalue_handle's are immutable - we would end up referring
   // do a different object.
   if (x.imag() != 0.0)
   {
      pvalue_ptr<StateComponent> s = this->base_begin_()->load();
      *s.mutate() *= std::exp(std::complex<double>(0,x.imag()));
      *this->base_begin_() = s;
   }
}

void InfiniteWavefunctionRight::normalize()
{
   LogAmplitude = 0;
}

InfiniteWavefunctionRight& operator*=(InfiniteWavefunctionRight& psi, std::complex<double> x)
{
   psi.scale(x);
   return psi;
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
InfiniteWavefunctionRight::SetDefaultAttributes(AttributeList& A) const
{
   A["WavefunctionType"] = "InfiniteRightCanonical";
   A["UnitCellSize"] = this->size();
   A["QShift"] = this->qshift();
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

void inplace_reflect(InfiniteWavefunctionRight& Psi)
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

   ProductMPO Str = StringOp * ProductMPO::make_identity(StringOp.LocalBasis1List(), Sector);

   StateComponent Init = MakeRandomStateComponent(Str.Basis2(), x.Basis2(), y.Basis2());

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
