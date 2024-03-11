// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// wavefunction/finitewavefunctionleft.cpp
//
// Copyright (C) 2017 Ian McCulloch <ian@qusim.net>
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

#include "finitewavefunctionleft.h"
#include "tensor/tensor_eigen.h"

std::string const FiniteWavefunctionLeft::Type = "FiniteWavefunctionLeft";

FiniteWavefunctionLeft
FiniteWavefunctionLeft::ConstructFromRightOrthogonal(LinearWavefunction Psi,
						     std::complex<double> a,
						     int Verbose)
{
   CHECK(Psi.Basis2().is_vacuum());
   CHECK(Psi.Basis1().size() == 1);

   if (Verbose > 0)
      std::cout << "Constructing left ortho matrices..." << std::endl;

   FiniteWavefunctionLeft Result;
   RealDiagonalOperator D(Psi.Basis1());
   double Norm = norm_frob(a);
   D(0,0) = LinearAlgebra::DiagonalMatrix<double>(1,1, Norm);

   Result.push_back_lambda(D);
   Result.setBasis1(D.Basis1());

   MatrixOperator M = D;
   MatrixOperator U, Vh;
   int n = 0;
   for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I, ++n)
   {
      if (Verbose > 1)
	 std::cout << "orthogonalizing site " << n << std::endl;
      StateComponent A = prod(M, *I);
      M = ExpandBasis2(A);
      SingularValueDecomposition(M, U, D, Vh);
      Result.push_back(prod(A, U));
      Result.push_back_lambda(D);
      M = D*Vh;
      ++n;
   }
   // The final Vh is a 1x1 unitary matrix.  To preserve
   // the phase of the wavefunction we need to incorporate it
   // back into the wavefunction.  Also incorporate the phase part of a
   Result.set_back(prod(Result.get_back(), Vh * (a / Norm)));
   Result.setBasis2(D.Basis2());

   if (Verbose > 0)
      std::cout << "Finished constructing canonical wavefunction." << std::endl;

   return Result;
}

FiniteWavefunctionLeft
FiniteWavefunctionLeft::Construct(LinearWavefunction Psi,
				  int Verbose)
{
   CHECK(Psi.Basis2().is_vacuum());
   CHECK(Psi.Basis1().size() == 1);

   if (Verbose > 0)
      std::cout << "Constructing canonical wavefunction..." << std::endl;

   MatrixOperator M = right_orthogonalize(Psi, Verbose-1);

   return ConstructFromRightOrthogonal(Psi, trace(M), Verbose);
}

void
FiniteWavefunctionLeft::SetDefaultAttributes(AttributeList& A) const
{
   A["WavefunctionType"] = "Finite";
   A["Size"] = this->size();
   A["TransformsAs"] = this->TransformsAs();
}


PStream::VersionTag
FiniteWavefunctionLeft::VersionT(1);

PStream::ipstream&
operator>>(PStream::ipstream& in, FiniteWavefunctionLeft& Psi)
{
   int Version = in.read<int>();
   if (Version > 1)
   {
      PANIC("FiniteWavefunctionLeft version is more recent than this program.");
   }

   Psi.CanonicalWavefunctionBase::ReadStream(in);

   return in;
}

PStream::opstream&
operator<<(PStream::opstream& out, FiniteWavefunctionLeft const& Psi)
{
   out << FiniteWavefunctionLeft::VersionT.default_version();
   Psi.CanonicalWavefunctionBase::WriteStream(out);

   return out;
}

std::complex<double>
overlap(FiniteWavefunctionLeft const& Psi1, FiniteWavefunctionLeft const& Psi2)
{
   if (Psi1.empty() || Psi2.empty())
      return 0.0;
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   if (Psi1.TransformsAs() != Psi2.TransformsAs())
      return 0.0;

   MatrixOperator E = scalar_prod(herm(Psi1.lambda_l()), Psi2.lambda_l());
   FiniteWavefunctionLeft::const_mps_iterator I1 = Psi1.begin(), I2 = Psi2.begin();
   while (I1 != Psi1.end())
   {
      E = operator_prod(herm(*I1), E, *I2);
      ++I1;
      ++I2;
   }
   return trace(E);
}

std::complex<double>
overlap_conj(FiniteWavefunctionLeft const& Psi1, FiniteWavefunctionLeft const& Psi2)
{
   if (Psi1.empty() || Psi2.empty())
      return 0.0;
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   if (Psi1.TransformsAs() != Psi2.TransformsAs())
      return 0.0;

   MatrixOperator E = scalar_prod(herm(Psi1.lambda_l()), Psi2.lambda_l());
   FiniteWavefunctionLeft::const_mps_iterator I1 = Psi1.begin(), I2 = Psi2.begin();
   while (I1 != Psi1.end())
   {
      E = operator_prod(herm(*I1), E, conj(*I2));
      ++I1;
      ++I2;
   }
   return trace(E);
}

std::complex<double>
expectation(FiniteWavefunctionLeft const& Psi1,
            BasicFiniteMPO const& M,
            FiniteWavefunctionLeft const& Psi2)
{
   if (Psi1.empty() || Psi2.empty())
      return 0.0;
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   CHECK(M.is_irreducible());
   CHECK_EQUAL(Psi1.size(), M.size());
   if (Psi1.TransformsAs() != Psi2.TransformsAs())
      return 0.0;

   StateComponent E(M.Basis1(), Psi1.Basis1(), Psi2.Basis1());
   E[0] = scalar_prod(herm(Psi1.lambda_l()), Psi2.lambda_l());
   FiniteWavefunctionLeft::const_mps_iterator I1 = Psi1.begin(), I2 = Psi2.begin();
   BasicFiniteMPO::const_iterator I = M.begin();
   while (I1 != Psi1.end())
   {
      E = contract_from_left(*I, herm(*I1), E, *I2);
      ++I1;
      ++I2;
      ++I;
   }
   return trace(E[0]);
}

double norm_2(FiniteWavefunctionLeft const& Psi)
{
   return norm_frob(Psi.lambda_l());
}

double norm_2_sq(FiniteWavefunctionLeft const& Psi)
{
   return norm_frob_sq(Psi.lambda_l());
}

FiniteWavefunctionLeft operator*(double a, FiniteWavefunctionLeft x)
{
   x *= a;
   return x;
}

FiniteWavefunctionLeft operator*(FiniteWavefunctionLeft x, double a)
{
   x *= a;
   return x;
}

FiniteWavefunctionLeft operator*(std::complex<double> a, FiniteWavefunctionLeft x)
{
   x *= a;
   return x;
}

FiniteWavefunctionLeft operator*(FiniteWavefunctionLeft x, std::complex<double> a)
{
   x *= a;
   return x;
}

FiniteWavefunctionLeft& operator*=(FiniteWavefunctionLeft& psi, double a)
{
   for (auto I = psi.lambda_begin_(); I != psi.lambda_end_(); ++I)
   {
      *I *= a;
   }
   return psi;
}

FiniteWavefunctionLeft& operator*=(FiniteWavefunctionLeft& psi, std::complex<double> a)
{
   double aNorm = norm_frob(a);
   psi *= aNorm;
   if (aNorm != 0.0)
      psi.set(0, psi[0] * (a / aNorm));
   return psi;
}

FiniteWavefunctionLeft conj(FiniteWavefunctionLeft Psi)
{
   inplace_conj(Psi);
   return Psi;
}

void inplace_conj(FiniteWavefunctionLeft& Psi)
{
   for (auto I = Psi.begin_(); I != Psi.end_(); ++I)
   {
      *I = conj(*I);
   }
}

FiniteWavefunctionLeft reflect(FiniteWavefunctionLeft Psi)
{
   inplace_reflect(Psi);
   return Psi;
}

void inplace_reflect(FiniteWavefunctionLeft& Psi)
{
   PANIC("inplace_reflect(FiniteWavefunctionLeft) is not yet implemented!");
}
