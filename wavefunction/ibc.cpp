// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/ibc.cpp
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

#include "ibc.h"
#include "tensor/tensor_eigen.h"

//
// WavefunctionSectionLeft
//

// Streaming version:
//
// Version 1:
// Base class CanonicalWavefunctionBase

PStream::VersionTag
WavefunctionSectionLeft::VersionT(1);

WavefunctionSectionLeft::WavefunctionSectionLeft()
{
}

std::string const WavefunctionSectionLeft::Type = "WavefunctionSectionLeft";

WavefunctionSectionLeft::WavefunctionSectionLeft(InfiniteWavefunctionLeft const& Psi)
   : CanonicalWavefunctionBase(Psi)
{
   LeftU_ = MatrixOperator::make_identity(this->Basis1());
   RightU_ = MatrixOperator::make_identity(this->Basis2());
   this->check_structure();
}

WavefunctionSectionLeft::WavefunctionSectionLeft(MatrixOperator const& C)
   : CanonicalWavefunctionBase(C.Basis1(), C.Basis2())
{
   RealDiagonalOperator Lambda;
   SingularValueDecompositionFull(C, LeftU_, Lambda, RightU_);
   this->push_back_lambda(Lambda);
}

PStream::opstream& operator<<(PStream::opstream& out, WavefunctionSectionLeft const& Psi)
{
   out << WavefunctionSectionLeft::VersionT.default_version();
   Psi.CanonicalWavefunctionBase::WriteStream(out);
   out << Psi.LeftU_;
   out << Psi.RightU_;
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, WavefunctionSectionLeft& Psi)
{
   PStream::VersionSentry Sentry(in, WavefunctionSectionLeft::VersionT, in.read<int>());
   if (Sentry.version() != 1)
   {
      PANIC("This program is too old to read this wavefunction, expected version = 1")(Sentry.version());
   }
   Psi.CanonicalWavefunctionBase::ReadStream(in);
   in >> Psi.LeftU_;
   in >> Psi.RightU_;
   return in;
}

void
inplace_reflect(WavefunctionSectionLeft& Psi)
{
   PANIC("Reflect() not yet implemented for WavefunctionSectionLeft");
}

void
inplace_conj(WavefunctionSectionLeft& Psi)
{
   for (WavefunctionSectionLeft::mps_iterator I = Psi.begin_(); I != Psi.end_(); ++I)
   {
      *I = conj(*I);
   }
}

WavefunctionSectionLeft
WavefunctionSectionLeft::ConstructFromLeftOrthogonal(LinearWavefunction const& Psi,
                                                     MatrixOperator const& Lambda,
                                                     int Verbose)
{
   return ConstructFromLeftOrthogonal(std::move(LinearWavefunction(Psi)), Lambda, Verbose);
}

// version of ConstructFromLeftOrthogonal with move semantics on Psi
WavefunctionSectionLeft
WavefunctionSectionLeft::ConstructFromLeftOrthogonal(LinearWavefunction&& Psi,
                                                     MatrixOperator const& Lambda,
                                                     int Verbose)
{
   WavefunctionSectionLeft Result;
   if (Verbose > 0)
   {
      std::cout << "Constructing canonical wavefunction..." << std::endl;
      std::cout << "Constructing right ortho matrices..." << std::endl;
   }

   MatrixOperator M = right_orthogonalize(Psi, Lambda, Verbose-1);

   MatrixOperator U;
   RealDiagonalOperator D;
   MatrixOperator Vh;
   SingularValueDecomposition(M, U, D, Vh);

   Result.LeftU_ = U;
   Result.push_back_lambda(D);
   Result.setBasis1(D.Basis1());

   M = D*Vh;

   if (Verbose > 0)
      std::cout << "Constructing left ortho matrices..." << std::endl;


   int n = 0;
   while (!Psi.empty())
   {
      if (Verbose > 1)
         std::cout << "orthogonalizing site " << n << std::endl;
      StateComponent A = prod(M, Psi.get_front());
      Psi.pop_front();
      std::tie(D, Vh) = OrthogonalizeBasis2(A);
      Result.push_back(A);
      Result.push_back_lambda(D);
      M = D*Vh;
      ++n;
   }
   Result.setBasis2(D.Basis2());
   Result.RightU_ = Vh;

   if (Verbose > 0)
      std::cout << "Finished constructing canonical wavefunction." << std::endl;

   return Result;
}

std::pair<LinearWavefunction, MatrixOperator>
get_left_canonical(WavefunctionSectionLeft const& Psi)
{
   LinearWavefunction PsiLinear(Psi.base_begin(), Psi.base_end());
   MatrixOperator Lambda = Psi.lambda_r();
   // incorporate the U matrices
   Lambda = Lambda*Psi.RightU();
   PsiLinear.set_front(prod(Psi.LeftU(), PsiLinear.get_front()));
   return std::make_pair(PsiLinear, Lambda);
}

void
WavefunctionSectionLeft::check_structure() const
{
   this->CanonicalWavefunctionBase::check_structure();
   CHECK_EQUAL(this->Basis1(), LeftU_.Basis2());
   CHECK_EQUAL(this->Basis2(), RightU_.Basis1());
}

//
// IBCWavefunction
//

// Streaming versions:
//
// Version 1:
// int WindowLeftSites
// int WindowRightSites
// int WindowOffset
// InfiniteWavefunctionLeft Left
// WavefunctionSectionLeft Window
// InfiniteWavefunctionRight Right

PStream::VersionTag
IBCWavefunction::VersionT(1);

std::string const IBCWavefunction::Type = "IBCWavefunction";

IBCWavefunction::IBCWavefunction()
   : WindowLeftSites(0), WindowRightSites(0), WindowOffset(0)
{
}

IBCWavefunction::IBCWavefunction(InfiniteWavefunctionLeft const& Left_,
                                 WavefunctionSectionLeft const& Window_,
                                 InfiniteWavefunctionRight const& Right_,
                                 int Offset)
   : WindowLeftSites(0), WindowRightSites(0), WindowOffset(Offset),
     Left(Left_), Window(Window_), Right(Right_)
{
}

IBCWavefunction::IBCWavefunction(InfiniteWavefunctionLeft const& Left_,
                                 WavefunctionSectionLeft const& Window_,
                                 InfiniteWavefunctionRight const& Right_,
                                 int Offset,
                                 int WindowLeft,
                                 int WindowRight)
   : WindowLeftSites(WindowLeft), WindowRightSites(WindowRight), WindowOffset(Offset),
     Left(Left_), Window(Window_), Right(Right_)
{
}

void IBCWavefunction::check_structure() const
{
   Left.check_structure();
   Window.check_structure();
   Right.check_structure();
   if (!Left.empty())
   {
      CHECK_EQUAL(Left.Basis2(), Window.Basis1());
   }
   if (!Right.empty())
   {
      CHECK_EQUAL(Window.Basis2(), Right.Basis1());
   }
}

void IBCWavefunction::debug_check_structure() const
{
   Left.debug_check_structure();
   Window.debug_check_structure();
   Right.debug_check_structure();
   if (!Left.empty())
   {
      DEBUG_CHECK_EQUAL(Left.Basis2(), Window.Basis1());
   }
   if (!Right.empty())
   {
      DEBUG_CHECK_EQUAL(Window.Basis2(), Right.Basis1());
   }
}

PStream::opstream& operator<<(PStream::opstream& out, IBCWavefunction const& Psi)
{
   out << IBCWavefunction::VersionT.default_version();
   out << Psi.WindowLeftSites;
   out << Psi.WindowRightSites;
   out << Psi.WindowOffset;
   out << Psi.Left;
   out << Psi.Window;
   out << Psi.Right;
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, IBCWavefunction& Psi)
{
   int Version = in.read<int>();
   PStream::VersionSentry Sentry(in, IBCWavefunction::VersionT, Version);

   if (Version != 1)
   {
      PANIC("This program is too old to read this wavefunction, expected version = 1")(Version);
   }

   in >> Psi.WindowLeftSites;
   in >> Psi.WindowRightSites;
   in >> Psi.WindowOffset;
   in >> Psi.Left;
   in >> Psi.Window;
   in >> Psi.Right;

   return in;
}

void
inplace_reflect(IBCWavefunction& Psi)
{
   InfiniteWavefunctionRight Temp = reflect(Psi.Left);
   Psi.Left = reflect(Psi.Right);
   Psi.Right = Temp;
   inplace_reflect(Psi.Window);
}

void
inplace_conj(IBCWavefunction& Psi)
{
   inplace_conj(Psi.Left);
   inplace_conj(Psi.Window);
   inplace_conj(Psi.Right);
}

void
IBCWavefunction::SetDefaultAttributes(AttributeList& A) const
{
   A["WavefunctionType"] = "IBC";
   A["WindowSize"] = this->window_size();
   A["WindowOffset"] = this->window_offset();
   A["LeftUnitCellSize"] = Left.size();
   A["RightUnitCellSize"] = Left.size();
}
