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

// The tolerance for the left/right boundary overlaps for the overlap of two general IBCs.
double const OverlapTol = 1e-12;

// The tolerance of the trace of the left/right boundary eigenvectors for
// fixing the phase when calculating the overlap of two general IBCs.
double const TraceTol = 1e-8;

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

   Psi.LeftU_ = conj(Psi.LeftU_);
   Psi.RightU_ = conj(Psi.RightU_);
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

// This class provides an "iterator" which runs over the A-matrices in an IBC
// wavefunction, where the left and right boundaries are composed of left- and
// right-orthogonal matrices respectively, the window is formed of
// left-orthogonal matrices. The orthogonality centre is incorporated into the
// first site in the right semi-infinite boundary.
class ConstIBCIterator
{
   public:
      ConstIBCIterator(IBCWavefunction const& Psi_, int Index_)
         : Psi(Psi_), Index(Index_)
      {
         PsiLeft = Psi.Left;
         PsiRight = Psi.Right;

         WindowLeftIndex = Psi.window_offset();
         WindowRightIndex = Psi.window_size() + Psi.window_offset() - 1;

         if (Index < WindowLeftIndex)
         {
            int IndexDiff = WindowLeftIndex - Index;

            for (int i = 0; i < (Psi.WindowLeftSites + IndexDiff) / PsiLeft.size(); ++i)
               inplace_qshift(PsiLeft, PsiLeft.qshift());

            C = PsiLeft.end();
            for (int i = 0; i < (Psi.WindowLeftSites + IndexDiff) % PsiLeft.size(); ++i)
               --C;
            
            if (C == PsiLeft.end())
            {
               inplace_qshift(PsiLeft, adjoint(PsiLeft.qshift()));
               C = PsiLeft.begin();
            }
         }
         else if (Index > WindowRightIndex)
         {
            int IndexDiff = Index - WindowRightIndex - 1;

            for (int i = 0; i < (Psi.WindowRightSites + IndexDiff) / PsiRight.size(); ++i)
               inplace_qshift(PsiRight, adjoint(PsiRight.qshift()));

            C = PsiRight.begin();
            for (int i = 0; i < (Psi.WindowRightSites + IndexDiff) % PsiRight.size(); ++i)
               ++C;
         }
         else
         {
            int IndexDiff = Index - WindowLeftIndex;

            C = Psi.Window.begin();
            for (int i = 0; i < IndexDiff; ++i)
               ++C;
         }
      }

      StateComponent operator*() const
      {
         StateComponent Result = *C;

         if (Index == WindowRightIndex + 1)
            Result = Psi.Window.lambda_r() * Psi.Window.RightU() * Result;

         if (Index == WindowLeftIndex)
            Result = Psi.Window.LeftU() * Result;

         return Result;
      }

      ConstIBCIterator& operator++()
      {
         ++Index;
         if (Index < WindowLeftIndex)
         {
            ++C;
            if (C == PsiLeft.end())
            {
               inplace_qshift(PsiLeft, adjoint(PsiLeft.qshift()));
               C = PsiLeft.begin();
            }
         }
         else if (Index == WindowRightIndex + 1)
         {
            C = PsiRight.begin();
            for (int i = 0; i < Psi.WindowRightSites; ++i)
               ++C;
         }
         else if (Index == WindowLeftIndex)
            C = Psi.Window.begin();
         else if (Index <= WindowRightIndex)
            ++C;
         else if (Index > WindowRightIndex + 1)
         {
            ++C;
            if (C == PsiRight.end())
            {
               inplace_qshift(PsiRight, adjoint(PsiRight.qshift()));
               C = PsiRight.begin();
            }
         }

         return *this;
      }

   private:
      IBCWavefunction Psi;
      InfiniteWavefunctionLeft PsiLeft;
      InfiniteWavefunctionRight PsiRight;
      CanonicalWavefunctionBase::const_mps_iterator C;
      int Index;
      int WindowLeftIndex;
      int WindowRightIndex;
};

std::complex<double>
expectation(IBCWavefunction const& Psi, UnitCellMPO Op, int Verbose)
{
   // We choose IndexLeft/IndexRight such that it is the first/last site in the
   // operator unit cell, in order for Op.ExtendToCover to work correctly.
   int IndexLeft = std::min(Psi.window_offset() - ((Psi.Left.size() - Psi.WindowLeftSites) % Op.unit_cell_size()),
                            Op.offset());
   int IndexRight = std::max(Psi.window_size() + Psi.window_offset() + ((Psi.Right.size() - Psi.WindowRightSites - 1) % Op.unit_cell_size()),
                             Op.size() + Op.offset() - 1);

   if (Verbose > 0)
      std::cout << "Calculating IBC expectation value over sites " << IndexLeft << " to " << IndexRight
                << " (" << IndexRight - IndexLeft + 1 << " sites total)" << std::endl;

   Op.ExtendToCover(IndexRight - IndexLeft + 1, IndexLeft);

   BasicFiniteMPO M = Op.MPO();

   ConstIBCIterator C = ConstIBCIterator(Psi, IndexLeft);

   MatrixOperator I = MatrixOperator::make_identity((*C).Basis1());
   StateComponent E(M.Basis1(), I.Basis1(), I.Basis2());
   E.front() = I;

   auto W = M.begin();

   for (int i = IndexLeft; i <= IndexRight; ++i)
   {
      if (Verbose > 2)
         std::cout << "Site " << i << std::endl;

      E = contract_from_left(*W, herm(*C), E, *C);
      ++C, ++W;
   }

   return trace(E[0]);
}

std::complex<double>
overlap_simple(IBCWavefunction const& Psi1, IBCWavefunction const& Psi2, int Verbose)
{
   int IndexLeft = std::min(Psi1.window_offset(), Psi2.window_offset());
   int IndexRight = std::max(Psi1.window_size() + Psi1.window_offset(),
                             Psi2.window_size() + Psi2.window_offset());

   if (Verbose > 0)
      std::cout << "Calculating IBC expectation value over sites " << IndexLeft << " to " << IndexRight
                << " (" << IndexRight - IndexLeft + 1 << " sites total)" << std::endl;

   ConstIBCIterator C1 = ConstIBCIterator(Psi1, IndexLeft);
   ConstIBCIterator C2 = ConstIBCIterator(Psi2, IndexLeft);

   CHECK((*C1).Basis1() == (*C2).Basis1());

   MatrixOperator I = MatrixOperator::make_identity((*C1).Basis1());
   MatrixOperator E = scalar_prod(herm(I), I);

   for (int i = IndexLeft; i <= IndexRight; ++i)
   {
      if (Verbose > 2)
         std::cout << "Site " << i << std::endl;

      E = operator_prod(herm(*C1), E, *C2);
      ++C1, ++C2;
   }

   return trace(E);
}

std::complex<double>
overlap(IBCWavefunction const& Psi1, IBCWavefunction const& Psi2, int Verbose)
{
   CHECK_EQUAL(Psi1.Left.size(), Psi2.Left.size());
   CHECK_EQUAL(Psi1.Right.size(), Psi2.Right.size());

   int IndexLeft1 = Psi1.window_offset() - ((Psi1.Left.size() - Psi1.WindowLeftSites) % Psi1.Left.size());
   int IndexLeft2 = Psi2.window_offset() - ((Psi2.Left.size() - Psi2.WindowLeftSites) % Psi2.Left.size());
   int IndexLeft = std::min(IndexLeft1, IndexLeft2);

   int IndexRight1 = Psi1.window_size() + Psi1.window_offset() + ((Psi1.Right.size() - Psi1.WindowRightSites - 1) % Psi1.Right.size());
   int IndexRight2 = Psi2.window_size() + Psi2.window_offset() + ((Psi2.Right.size() - Psi2.WindowRightSites - 1) % Psi2.Right.size());
   int IndexRight = std::max(IndexRight1, IndexRight2);

   // Before calculating the overlap, we must find the left/right eigenvectors
   // of the left/right boundary transfer matrices.

   // Get the identity quantum number.
   QuantumNumber QI(Psi1.GetSymmetryList());

   // Ensure that the semi-infinite boundaries have the same quantum numbers.
   QuantumNumber QShiftLeft1 = QI;
   for (int i = 0; i < (IndexLeft1 - IndexLeft) / Psi1.Left.size(); ++i)
      QShiftLeft1 = delta_shift(QShiftLeft1, Psi1.Left.qshift());

   QuantumNumber QShiftLeft2 = QI;
   for (int i = 0; i < (IndexLeft2 - IndexLeft) / Psi2.Left.size(); ++i)
      QShiftLeft2 = delta_shift(QShiftLeft2, Psi2.Left.qshift());

   InfiniteWavefunctionLeft Psi1Left = Psi1.Left;
   inplace_qshift(Psi1Left, QShiftLeft1);

   InfiniteWavefunctionLeft Psi2Left = Psi2.Left;
   inplace_qshift(Psi2Left, QShiftLeft2);

   // Calculate the left eigenvector of the left semi-infinite boundary.
   std::complex<double> OverlapL;
   StateComponent ESC;
   std::tie(OverlapL, ESC) = overlap(Psi1Left, Psi2Left, QI);

   // Check that the eigenvalue has magnitude 1.
   //TRACE(OverlapL);
   if (std::abs(std::abs(OverlapL) - 1.0) > OverlapTol)
      WARNING("The overlap of the left boundaries is below threshold.")(OverlapL);

   MatrixOperator E = delta_shift(ESC.front(), adjoint(Psi1.Left.qshift()));

   QuantumNumber QShiftRight1 = QI;
   for (int i = 0; i < (IndexRight - IndexRight1) / Psi1.Right.size(); ++i)
      QShiftRight1 = delta_shift(QShiftRight1, adjoint(Psi1.Right.qshift()));

   QuantumNumber QShiftRight2 = QI;
   for (int i = 0; i < (IndexRight - IndexRight2) / Psi2.Right.size(); ++i)
      QShiftRight2 = delta_shift(QShiftRight2, adjoint(Psi2.Right.qshift()));

   InfiniteWavefunctionRight Psi1Right = Psi1.Right;
   inplace_qshift(Psi1Right, QShiftRight1);

   InfiniteWavefunctionRight Psi2Right = Psi2.Right;
   inplace_qshift(Psi2Right, QShiftRight2);

   // Calculate the right eigenvector of the right semi-infinite boundary.
   std::complex<double> OverlapR;
   StateComponent FSC;
   std::tie(OverlapR, FSC) = overlap(Psi1Right, Psi2Right, QI);

   // Check that the eigenvalue has magnitude 1.
   //TRACE(OverlapR);
   if (std::abs(std::abs(OverlapR) - 1.0) > OverlapTol)
      WARNING("The overlap of the right boundaries is below threshold.")(OverlapR);

   MatrixOperator F = FSC.front();

   // Rescale E and F by the bond dimension.
   E *= std::sqrt(std::min(E.Basis1().total_degree(), E.Basis2().total_degree()));
   F *= std::sqrt(std::min(F.Basis1().total_degree(), F.Basis2().total_degree()));

   if (E.Basis1() == E.Basis2() && F.Basis1() == F.Basis2())
   {
      //E *= std::sqrt(E.Basis1().total_degree());
      //F *= std::sqrt(F.Basis1().total_degree());

      // Remove spurious phase from E and F by setting the phase of the trace
      // to be zero (but only if the trace is nonzero).

      std::complex<double> ETrace = trace(E);

      if (std::abs(ETrace) > TraceTol)
         E *= std::conj(ETrace) / std::abs(ETrace);
      else
         WARNING("The trace of E is below threshold, so the overlap will have a spurious phase contribution.")(ETrace);

      std::complex<double> FTrace = trace(F);

      if (std::abs(FTrace) > TraceTol)
         F *= std::conj(FTrace) / std::abs(FTrace);
      else
         WARNING("The trace of F is below threshold, so the overlap will have a spurious phase contribution.")(FTrace);
   }
   else
      WARNING("Psi1 and Psi2 have different boundary bases, so the overlap will have a spurious phase contribution.");

   return overlap(Psi1, Psi2, E, F, Verbose);
}

std::complex<double>
overlap(IBCWavefunction const& Psi1, IBCWavefunction const& Psi2,
        MatrixOperator const& E_, MatrixOperator const& F_, int Verbose)
{
   CHECK_EQUAL(Psi1.Left.size(), Psi2.Left.size());
   CHECK_EQUAL(Psi1.Right.size(), Psi2.Right.size());

   int IndexLeft1 = Psi1.window_offset() - ((Psi1.Left.size() - Psi1.WindowLeftSites) % Psi1.Left.size());
   int IndexLeft2 = Psi2.window_offset() - ((Psi2.Left.size() - Psi2.WindowLeftSites) % Psi2.Left.size());
   int IndexLeft = std::min(IndexLeft1, IndexLeft2);

   int IndexRight1 = Psi1.window_size() + Psi1.window_offset() + ((Psi1.Right.size() - Psi1.WindowRightSites - 1) % Psi1.Right.size());
   int IndexRight2 = Psi2.window_size() + Psi2.window_offset() + ((Psi2.Right.size() - Psi2.WindowRightSites - 1) % Psi2.Right.size());
   int IndexRight = std::max(IndexRight1, IndexRight2);

   if (Verbose > 0)
      std::cout << "Calculating IBC overlap over sites " << IndexLeft << " to " << IndexRight
                << " (" << IndexRight - IndexLeft + 1 << " sites total)" << std::endl;

   MatrixOperator E = E_;
   MatrixOperator F = F_;

   // Calculate the overlap.
   ConstIBCIterator C1 = ConstIBCIterator(Psi1, IndexLeft);
   ConstIBCIterator C2 = ConstIBCIterator(Psi2, IndexLeft);

   for (int i = IndexLeft; i <= IndexRight; ++i)
   {
      if (Verbose > 2)
         std::cout << "Site " << i << std::endl;

      E = operator_prod(herm(*C1), E, *C2);
      ++C1, ++C2;
   }

   return inner_prod(F, E);
}
