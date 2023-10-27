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
#include "mpwavefunction.h"
#include "tensor/tensor_eigen.h"
#include "mp-algorithms/transfer.h"

// The tolerance of the trace of the left/right boundary eigenvectors for
// fixing the phase when calculating the overlap of two general IBCs.
double const TraceTol = 1e-8;

//
// WavefunctionSectionLeft
//

// Streaming version:
//
// Version 1:
// Base class CanonicalWavefunctionBase
// LeftU
// RightU

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
//
// Version 2:
// int WindowLeftSites
// int WindowRightSites
// int WindowOffset
// string WavefunctionLeftFile
// string WavefunctionRightFile
// InfiniteWavefunctionLeft Left (only if WavefunctionLeftFile is empty)
// WavefunctionSectionLeft Window
// InfiniteWavefunctionRight Right (only if WavefunctionRightFile is empty)
//
// Version 3:
// int WindowLeftSites
// int WindowRightSites
// int WindowOffset
// string WavefunctionLeftFile
// string WavefunctionRightFile
// InfiniteWavefunctionLeft Left (only if WavefunctionLeftFile is empty)
// WavefunctionSectionLeft Window
// InfiniteWavefunctionRight Right (only if WavefunctionRightFile is empty)
// QuantumNumber LeftQShift
// QuantumNumber RightQShift

PStream::VersionTag
IBCWavefunction::VersionT(3);

std::string const IBCWavefunction::Type = "IBCWavefunction";

IBCWavefunction::IBCWavefunction()
   : WindowLeftSites(0), WindowRightSites(0), WindowOffset(0)
{
}

IBCWavefunction::IBCWavefunction(InfiniteWavefunctionLeft const& Left_,
                                 WavefunctionSectionLeft const& Window_,
                                 InfiniteWavefunctionRight const& Right_,
                                 int Offset,
                                 int WindowLeft,
                                 int WindowRight)
   : WindowLeftSites(WindowLeft), WindowRightSites(WindowRight), WindowOffset(Offset),
     Left(Left_), Window(Window_), Right(Right_),
     LeftQShift(QuantumNumber(Left.GetSymmetryList())), RightQShift(QuantumNumber(Right.GetSymmetryList()))
{
}

IBCWavefunction::IBCWavefunction(InfiniteWavefunctionLeft const& Left_,
                                 WavefunctionSectionLeft const& Window_,
                                 InfiniteWavefunctionRight const& Right_,
                                 QuantumNumbers::QuantumNumber LeftQShift_,
                                 QuantumNumbers::QuantumNumber RightQShift_,
                                 int Offset,
                                 int WindowLeft,
                                 int WindowRight)
   : WindowLeftSites(WindowLeft), WindowRightSites(WindowRight), WindowOffset(Offset),
     Left(Left_), Window(Window_), Right(Right_),
     LeftQShift(LeftQShift_), RightQShift(RightQShift_)
{
}

void IBCWavefunction::check_structure() const
{
   Left.check_structure();
   Window.check_structure();
   Right.check_structure();
   if (!Left.empty())
   {
      CHECK_EQUAL(delta_shift(Left.Basis2(), LeftQShift), Window.LeftU().Basis1());
   }
   if (!Right.empty())
   {
      CHECK_EQUAL(Window.RightU().Basis2(), delta_shift(Right.Basis1(), RightQShift));
   }
}

void IBCWavefunction::debug_check_structure() const
{
   Left.debug_check_structure();
   Window.debug_check_structure();
   Right.debug_check_structure();
   if (!Left.empty())
   {
      DEBUG_CHECK_EQUAL(delta_shift(Left.Basis2(), LeftQShift), Window.LeftU().Basis1());
   }
   if (!Right.empty())
   {
      DEBUG_CHECK_EQUAL(Window.RightU().Basis2(), delta_shift(Right.Basis1(), RightQShift));
   }
}

PStream::opstream& operator<<(PStream::opstream& out, IBCWavefunction const& Psi)
{
   out << IBCWavefunction::VersionT.default_version();
   out << Psi.WindowLeftSites;
   out << Psi.WindowRightSites;
   out << Psi.WindowOffset;
   out << Psi.WavefunctionLeftFile;
   out << Psi.WavefunctionRightFile;

   if (Psi.WavefunctionLeftFile.empty())
      out << Psi.Left;

   out << Psi.Window;

   if (Psi.WavefunctionRightFile.empty())
      out << Psi.Right;

   out << Psi.LeftQShift;
   out << Psi.RightQShift;

   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, IBCWavefunction& Psi)
{
   int Version = in.read<int>();
   PStream::VersionSentry Sentry(in, IBCWavefunction::VersionT, Version);

   in >> Psi.WindowLeftSites;
   in >> Psi.WindowRightSites;
   in >> Psi.WindowOffset;
   if (Version == 1)
   {
      in >> Psi.Left;
      in >> Psi.Window;
      in >> Psi.Right;
   }
   else if (Version == 2 || Version == 3)
   {
      in >> Psi.WavefunctionLeftFile;
      in >> Psi.WavefunctionRightFile;
      if (Psi.WavefunctionLeftFile.empty())
         in >> Psi.Left;
      else
      {
         pvalue_ptr<MPWavefunction> PsiLeft = pheap::ImportHeap(Psi.WavefunctionLeftFile);
         Psi.Left = PsiLeft->get<InfiniteWavefunctionLeft>();
      }
      in >> Psi.Window;
      if (Psi.WavefunctionRightFile.empty())
         in >> Psi.Right;
      else
      {
         pvalue_ptr<MPWavefunction> PsiRight = pheap::ImportHeap(Psi.WavefunctionRightFile);
         Psi.Right = PsiRight->get<InfiniteWavefunctionRight>();
      }
      if (Version == 3)
      {
         in >> Psi.LeftQShift;
         in >> Psi.RightQShift;
      }
   }
   else
   {
      PANIC("This program is too old to read this wavefunction, expected version <= 3")(Version);
   }

   if (Version < 3)
   {
      if (!Psi.left().empty())
         Psi.LeftQShift = QuantumNumber(Psi.Left.GetSymmetryList());
      if (!Psi.right().empty())
         Psi.RightQShift = QuantumNumber(Psi.Right.GetSymmetryList());
   }

   return in;
}

void
inplace_reflect(IBCWavefunction& Psi)
{
   InfiniteWavefunctionRight Temp = reflect(Psi.Left);
   Psi.Left = reflect(Psi.Right);
   Psi.Right = Temp;
   inplace_reflect(Psi.Window);
   std::swap(Psi.LeftQShift, Psi.RightQShift);
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
   A["RightUnitCellSize"] = Right.size();
   A["LeftFilename"] = this->get_left_filename();
   A["RightFilename"] = this->get_right_filename();
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
         PsiLeft = Psi.left();
         PsiRight = Psi.right();

         inplace_qshift(PsiLeft, Psi.left_qshift());
         inplace_qshift(PsiRight, Psi.right_qshift());

         WindowLeftIndex = Psi.window_offset();
         WindowRightIndex = Psi.window_size() + Psi.window_offset() - 1;

         if (Index < WindowLeftIndex)
         {
            int IndexDiff = WindowLeftIndex - Index;

            for (int i = 0; i < (Psi.window_left_sites() + IndexDiff) / PsiLeft.size(); ++i)
               inplace_qshift(PsiLeft, PsiLeft.qshift());

            C = PsiLeft.end();
            for (int i = 0; i < (Psi.window_left_sites() + IndexDiff) % PsiLeft.size(); ++i)
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

            for (int i = 0; i < (Psi.window_right_sites() + IndexDiff) / PsiRight.size(); ++i)
               inplace_qshift(PsiRight, adjoint(PsiRight.qshift()));

            C = PsiRight.begin();
            for (int i = 0; i < (Psi.window_right_sites() + IndexDiff) % PsiRight.size(); ++i)
               ++C;
         }
         else
         {
            int IndexDiff = Index - WindowLeftIndex;

            C = Psi.window().begin();
            for (int i = 0; i < IndexDiff; ++i)
               ++C;
         }
      }

      StateComponent operator*() const
      {
         StateComponent Result = *C;

         if (Index == WindowRightIndex + 1)
            Result = Psi.window().lambda_r() * Psi.window().RightU() * Result;

         if (Index == WindowLeftIndex)
            Result = Psi.window().LeftU() * Result;

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
            for (int i = 0; i < Psi.window_right_sites(); ++i)
               ++C;
         }
         else if (Index == WindowLeftIndex)
            C = Psi.window().begin();
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
      const IBCWavefunction Psi;
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
   int IndexLeft = std::min(Psi.window_offset() - ((Psi.left().size() - Psi.window_left_sites()) % Op.unit_cell_size()),
                            Op.offset());
   int IndexRight = std::max(Psi.window_size() + Psi.window_offset() + ((Psi.right().size() - Psi.window_right_sites() - 1) % Op.unit_cell_size()),
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

std::tuple<StateComponent, StateComponent>
get_boundary_transfer_eigenvectors(IBCWavefunction const& Psi1, ProductMPO const& StringOp,
                                   IBCWavefunction const& Psi2, double UnityEpsilon, int Verbose)
{
   CHECK_EQUAL(Psi1.left().size(), Psi2.left().size());
   CHECK_EQUAL(Psi1.right().size(), Psi2.right().size());

   int LeftSize = Psi1.left().size();
   int RightSize = Psi1.right().size();

   int IndexLeft1 = Psi1.window_offset() - ((LeftSize - Psi1.window_left_sites()) % LeftSize);
   int IndexLeft2 = Psi2.window_offset() - ((LeftSize - Psi2.window_left_sites()) % LeftSize);
   int IndexLeft = std::min(IndexLeft1, IndexLeft2);

   int IndexRight1 = Psi1.window_size() + Psi1.window_offset() + ((RightSize - Psi1.window_right_sites() - 1) % RightSize);
   int IndexRight2 = Psi2.window_size() + Psi2.window_offset() + ((RightSize - Psi2.window_right_sites() - 1) % RightSize);
   int IndexRight = std::max(IndexRight1, IndexRight2);

   // Ensure that the semi-infinite boundaries have the same quantum numbers.
   InfiniteWavefunctionLeft Psi1Left = Psi1.left();
   inplace_qshift(Psi1Left, Psi1.left_qshift());
   for (int i = 0; i < (IndexLeft1 - IndexLeft) / LeftSize; ++i)
      inplace_qshift(Psi1Left, Psi1Left.qshift());
   if (Psi1.window_left_sites() > 0)
      inplace_qshift(Psi1Left, Psi1Left.qshift());

   InfiniteWavefunctionLeft Psi2Left = Psi2.left();
   inplace_qshift(Psi2Left, Psi2.left_qshift());
   for (int i = 0; i < (IndexLeft2 - IndexLeft) / LeftSize; ++i)
      inplace_qshift(Psi2Left, Psi2Left.qshift());
   if (Psi2.window_left_sites() > 0)
      inplace_qshift(Psi2Left, Psi2Left.qshift());

   // Calculate the left eigenvector of the left semi-infinite boundary.
   std::complex<double> OverlapL;
   MatrixOperator EL, ER;
   std::tie(OverlapL, EL, ER) = get_transfer_eigenpair(get_left_canonical(Psi1Left).first, get_left_canonical(Psi2Left).first, Psi1Left.qshift(), StringOp);

   // Check that the eigenvalue has magnitude 1.
   if (std::abs(std::abs(OverlapL) - 1.0) > UnityEpsilon)
      WARNING("The overlap of the left boundaries is below threshold.")(OverlapL)(std::abs(OverlapL));

   // Normalize OverlapL.
   OverlapL = OverlapL / std::abs(OverlapL);

   // Normalize EL s.t. the sum of the singular values of ER = 1.
   MatrixOperator U, Vh;
   RealDiagonalOperator D;
   SingularValueDecomposition(ER, U, D, Vh);

   ER *= 1.0 / trace(D);
   EL *= 1.0 / inner_prod(delta_shift(ER, Psi1Left.qshift()), EL);

   StateComponent E(StringOp.Basis1(), EL.Basis1(), EL.Basis2());
   E.front() = EL;

   E = delta_shift(E, adjoint(Psi1Left.qshift()));

   InfiniteWavefunctionRight Psi1Right = Psi1.right();
   inplace_qshift(Psi1Right, Psi1.right_qshift());
   for (int i = 0; i < (IndexRight - IndexRight1) / RightSize; ++i)
      inplace_qshift(Psi1Right, adjoint(Psi1Right.qshift()));

   InfiniteWavefunctionRight Psi2Right = Psi2.right();
   inplace_qshift(Psi2Right, Psi2.right_qshift());
   for (int i = 0; i < (IndexRight - IndexRight2) / RightSize; ++i)
      inplace_qshift(Psi2Right, adjoint(Psi2Right.qshift()));

   // Calculate the right eigenvector of the right semi-infinite boundary.
   std::complex<double> OverlapR;
   MatrixOperator FL, FR;
   std::tie(OverlapR, FL, FR) = get_transfer_eigenpair(get_right_canonical(Psi1Right).second, get_right_canonical(Psi2Right).second, Psi1Right.qshift(), StringOp);

   // Check that the eigenvalue has magnitude 1.
   if (std::abs(std::abs(OverlapR) - 1.0) > UnityEpsilon)
      WARNING("The overlap of the right boundaries is below threshold.")(OverlapR)(std::abs(OverlapR));

   // Normalize OverlapR.
   OverlapR = OverlapR / std::abs(OverlapR);

   // Normalize FR s.t. the sum of the singular values of FL = 1.
   SingularValueDecomposition(FL, U, D, Vh);

   FL *= 1.0 / trace(D);
   FR *= 1.0 / inner_prod(delta_shift(FR, Psi1Right.qshift()), FL);

   StateComponent F(StringOp.Basis1(), FR.Basis1(), FR.Basis2());
   F.front() = FR;

   // Remove spurious phase from E and F by setting the phase of the trace
   // to be zero (but only if the trace is nonzero).
   // This will only work if E and F are square.
   // This only ensures that any overlaps calculated using the same E and F
   // matrices will have the same global phase: calculations using different
   // boundary transfer matrices may end up having a nontrivial global phase
   // shift.
   // TODO: Figure out a better method to fix the phase.
   if (E.Basis1() == E.Basis2() && F.Basis1() == F.Basis2() && E.size() == 1 && F.size() == 1)
   {
      std::complex<double> ETrace = trace(E.front());

      if (std::abs(ETrace) > TraceTol)
         E *= std::conj(ETrace) / std::abs(ETrace);
      else
         WARNING("The trace of E is below threshold, so the overlap will have a spurious phase contribution.")(ETrace);

      std::complex<double> FTrace = trace(F.front());

      if (std::abs(FTrace) > TraceTol)
         F *= std::conj(FTrace) / std::abs(FTrace);
      else
         WARNING("The trace of F is below threshold, so the overlap will have a spurious phase contribution.")(FTrace);
   }
   else
      WARNING("Psi1 and Psi2 have different boundary bases, so the overlap will have a spurious phase contribution.");

   // Compensate for the boundary contributions to the phase for the boundaries
   // which have been incorporated into the window.
   E *= std::pow(OverlapL, IndexLeft / LeftSize);
   F *= std::pow(OverlapR, IndexRight / RightSize);

   return std::make_tuple(E, F);
}

std::complex<double>
overlap(IBCWavefunction const& Psi1, ProductMPO const& StringOp, IBCWavefunction const& Psi2, double UnityEpsilon, int Verbose)
{
   StateComponent E, F;
   std::tie(E, F) = get_boundary_transfer_eigenvectors(Psi1, StringOp, Psi2, UnityEpsilon, Verbose);
   return overlap(Psi1, StringOp, Psi2, E, F, Verbose);
}

std::complex<double>
overlap(IBCWavefunction const& Psi1, IBCWavefunction const& Psi2, double UnityEpsilon, int Verbose)
{
   return overlap(Psi1, ProductMPO::make_identity(ExtractLocalBasis(Psi2.left())), Psi2, UnityEpsilon, Verbose);
}

std::complex<double>
overlap(IBCWavefunction const& Psi1, ProductMPO const& StringOp, IBCWavefunction const& Psi2,
        StateComponent const& E_, StateComponent const& F_, int Verbose)
{
   CHECK_EQUAL(Psi1.left().size(), Psi2.left().size());
   CHECK_EQUAL(Psi1.right().size(), Psi2.right().size());

   int LeftSize = Psi1.left().size();
   int RightSize = Psi1.right().size();

   int IndexLeft1 = Psi1.window_offset() - ((LeftSize - Psi1.window_left_sites()) % LeftSize);
   int IndexLeft2 = Psi2.window_offset() - ((LeftSize - Psi2.window_left_sites()) % LeftSize);
   int IndexLeft = std::min(IndexLeft1, IndexLeft2);

   int IndexRight1 = Psi1.window_size() + Psi1.window_offset() + ((RightSize - Psi1.window_right_sites() - 1) % RightSize);
   int IndexRight2 = Psi2.window_size() + Psi2.window_offset() + ((RightSize - Psi2.window_right_sites() - 1) % RightSize);
   int IndexRight = std::max(IndexRight1, IndexRight2);

   if (Verbose > 0)
      std::cout << "Calculating IBC overlap over sites " << IndexLeft << " to " << IndexRight
                << " (" << IndexRight - IndexLeft + 1 << " sites total)" << std::endl;

   StateComponent E = E_;
   StateComponent F = F_;

   // Calculate the overlap.
   ConstIBCIterator C1 = ConstIBCIterator(Psi1, IndexLeft);
   ConstIBCIterator C2 = ConstIBCIterator(Psi2, IndexLeft);
   ProductMPO::const_iterator W = StringOp.begin();

   // Initialize the position of W.
   for (int i = 0; i < (IndexLeft % LeftSize + LeftSize) % LeftSize; ++i)
   {
      ++W;
      if (W == StringOp.end())
         W = StringOp.begin();
   }

   for (int i = IndexLeft; i <= IndexRight; ++i)
   {
      if (Verbose > 2)
         std::cout << "Site " << i << std::endl;

      E = contract_from_left(*W, herm(*C1), E, *C2);

      ++C1, ++C2, ++W;
      if (W == StringOp.end())
         W = StringOp.begin();
   }

   return inner_prod(F, E);
}

std::complex<double>
overlap(IBCWavefunction const& Psi1, IBCWavefunction const& Psi2,
        StateComponent const& E, StateComponent const& F, int Verbose)
{
   return overlap(Psi1, ProductMPO::make_identity(ExtractLocalBasis(Psi2.left())), Psi2, E, F, Verbose);
}
