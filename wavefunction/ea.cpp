// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// wavefunction/ea.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "ea.h"
#include "mpwavefunction.h"

//
// EAWavefunction
//

// Streaming versions:
//
// Version 1:
// InfiniteWavefunctionLeft Left
// std::vector<WavefunctionSectionLeft> WindowVec
// InfiniteWavefunctionRight Right
// std::complex<double> ExpIK
// std::complex<double> GSOverlap
//
// Version 2:
// string WavefunctionLeftFile
// string WavefunctionRightFile
// InfiniteWavefunctionLeft Left (only if WavefunctionLeftFile is empty)
// std::vector<WavefunctionSectionLeft> WindowVec
// InfiniteWavefunctionRight Right (only if WavefunctionRightFile is empty)
// std::complex<double> ExpIK
// std::complex<double> GSOverlap
//
// Version 3:
// string WavefunctionLeftFile
// string WavefunctionRightFile
// InfiniteWavefunctionLeft Left (only if WavefunctionLeftFile is empty)
// std::vector<WavefunctionSectionLeft> WindowVec
// InfiniteWavefunctionRight Right (only if WavefunctionRightFile is empty)
// int LeftIndex
// int RightIndex
// QuantumNumber LeftQShift
// QuantumNumber RightQShift
// std::complex<double> ExpIK
// std::complex<double> GSOverlap

PStream::VersionTag
EAWavefunction::VersionT(3);

std::string const EAWavefunction::Type = "EAWavefunction";

EAWavefunction::EAWavefunction()
   : LeftIndex(0), RightIndex(0), ExpIK(1.0), GSOverlap(0.0)
{
}

EAWavefunction::EAWavefunction(InfiniteWavefunctionLeft const& Left_,
                               std::vector<WavefunctionSectionLeft> const& WindowVec_,
                               InfiniteWavefunctionRight const& Right_,
                               int LeftIndex_,
                               int RightIndex_,
                               std::complex<double> ExpIK_,
                               std::complex<double> GSOverlap_)
   : Left(Left_), WindowVec(WindowVec_), Right(Right_),
     LeftIndex(LeftIndex_), RightIndex(RightIndex_),
     ExpIK(ExpIK_), GSOverlap(GSOverlap_),
     LeftQShift(QuantumNumber(Left.GetSymmetryList())), RightQShift(QuantumNumber(Right.GetSymmetryList()))
{
}

EAWavefunction::EAWavefunction(InfiniteWavefunctionLeft const& Left_,
                               std::vector<WavefunctionSectionLeft> const& WindowVec_,
                               InfiniteWavefunctionRight const& Right_,
                               QuantumNumber LeftQShift_,
                               QuantumNumber RightQShift_,
                               int LeftIndex_,
                               int RightIndex_,
                               std::complex<double> ExpIK_,
                               std::complex<double> GSOverlap_)
   : Left(Left_), WindowVec(WindowVec_), Right(Right_),
     LeftIndex(LeftIndex_), RightIndex(RightIndex_),
     ExpIK(ExpIK_), GSOverlap(GSOverlap_),
     LeftQShift(LeftQShift_), RightQShift(RightQShift_)
{
}

void EAWavefunction::check_structure() const
{
   Left.check_structure();
   for (auto Window : WindowVec)
      Window.check_structure();
   Right.check_structure();
   // TODO
#if 0
   if (!Left.empty())
   {
      CHECK_EQUAL(Left.Basis2(), Window.Basis1());
   }
   if (!Right.empty())
   {
      CHECK_EQUAL(Window.Basis2(), Right.Basis1());
   }
#endif
   CHECK_EQUAL(Left.size(), Right.size());
   CHECK_EQUAL(Left.size(), WindowVec.size());
   CHECK(LeftIndex >= 0)(LeftIndex < Left.size());
   CHECK(RightIndex >= 0)(RightIndex < Right.size());
}

void EAWavefunction::debug_check_structure() const
{
   Left.debug_check_structure();
   for (auto Window : WindowVec)
      Window.debug_check_structure();
   Right.debug_check_structure();
   // TODO
#if 0
   if (!Left.empty())
   {
      DEBUG_CHECK_EQUAL(Left.Basis2(), Window.Basis1());
   }
   if (!Right.empty())
   {
      DEBUG_CHECK_EQUAL(Window.Basis2(), Right.Basis1());
   }
#endif
   DEBUG_CHECK_EQUAL(Left.size(), Right.size());
   DEBUG_CHECK_EQUAL(Left.size(), WindowVec.size());
   DEBUG_CHECK(LeftIndex >= 0)(LeftIndex < Left.size());
   DEBUG_CHECK(RightIndex >= 0)(RightIndex < Right.size());
}

PStream::opstream& operator<<(PStream::opstream& out, EAWavefunction const& Psi)
{
   out << EAWavefunction::VersionT.default_version();

   out << Psi.WavefunctionLeftFile;
   out << Psi.WavefunctionRightFile;

   if (Psi.WavefunctionLeftFile.empty())
      out << Psi.Left;

   out << Psi.WindowVec;

   if (Psi.WavefunctionRightFile.empty())
      out << Psi.Right;

   out << Psi.LeftIndex;
   out << Psi.RightIndex;

   out << Psi.LeftQShift;
   out << Psi.RightQShift;

   out << Psi.ExpIK;
   out << Psi.GSOverlap;

   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, EAWavefunction& Psi)
{
   int Version = in.read<int>();
   PStream::VersionSentry Sentry(in, EAWavefunction::VersionT, Version);

   if (Version == 1)
   {
      in >> Psi.Left;
      in >> Psi.WindowVec;
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
      in >> Psi.WindowVec;
      if (Psi.WavefunctionRightFile.empty())
         in >> Psi.Right;
      else
      {
         pvalue_ptr<MPWavefunction> PsiRight = pheap::ImportHeap(Psi.WavefunctionRightFile);
         Psi.Right = PsiRight->get<InfiniteWavefunctionRight>();
      }
      if (Version == 3)
      {
         in >> Psi.LeftIndex;
         in >> Psi.RightIndex;

         in >> Psi.LeftQShift;
         in >> Psi.RightQShift;
      }
   }
   else
   {
      PANIC("This program is too old to read this wavefunction, expected version <= 3")(Version);
   }

   in >> Psi.ExpIK;
   in >> Psi.GSOverlap;

   if (Version < 3)
   {
      Psi.LeftIndex = 0;
      Psi.RightIndex = 0;

      if (!Psi.left().empty())
         Psi.LeftQShift = QuantumNumber(Psi.Left.GetSymmetryList());
      if (!Psi.right().empty())
         Psi.RightQShift = QuantumNumber(Psi.Right.GetSymmetryList());
   }

   return in;
}

// TODO: Test that this works properly.
void
inplace_reflect(EAWavefunction& Psi)
{
   InfiniteWavefunctionRight Temp = reflect(Psi.Left);
   Psi.Left = reflect(Psi.Right);
   Psi.Right = Temp;
   for (auto& Window : Psi.WindowVec)
      inplace_reflect(Window);
   std::swap(Psi.LeftQShift, Psi.RightQShift);
}

void
inplace_conj(EAWavefunction& Psi)
{
   inplace_conj(Psi.Left);
   for (auto& Window : Psi.WindowVec)
      inplace_conj(Window);
   inplace_conj(Psi.Right);
}

void
EAWavefunction::SetDefaultAttributes(AttributeList& A) const
{
   A["WavefunctionType"] = "EA";
   A["WindowSize"] = this->window_size();
   A["ExpIK"] = ExpIK;
   A["GSOverlap"] = GSOverlap;
   A["LeftUnitCellSize"] = Left.size();
   A["RightUnitCellSize"] = Right.size();
   A["LeftFilename"] = this->get_left_filename();
   A["RightFilename"] = this->get_right_filename();
}
