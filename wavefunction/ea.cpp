// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/ea.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "ea.h"

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

PStream::VersionTag
EAWavefunction::VersionT(1);

std::string const EAWavefunction::Type = "EAWavefunction";

EAWavefunction::EAWavefunction()
   : ExpIK(1.0), GSOverlap(0.0)
{
}

EAWavefunction::EAWavefunction(InfiniteWavefunctionLeft const& Left_,
                               std::vector<WavefunctionSectionLeft> const& WindowVec_,
                               InfiniteWavefunctionRight const& Right_,
                               std::complex<double> ExpIK_,
                               std::complex<double> GSOverlap_)
   : Left(Left_), WindowVec(WindowVec_), Right(Right_), ExpIK(ExpIK_), GSOverlap(GSOverlap_)
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
}

PStream::opstream& operator<<(PStream::opstream& out, EAWavefunction const& Psi)
{
   out << EAWavefunction::VersionT.default_version();
   out << Psi.Left;
   out << Psi.WindowVec;
   out << Psi.Right;
   out << Psi.ExpIK;
   out << Psi.GSOverlap;
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, EAWavefunction& Psi)
{
   int Version = in.read<int>();
   PStream::VersionSentry Sentry(in, EAWavefunction::VersionT, Version);

   if (Version != 1)
   {
      PANIC("This program is too old to read this wavefunction, expected version = 1")(Version);
   }

   in >> Psi.Left;
   in >> Psi.WindowVec;
   in >> Psi.Right;
   in >> Psi.ExpIK;
   in >> Psi.GSOverlap;

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
   A["RightUnitCellSize"] = Left.size();
}
