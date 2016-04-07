// -*- C++ -*-

#include "ibc.h"

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

WavefunctionSectionLeft::WavefunctionSectionLeft(InfiniteWavefunctionLeft const& Psi)
   : CanonicalWavefunctionBase(Psi)
{
   Psi.check_structure();
   this->check_structure();
}

PStream::opstream& operator<<(PStream::opstream& out, WavefunctionSectionLeft const& Psi)
{
   out << WavefunctionSectionLeft::VersionT.default_version();
   Psi.CanonicalWavefunctionBase::WriteStream(out);
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
