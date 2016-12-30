// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/mpwavefunction.cpp
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

#include "mpwavefunction.h"

// default version is 5.
//
// Versions 1 and 2 don't contain an AttributeList (although older versions of
// InfiniteWavefunction do, but it isn't used anywhere so no need to keep it).
// Version 4 introduces the HistoryLog.
//
// Version 1:
// No separate version number, this needs to be detected from the stream metadata version.
// In this case, the stream contains a InfiniteWavefunctionLeft in version 1 format
// (no separate version number)
//
// Version 2:
// InfiniteWavefunctionLeft in version 2 format (no separate version number)
//
// Version 3:
// boost::variant<InfiniteWavefunctionLeft>
// AttributeList
//
// Version 4:
// boost::variant<InfiniteWavefunctionLeft>
// AttributeList
// HistoryLog
//
// Version 5:
// boost::variant<InfiniteWavefunctionLeft, IBCWavefunction>
// AttributeList
// HistoryLog


PStream::VersionTag MPWavefunction::VersionT(5);

PStream::ipstream&
operator>>(PStream::ipstream& in, MPWavefunction& Psi)
{
   // old streams didn't have a version number, so hack around it with the metadata version

   int Version;

   PHeapFileSystem::ipheapstream* S = dynamic_cast<PHeapFileSystem::ipheapstream*>(&in);
   if (S && S->version() == 1)
   {
      // old file, need to set the version number manually
      Version = 1;
   }
   else
   {
      // new file, read the version number
      Version = in.read<int>();
   }

   PStream::VersionSentry Sentry(in, MPWavefunction::VersionT, Version);
   read_version(in, Psi, Sentry.version());
   return in;
}

void read_version(PStream::ipstream& in, MPWavefunction& Psi, int Version)
{
   DEBUG_TRACE("Reading MPWavefunction")(Version);
   if (Version == 1)
   {
      InfiniteWavefunctionLeft x;
      read_version(in, x, 1);
      Psi = x;
      return;
   }

   if (Version == 2)
   {
      InfiniteWavefunctionLeft x;
      read_version(in, x, 2);
      Psi = x;
      return;
   }

   if (Version <= 4)
   {
      boost::variant<InfiniteWavefunctionLeft> x;
      in >> x;
      Psi.Psi_ = x;
   }
   else if (Version == 5)
   {
      boost::variant<InfiniteWavefunctionLeft, IBCWavefunction> x;
      in >> x;
      Psi.Psi_ = x;
   }

   in >> Psi.Attr_;

   if (Version >= 4)
      in >> Psi.History_;

   if (Version > 5)
   {
      PANIC("Version of MPWavefunction is newer than this sofware!");
   }
}

PStream::opstream&
operator<<(PStream::opstream& out, MPWavefunction const& Psi)
{
   out << MPWavefunction::VersionT.default_version();
   out << Psi.Psi_;
   out << Psi.Attr_;
   out << Psi.History_;
   return out;
}

struct ApplyCheckStructure : boost::static_visitor<void>
{
   template <typename T>
   void operator()(T const& x) const
   {
      x.check_structure();
   }
};

void
MPWavefunction::check_structure() const
{
   boost::apply_visitor(ApplyCheckStructure(), this->Wavefunction());
}

void
MPWavefunction::debug_check_structure() const
{
#if !defined(NDEBUG)
   boost::apply_visitor(ApplyCheckStructure(), this->Wavefunction());
#endif
}

struct DoSetDefaultAttributes : public boost::static_visitor<void>
{
   DoSetDefaultAttributes(AttributeList& A_) : A(A_) {}

   template <typename T>
   void operator()(T const& Psi) const
   {
      Psi.SetDefaultAttributes(A);
   }

   AttributeList& A;
};

void
MPWavefunction::SetDefaultAttributes()
{
   boost::apply_visitor(DoSetDefaultAttributes(Attr_), Psi_);
}
