// -*- C++ -*-

#include "mpwavefunction.h"

// default version is 3.
//
// Versions 1 and 2 don't contain an AttributeList (although older versions of
// InfiniteWavefunction do, but it isn't used anywhere so no need to keep it)
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

PStream::VersionTag MPWavefunction::VersionT(3);

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
   }
   else if (Version == 2)
   {
      InfiniteWavefunctionLeft x;
      read_version(in, x, 2);
      Psi = x;
   }
   else if (Version == 3)
   {
      in >> Psi.Psi_;
      in >> Psi.Attr_;
   }
   else
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
   return out;
}
