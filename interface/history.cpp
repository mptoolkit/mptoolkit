// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// interface/history.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "history.h"
#include <iomanip>

std::string
HistoryEntry::MakeTimestamp(time_t Timestamp_)
{
   char s[200];
   int const max = 200;
   if (strftime(s, max, HISTORY_TIMESTAMP_FORMAT, localtime(&Timestamp_)) == 0)
      s[0] = '\0';
   return s;
}

HistoryEntry::HistoryEntry(std::string const& Str)
   : Entry(Str)
{
   Timestamp = MakeTimestamp(time(NULL));
}

HistoryEntry::HistoryEntry(std::string const& Str, time_t Timestamp_)
   : Entry(Str)
{
   Timestamp = MakeTimestamp(Timestamp_);
}

std::ostream& operator<<(std::ostream& out, HistoryEntry const& H)
{
   out << H.formatted();
   return out;
}

#if defined(USE_PSTREAM)
PStream::opstream& operator<<(PStream::opstream& out, HistoryEntry const& H)
{
   out << H.Timestamp;
   out << H.Entry;
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, HistoryEntry& H)
{
   in >> H.Timestamp;
   in >> H.Entry;
   return in;
}
#endif

//
// HistoryLog
//

void
HistoryLog::append(HistoryEntry const& H)
{
   History_.push_back(H);
}

void
HistoryLog::append(std::string const& Entry)
{
   History_.push_back(HistoryEntry(Entry));
}

HistoryEntry
HistoryLog::front() const
{ 
   if (History_.empty())
      return HistoryEntry("(empty history log)");
   return History_.front(); 
}

HistoryEntry
HistoryLog::back() const
{ 
   if (History_.empty())
      return HistoryEntry("(empty history log)");
   return History_.back(); 
}

void
HistoryLog::print(std::ostream& out) const
{
   for (const_iterator I = this->begin(); I != this->end(); ++I)
   {
      out << (*I) << "\n\n";
   }
   out << std::flush;
}

void
HistoryLog::print_newest_first(std::ostream& out) const
{
   for (const_iterator I = this->end(); I != this->begin(); --I)
   {
      const_iterator J = I;
      --J;
      out << (*J) << "\n\n";
   }
}

#if defined(USE_PSTREAM)
PStream::opstream& operator<<(PStream::opstream& out, HistoryLog const& Log)
{
   out << Log.History_;
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, HistoryLog& Log)
{
   in >> Log.History_;
   return in;
}

#endif
