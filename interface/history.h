// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// interface/history.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
//
// history.h
//
// A simple history log for the MPToolkit

#if !defined(MPTOOLKIT_INTERFACE_HISTORY_H)
#define MPTOOLKIT_INTERFACE_HISTORY_H

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <list>
#include <string>
#include <time.h>
#include <iostream>
#if defined(USE_PSTREAM)
#include "pstream/pstream.h"
#endif

#define HISTORY_TIMESTAMP_FORMAT "%a, %d %b %Y %T %z"

class HistoryEntry
{
   public:
      HistoryEntry() {}

      // compiler generated ctors are OK

      explicit HistoryEntry(std::string const& Str);

      HistoryEntry(std::string const& Str, time_t Tatestamp_);

      std::string timestamp() const { return Timestamp; }

      std::string entry() const { return Entry; }

      std::string formatted() const { return "Date: " + Timestamp + '\n' + Entry; }

      // Helper to construct a timestamp from a time_t, using the current local timezone
      static std::string MakeTimestamp(time_t Timestamp_);

   private:
      std::string Timestamp;
      std::string Entry;

#if defined(USE_PSTREAM)
      friend PStream::opstream& operator<<(PStream::opstream& out, HistoryEntry const& H);
      friend PStream::ipstream& operator>>(PStream::ipstream& in, HistoryEntry& H);
#endif
};

std::ostream& operator<<(std::ostream& out, HistoryEntry const& H);


class HistoryLog
{
   public:
      // compiler generated ctors are OK

      typedef std::list<HistoryEntry>::const_iterator const_iterator;

      // adds a new entry (assumed to be later in time than previous entries)
      void append(HistoryEntry const& H);

      // adds a new entry, with timestamp set to the current system time
      void append(std::string const& Entry);

      const_iterator begin() const { return History_.begin(); }
      const_iterator end() const { return History_.end(); }

      HistoryEntry front() const;
      HistoryEntry back() const;

      // prints the history log to the specified stream, oldest entry first
      void print(std::ostream& out) const;

      // prints the history log to the specified stream, newest entry first 
      void print_newest_first(std::ostream& out) const;

   private:
      std::list<HistoryEntry> History_;

#if defined(USE_PSTREAM)
      friend PStream::opstream& operator<<(PStream::opstream& out, HistoryLog const& Log);
      friend PStream::ipstream& operator>>(PStream::ipstream& in, HistoryLog& Log);
#endif
      
};

#endif
