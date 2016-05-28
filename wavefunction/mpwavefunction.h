// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/mpwavefunction.h
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
// MPWavefunction - the main wavefunction class for the MPToolkit
// The MPWavefunction is a holder that contains a variant for the concrete wavefunction types.
// The MPWavefunction also holds the wavefunction attributes.

#if !defined(MPTOOLKIT_WAVEFUNCTION_MPWAVEFUNCTION_H)
#define MPTOOLKIT_WAVEFUNCTION_MPWAVEFUNCTION_H

#include "infinitewavefunctionleft.h"
#include "ibc.h"
#include <boost/variant.hpp>
#include "pstream/pstream.h"
#include "pstream/variant.h"
#include "interface/attributes.h"
#include "interface/history.h"

typedef boost::variant<InfiniteWavefunctionLeft, IBCWavefunction> WavefunctionTypes;

class MPWavefunction
{
   public:
      MPWavefunction() {}

      MPWavefunction(WavefunctionTypes const& Psi) : Psi_(Psi) {}

      MPWavefunction(WavefunctionTypes const& Psi, AttributeList const& Attr) 
	 : Psi_(Psi), Attr_(Attr) {}

      // construct an MPWavefunction with the given attributes and history
      MPWavefunction(AttributeList const& Attr, HistoryLog const& Hist)
	 : Attr_(Attr), History_(Hist) {}

      MPWavefunction& operator=(MPWavefunction const& Psi2)
      { Psi_ = Psi2.Psi_; Attr_ = Psi2.Attr_; History_ = Psi2.History_; return *this; }

      MPWavefunction& operator=(WavefunctionTypes const& Psi2)
      { Psi_ = Psi2; return *this; }

      WavefunctionTypes& Wavefunction() { return Psi_; }

      WavefunctionTypes const& Wavefunction() const { return Psi_; }

      // Get the wavefunction, assuming it is of the specified type (will fail otherwise)
      template <typename T>
      T& get() { return boost::get<T>(Psi_); }

      template <typename T>
      T const& get() const { return boost::get<T>(Psi_); }

      // returns true if the MPWavefunction contains the given type
      template <typename T>
      bool is() const { return bool(boost::get<T>(&Psi_)); }

      AttributeList& Attributes() { return Attr_; }
      AttributeList const& Attributes() const { return Attr_; }

      // Set the default attributes for the wavefunction type
      void SetDefaultAttributes();

      // returns the history log
      HistoryLog const& History() const { return History_; }

      // appends an entry to the history log
      void AppendHistory(std::string const& s)
      { History_.append(s); }

      static PStream::VersionTag VersionT;

      void check_structure() const;

      void debug_check_structure() const;

   private:
      WavefunctionTypes Psi_;
      AttributeList Attr_;
      HistoryLog History_;

      friend PStream::ipstream& operator>>(PStream::ipstream& in, MPWavefunction& Psi);
      friend PStream::opstream& operator<<(PStream::opstream& out, MPWavefunction const& Psi);
      friend void read_version(PStream::ipstream& in, MPWavefunction& Psi, int Version);
};

#endif
