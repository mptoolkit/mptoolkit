// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/mpwavefunction.h
//
// Copyright (C) 2015-2020 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "infinitewavefunctionright.h"
#include "ibc.h"
#include "finitewavefunctionleft.h"
#include <boost/variant.hpp>
#include "pstream/pstream.h"
#include "pstream/variant.h"
#include "interface/attributes.h"
#include "interface/history.h"

typedef boost::variant<InfiniteWavefunctionLeft,
                       IBCWavefunction,
                       FiniteWavefunctionLeft,
                       InfiniteWavefunctionRight> WavefunctionTypes;

class InvalidWavefunction : public std::runtime_error
{
   public:
      explicit InvalidWavefunction(std::string const& Expected, std::string const& Suplied);
};

inline
InvalidWavefunction::InvalidWavefunction(std::string const& Expected, std::string const& Supplied)
   : std::runtime_error("Wavefunction type is not correct, expected " + Expected + " but was given "
                        + Supplied)
{
}

class MPWavefunction
{
   public:
      MPWavefunction() : Version_(0) {}

      MPWavefunction(WavefunctionTypes const& Psi) : Psi_(Psi), Version_(0) {}

      MPWavefunction(WavefunctionTypes const& Psi, AttributeList const& Attr)
         : Psi_(Psi), Attr_(Attr), Version_(0) {}

      // construct an MPWavefunction with the given attributes and history
      MPWavefunction(AttributeList const& Attr, HistoryLog const& Hist)
         : Attr_(Attr), History_(Hist), Version_(0) {}

      MPWavefunction& operator=(MPWavefunction const& Psi2)
      { Psi_ = Psi2.Psi_; Attr_ = Psi2.Attr_; History_ = Psi2.History_; Version_ = 0; return *this; }

      MPWavefunction& operator=(WavefunctionTypes const& Psi2)
      { Psi_ = Psi2; Version_ = 0; return *this; }

      WavefunctionTypes& Wavefunction() { return Psi_; }

      WavefunctionTypes const& Wavefunction() const { return Psi_; }

      // Get the wavefunction, assuming it is of the specified type (will fail otherwise)
      template <typename T>
      T& get();

      template <typename T>
      T const& get() const;

      // returns true if the MPWavefunction contains the given type
      template <typename T>
      bool is() const { return bool(boost::get<T>(&Psi_)); }

      // returns a string indicating the type of wavefunction
      std::string Type() const;

      AttributeList& Attributes() { return Attr_; }
      AttributeList const& Attributes() const { return Attr_; }

      // Set the default attributes for the wavefunction type
      void SetDefaultAttributes();

      // returns the history log
      HistoryLog const& History() const { return History_; }

      // appends an entry to the history log
      void AppendHistoryCommand(std::string const& s)
      { History_.append_command(s); }

      void AppendHistoryNote(std::string const& s)
      { History_.append_note(s); }

      static PStream::VersionTag VersionT;

      void check_structure() const;

      void debug_check_structure() const;

      // read-only; version number read from the stream.  Will be zero for
      // wavefunctions that are not read directly from a stream.
      int version() const { return Version_; }

   private:
      WavefunctionTypes Psi_;
      AttributeList Attr_;
      HistoryLog History_;
      int Version_;

      friend PStream::ipstream& operator>>(PStream::ipstream& in, MPWavefunction& Psi);
      friend PStream::opstream& operator<<(PStream::opstream& out, MPWavefunction const& Psi);
      friend void read_version(PStream::ipstream& in, MPWavefunction& Psi, int Version);
};

template <typename T>
T& MPWavefunction::get()
{
   T* Result = boost::get<T>(&Psi_);
   if (!Result)
   {
      throw InvalidWavefunction(T::Type, this->Type());
   }
   return *Result;
}

template <typename T>
T const& MPWavefunction::get() const
{
   T const* Result = boost::get<T>(&Psi_);
   if (!Result)
   {
      throw InvalidWavefunction(T::Type, this->Type());
   }
   return *Result;
}

#endif
